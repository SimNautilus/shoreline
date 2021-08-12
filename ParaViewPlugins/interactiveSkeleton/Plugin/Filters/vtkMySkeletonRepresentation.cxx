/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMySkeletonRepresentation.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMySkeletonRepresentation.h"

#include "vtkActor.h"
#include "vtkAssemblyPath.h"
#include "vtkBoundingBox.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellPicker.h"
#include "vtkConeSource.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkInteractorObserver.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPickingManager.h"
#include "vtkPlaneSource.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkDoubleArray.h"
#include "vtkVector.h"
#include "vtkIdList.h"

#include "vtkMySkeletonSource.h"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <cmath>

//------------------------------------------------------------------------------
vtkStandardNewMacro(vtkMySkeletonRepresentation);
//------------------------------------------------------------------------------
vtkMySkeletonRepresentation::vtkMySkeletonRepresentation()
{
  //////////////////////////////////////////////////////////////////////////////
  // Read the representation of the widget from file
  vtkNew<vtkXMLPolyDataReader> reader;
  std::string cache_path_skel = std::string(CACHE_PATH) + std::string("/skeleton.vtp");
  reader->SetFileName( cache_path_skel.c_str() );
  reader->Update();

  this->Skeleton = vtkMySkeletonSource::New();
  this->Skeleton->SetSkeleton(reader->GetOutput());
  this->Skeleton->Modified();
  this->Skeleton->Update();

  vtkNew<vtkXMLPolyDataWriter> writer;
  std::string cache_path_skelMod = std::string(CACHE_PATH) + std::string("/skeletonMod.vtp");
  writer->SetFileName ( cache_path_skelMod.c_str() );
  writer->SetInputData ( this->Skeleton->GetSkeleton() );
  writer->Write();

  //////////////////////////////////////////////////////////////////////////////
  // initialize the Eigen components that compute the bendy behavior
  vtkPolyData* sklWidget = this->Skeleton->GetSkeleton();
  const vtkIdType nVerts = sklWidget->GetNumberOfPoints();
  const vtkIdType nEdges = sklWidget->GetNumberOfLines();
  this->b_vertexDegree.resize(nVerts);
  this->b_fixedVerts.resize(nVerts);
  this->b_bcVerts.resize(nVerts);

  // build the neighbor map and the vertex degree map
  for(int i = 0; i < nEdges; ++i) {
    vtkSmartPointer<vtkIdList> pts = sklWidget->GetCell(i)->GetPointIds();
    vtkIdType nPts = sklWidget->GetCell(i)->GetNumberOfPoints();

    std::list<int> neighbors;

    for(int p = 0; p < nPts; ++p) {
      int vInd = pts->GetId(p);

      // get neighbors about point vInd
      for(int n = 0; n < nPts; ++n) {
        if(n == p) continue;
        neighbors.push_back(pts->GetId(n));
      }

      // check if the handle vInd is in the map, if not add it
      if(this->b_neighborsMap.count(vInd)) {
        this->b_neighborsMap[vInd].merge(neighbors);
        this->b_neighborsMap[vInd].sort();
        this->b_neighborsMap[vInd].unique();
      } else {
        this->b_neighborsMap[vInd] = neighbors;
      }

      this->b_vertexDegree[vInd] += 1;
    }
  }

  // collect the vertices to fix and mark them as bcs
  for(int i = 0; i < nVerts; ++i) {
    const int degree = this->b_vertexDegree[i];
    if(degree != 2) {
      this->b_fixedVerts[i] = 1; // identify a fixed vertex
      this->b_bcVerts[i] = 1; // identify a boundary condition
    }
  }

  // initialize the eigen data
  this->b_matL = Eigen::MatrixXd::Zero(nVerts, nVerts);
  this->b_matM = Eigen::MatrixXd::Zero(nVerts, nVerts);

  this->b_vecb = Eigen::MatrixXd::Zero(nVerts,3);
  this->b_vecx = Eigen::MatrixXd::Zero(nVerts,3);
  this->b_arcLength = Eigen::VectorXd::Zero(nVerts);

  //////////////////////////////////////////////////////////////////////////////
  this->LastEventPosition[0] = VTK_DOUBLE_MAX;
  this->LastEventPosition[1] = VTK_DOUBLE_MAX;
  this->LastEventPosition[2] = VTK_DOUBLE_MAX;

  this->Bounds[0] = VTK_DOUBLE_MAX;
  this->Bounds[1] = -VTK_DOUBLE_MAX;
  this->Bounds[2] = VTK_DOUBLE_MAX;
  this->Bounds[3] = -VTK_DOUBLE_MAX;
  this->Bounds[4] = VTK_DOUBLE_MAX;
  this->Bounds[5] = -VTK_DOUBLE_MAX;

  this->HandleSize = 5.0;

  this->InteractionState = vtkMySkeletonRepresentation::Outside;
  this->ProjectToPlane = 0;   // default off
  this->ProjectionNormal = 0; // default YZ not used
  this->ProjectionPosition = 0.0;
  this->PlaneSource = nullptr;
  this->Closed = 0;

  // Build the representation of the widget

  this->DirectionalLine = false;

  // Create the handles along a straight line within the bounds of a unit cube
  this->NumberOfHandles = this->Skeleton->GetNumberOfPoints();
  this->Handle = new vtkActor*[this->NumberOfHandles];
  this->HandleGeometry = new HandleSource*[this->NumberOfHandles];

  double pt[3];
  vtkPoints* points = this->Skeleton->GetPoints();
  for (int i = 0; i < this->NumberOfHandles; ++i)
  {
    this->HandleGeometry[i] = HandleSource::New();
    vtkPolyDataMapper* handleMapper = vtkPolyDataMapper::New();
    handleMapper->SetInputConnection(this->HandleGeometry[i]->GetOutputPort());
    this->Handle[i] = vtkActor::New();
    this->Handle[i]->SetMapper(handleMapper);
    handleMapper->Delete();
    points->GetPoint(i, pt);
    this->HandleGeometry[i]->SetCenter(pt);
  }

  this->LineActor = vtkActor::New();

  // Default bounds to get started
  double bounds[6] = { -0.5, 0.5, -0.5, 0.5, -0.5, 0.5 };

  // Initial creation of the widget, serves to initialize it
  this->PlaceFactor = 1.0;
  this->PlaceWidget(bounds);

  // Manage the picking stuff
  this->HandlePicker = vtkCellPicker::New();
  this->HandlePicker->SetTolerance(0.005);

  for (int i = 0; i < this->NumberOfHandles; ++i)
  {
    this->HandlePicker->AddPickList(this->Handle[i]);
  }
  this->HandlePicker->PickFromListOn();

  this->LinePicker = vtkCellPicker::New();
  this->LinePicker->SetTolerance(0.01);
  this->LinePicker->AddPickList(this->LineActor);
  this->LinePicker->PickFromListOn();

  this->LastPickPosition[0] = VTK_DOUBLE_MAX;
  this->LastPickPosition[1] = VTK_DOUBLE_MAX;
  this->LastPickPosition[2] = VTK_DOUBLE_MAX;

  this->CurrentHandle = nullptr;
  this->CurrentHandleIndex = -1;
  this->FirstSelected = true;

  this->Transform = vtkTransform::New();

  // Set up the initial properties
  this->HandleProperty = nullptr;
  this->SelectedHandleProperty = nullptr;
  this->LineProperty = nullptr;
  this->SelectedLineProperty = nullptr;
  this->CreateDefaultProperties();

  this->Centroid[0] = 0.0;
  this->Centroid[1] = 0.0;
  this->Centroid[2] = 0.0;

  this->TranslationAxis = Axis::NONE;
  this->Visibility = true;

  vtkPolyDataMapper* lineMapper = vtkPolyDataMapper::New();
  lineMapper->SetInputConnection(this->Skeleton->GetOutputPort());
  lineMapper->SetResolveCoincidentTopologyToPolygonOffset();

  this->LineActor->SetMapper(lineMapper);
  lineMapper->Delete();

}

//------------------------------------------------------------------------------
vtkMySkeletonRepresentation::~vtkMySkeletonRepresentation()
{
  this->LineActor->Delete();

  for (int i = 0; i < this->NumberOfHandles; ++i)
  {
    this->HandleGeometry[i]->Delete();
    this->Handle[i]->Delete();
  }
  delete[] this->Handle;
  delete[] this->HandleGeometry;

  this->HandlePicker->Delete();
  this->LinePicker->Delete();

  if (this->HandleProperty)
  {
    this->HandleProperty->Delete();
  }
  if (this->SelectedHandleProperty)
  {
    this->SelectedHandleProperty->Delete();
  }
  if (this->LineProperty)
  {
    this->LineProperty->Delete();
  }
  if (this->SelectedLineProperty)
  {
    this->SelectedLineProperty->Delete();
  }

  this->Transform->Delete();
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::SetDirectionalLine(bool val)
{
  if (this->DirectionalLine == val)
  {
    return;
  }

  this->DirectionalLine = val;
  this->Modified();

  if (this->NumberOfHandles < 2)
  {
    return;
  }

  this->HandleGeometry[this->NumberOfHandles - 1]->SetUseSphere(true);
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::SetClosed(vtkTypeBool closed)
{
  if (this->Closed == closed)
  {
    return;
  }
  this->Closed = closed;

  this->BuildRepresentation();
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::RegisterPickers()
{
  vtkPickingManager* pm = this->GetPickingManager();
  if (!pm)
  {
    return;
  }
  pm->AddPicker(this->HandlePicker, this);
  pm->AddPicker(this->LinePicker, this);
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::SetHandlePosition(int handle, double x, double y, double z)
{
  if (handle < 0 || handle >= this->NumberOfHandles)
  {
    vtkErrorMacro(<< "vtkMySkeletonRepresentation: handle index out of range.");
    return;
  }
  this->HandleGeometry[handle]->SetCenter(x, y, z);
  this->HandleGeometry[handle]->Update();
  if (this->ProjectToPlane)
  {
    this->ProjectPointsToPlane();
  }
  this->BuildRepresentation();
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::SetHandlePosition(int handle, double xyz[3])
{
  this->SetHandlePosition(handle, xyz[0], xyz[1], xyz[2]);
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::GetHandlePosition(int handle, double xyz[3])
{
  if (handle < 0 || handle >= this->NumberOfHandles)
  {
    vtkErrorMacro(<< "vtkMySkeletonRepresentation: handle index out of range.");
    return;
  }
  this->HandleGeometry[handle]->GetCenter(xyz);
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::GetHandlePosition(int handle, double& x, double& y, double& z)
{
  double xyz[3];
  if (handle < 0 || handle >= this->NumberOfHandles)
  {
    vtkErrorMacro(<< "vtkMySkeletonRepresentation: handle index out of range.");
    return;
  }
  this->HandleGeometry[handle]->GetCenter(xyz);
  x = xyz[0];
  y = xyz[1];
  z = xyz[2];
}

//------------------------------------------------------------------------------
double* vtkMySkeletonRepresentation::GetHandlePosition(int handle)
{
  if (handle < 0 || handle >= this->NumberOfHandles)
  {
    vtkErrorMacro(<< "vtkMySkeletonRepresentation: handle index out of range.");
    return nullptr;
  }

  return this->HandleGeometry[handle]->GetCenter();
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::ProjectPointsToPlane()
{
  if (this->ProjectionNormal == VTK_PROJECTION_OBLIQUE)
  {
    if (this->PlaneSource != nullptr)
    {
      this->ProjectPointsToObliquePlane();
    }
    else
    {
      vtkGenericWarningMacro(<< "Set the plane source for oblique projections...");
    }
  }
  else
  {
    this->ProjectPointsToOrthoPlane();
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::ProjectPointsToObliquePlane()
{
  double o[3];
  double u[3];
  double v[3];

  this->PlaneSource->GetPoint1(u);
  this->PlaneSource->GetPoint2(v);
  this->PlaneSource->GetOrigin(o);

  int i;
  for (i = 0; i < 3; ++i)
  {
    u[i] = u[i] - o[i];
    v[i] = v[i] - o[i];
  }
  vtkMath::Normalize(u);
  vtkMath::Normalize(v);

  double o_dot_u = vtkMath::Dot(o, u);
  double o_dot_v = vtkMath::Dot(o, v);
  double fac1;
  double fac2;
  double ctr[3];
  for (i = 0; i < this->NumberOfHandles; ++i)
  {
    this->HandleGeometry[i]->GetCenter(ctr);
    fac1 = vtkMath::Dot(ctr, u) - o_dot_u;
    fac2 = vtkMath::Dot(ctr, v) - o_dot_v;
    ctr[0] = o[0] + fac1 * u[0] + fac2 * v[0];
    ctr[1] = o[1] + fac1 * u[1] + fac2 * v[1];
    ctr[2] = o[2] + fac1 * u[2] + fac2 * v[2];
    this->HandleGeometry[i]->SetCenter(ctr);
    this->HandleGeometry[i]->Update();
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::ProjectPointsToOrthoPlane()
{
  double ctr[3];
  for (int i = 0; i < this->NumberOfHandles; ++i)
  {
    this->HandleGeometry[i]->GetCenter(ctr);
    ctr[this->ProjectionNormal] = this->ProjectionPosition;
    this->HandleGeometry[i]->SetCenter(ctr);
    this->HandleGeometry[i]->Update();
  }
}

//------------------------------------------------------------------------------
int vtkMySkeletonRepresentation::GetHandleIndex(vtkProp* prop)
{
  auto iter =
    std::find(this->Handle, this->Handle + this->NumberOfHandles, static_cast<vtkActor*>(prop));
  return (iter != this->Handle + NumberOfHandles)
    ? static_cast<int>(std::distance(this->Handle, iter))
    : -1;
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::SetCurrentHandleIndex(int index)
{
  if (index < -1 || index >= this->NumberOfHandles)
  {
    index = -1;
  }

  if (index != this->CurrentHandleIndex)
  {
    this->CurrentHandleIndex = index;
    this->HighlightHandle(index == -1 ? nullptr : this->Handle[index]);
  }
}

//------------------------------------------------------------------------------
int vtkMySkeletonRepresentation::HighlightHandle(vtkProp* prop)
{
  // First unhighlight anything picked
  if (this->CurrentHandle)
  {
    this->CurrentHandle->SetProperty(this->HandleProperty);
  }

  this->CurrentHandle = static_cast<vtkActor*>(prop);

  if (this->CurrentHandle)
  {
    this->CurrentHandle->SetProperty(this->SelectedHandleProperty);
    return this->GetHandleIndex(prop);
  }
  return -1;
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::HighlightLine(int highlight)
{
  if (highlight)
  {
    this->LineActor->SetProperty(this->SelectedLineProperty);
  }
  else
  {
    this->LineActor->SetProperty(this->LineProperty);
  }
}

//------------------------------------------------------------------------------
inline double dist(double x0[3], double x1[3]) {
  return sqrt((x0[0]-x1[0])*(x0[0]-x1[0])
            + (x0[1]-x1[1])*(x0[1]-x1[1])
            + (x0[2]-x1[2])*(x0[2]-x1[2]));
}

void vtkMySkeletonRepresentation::ApplySkeletonDisplacement()
{
  // loop over handles and apply deformation
  double newCtr[3];
  for(int i = 0; i < this->NumberOfHandles; ++i) {
    double* ctr = this->HandleGeometry[i]->GetCenter();
    newCtr[0] = ctr[0] + this->b_vecx(i,0);
    newCtr[1] = ctr[1] + this->b_vecx(i,1);
    newCtr[2] = ctr[2] + this->b_vecx(i,2);
    this->HandleGeometry[i]->SetCenter(newCtr);
    this->HandleGeometry[i]->Update();
  }
}

void vtkMySkeletonRepresentation::MovePoint(double* p1, double* p2)
{
  if (this->CurrentHandleIndex < 0 || this->CurrentHandleIndex >= this->NumberOfHandles)
  {
    vtkGenericWarningMacro(<< "Poly line handle index out of range.");
    return;
  }

  // Get the motion vector
  double v[3] = { 0, 0, 0 };
  // Move the center of the handle along the motion vector
  if (this->TranslationAxis == Axis::NONE)
  {
    v[0] = p2[0] - p1[0];
    v[1] = p2[1] - p1[1];
    v[2] = p2[2] - p1[2];
  }
  // Translation restriction handling
  else
  {
    v[this->TranslationAxis] = p2[this->TranslationAxis] - p1[this->TranslationAxis];
  }

  // compute bending of skeleton
  vtkPolyData* sklWidget = this->Skeleton->GetSkeleton();
  const vtkIdType nVerts = sklWidget->GetNumberOfPoints();
  const vtkIdType nEdges = sklWidget->GetNumberOfLines();

  // zero out the rhs
  this->b_vecb.fill(0.0);
  this->b_vecx.fill(0.0);
  this->b_matM.fill(0.0);
  this->b_matL.fill(0.0);
  std::fill(this->b_bcVerts.begin(), this->b_bcVerts.end(), 0);

  // apply displacement to grabbed handle
  this->b_vecb(this->CurrentHandleIndex,0) = v[0];
  this->b_vecb(this->CurrentHandleIndex,1) = v[1];
  this->b_vecb(this->CurrentHandleIndex,2) = v[2];
  this->b_bcVerts[this->CurrentHandleIndex] = 1;

  // apply displacements to neighbors
  for(auto neighbor : this->b_neighborsMap[this->CurrentHandleIndex]) {
    this->b_bcVerts[neighbor] = 1;
    this->b_vecb(neighbor,0) = v[0];
    this->b_vecb(neighbor,1) = v[1];
    this->b_vecb(neighbor,2) = v[2];
  }

  // fix the fixed nodes
  for(int i = 0; i < nVerts; ++i) {
    if(this->b_fixedVerts[i])
      this->b_bcVerts[i] = 1;
  }
  // compute arc lengths of skeleton bones
  double pt0[3];
  double pt1[3];
  const double mass = 0.01;
  for(int i = 0; i < nEdges; ++i) {
    vtkSmartPointer<vtkIdList> pts = sklWidget->GetCell(i)->GetPointIds();
    int v0 = pts->GetId(0);
    int v1 = pts->GetId(1);
    sklWidget->GetPoint(v0, pt0);
    sklWidget->GetPoint(v1, pt1);
    const double length = dist(pt0, pt1);

    if(!this->b_bcVerts[v0]) {
      this->b_arcLength(v0) += mass/length;
    } else {
      this->b_arcLength(v0) = 1;
    }
    if(!this->b_bcVerts[v1]) {
      this->b_arcLength(v1) += mass/length;
    } else {
      this->b_arcLength(v1) = 1;
    }
  }

  this->b_matM = this->b_arcLength.asDiagonal();

  // form and assemble the laplacian matrix over the skeleton
  for(int e = 0; e < nEdges; ++e) {
    vtkSmartPointer<vtkIdList> pts = sklWidget->GetCell(e)->GetPointIds();
    int i = pts->GetId(0);
    int j = pts->GetId(1);
    sklWidget->GetPoint(i, pt0);
    sklWidget->GetPoint(j, pt1);
    double length = dist(pt0, pt1);
    if(this->b_bcVerts[i]) {
      this->b_matL(i, i) = 1;
    } else {
      this->b_matL(i, i) += 1 /length;
      this->b_matL(i, j) = -1/length;
    }

    if(this->b_bcVerts[j]) {
      this->b_matL(j, j) = 1;
    } else {
      this->b_matL(j, j) += 1 / length;
      this->b_matL(j, i) = -1/ length;
    }
  }

  // solve for the skeleton deformation
  this->b_dec.compute(this->b_matL * this->b_matM * this->b_matL);
  this->b_vecx = this->b_dec.solve(this->b_vecb);

  ApplySkeletonDisplacement();
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::Translate(double* p1, double* p2)
{
  // Get the motion vector
  double v[3] = { 0, 0, 0 };
  // Move the center of the handle along the motion vector
  if (this->TranslationAxis == Axis::NONE)
  {
    v[0] = p2[0] - p1[0];
    v[1] = p2[1] - p1[1];
    v[2] = p2[2] - p1[2];
  }
  // Translation restriction handling
  else
  {
    // this->TranslationAxis in [0,2]
    assert(this->TranslationAxis > -1 && this->TranslationAxis < 3 &&
      "this->TranslationAxis shoud be in [0,2]");
    v[this->TranslationAxis] = p2[this->TranslationAxis] - p1[this->TranslationAxis];
  }

  double newCtr[3];
  for (int i = 0; i < this->NumberOfHandles; ++i)
  {
    double* ctr = this->HandleGeometry[i]->GetCenter();
    for (int j = 0; j < 3; ++j)
    {
      newCtr[j] = ctr[j] + v[j];
    }
    this->HandleGeometry[i]->SetCenter(newCtr);
    this->HandleGeometry[i]->Update();
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::Scale(double* p1, double* p2, int vtkNotUsed(X), int Y)
{
  // Get the motion vector
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  double center[3] = { 0.0, 0.0, 0.0 };
  double avgdist = 0.0;
  double* prevctr = this->HandleGeometry[0]->GetCenter();
  double* ctr;

  center[0] += prevctr[0];
  center[1] += prevctr[1];
  center[2] += prevctr[2];

  int i;
  for (i = 1; i < this->NumberOfHandles; ++i)
  {
    ctr = this->HandleGeometry[i]->GetCenter();
    center[0] += ctr[0];
    center[1] += ctr[1];
    center[2] += ctr[2];
    avgdist += sqrt(vtkMath::Distance2BetweenPoints(ctr, prevctr));
    prevctr = ctr;
  }

  avgdist /= this->NumberOfHandles;

  center[0] /= this->NumberOfHandles;
  center[1] /= this->NumberOfHandles;
  center[2] /= this->NumberOfHandles;

  // Compute the scale factor
  double sf = vtkMath::Norm(v) / avgdist;
  if (Y > this->LastEventPosition[1])
  {
    sf = 1.0 + sf;
  }
  else
  {
    sf = 1.0 - sf;
  }

  // Move the handle points
  double newCtr[3];
  for (i = 0; i < this->NumberOfHandles; ++i)
  {
    ctr = this->HandleGeometry[i]->GetCenter();
    for (int j = 0; j < 3; ++j)
    {
      newCtr[j] = sf * (ctr[j] - center[j]) + center[j];
    }
    this->HandleGeometry[i]->SetCenter(newCtr);
    this->HandleGeometry[i]->Update();
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::Spin(double* p1, double* p2, double* vpn)
{
  // Mouse motion vector in world space
  double v[3];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  // Axis of rotation
  double axis[3] = { 0.0, 0.0, 0.0 };

  if (this->ProjectToPlane)
  {
    if (this->ProjectionNormal == VTK_PROJECTION_OBLIQUE)
    {
      if (this->PlaneSource != nullptr)
      {
        double* normal = this->PlaneSource->GetNormal();
        axis[0] = normal[0];
        axis[1] = normal[1];
        axis[2] = normal[2];
        vtkMath::Normalize(axis);
      }
      else
      {
        axis[0] = 1.;
      }
    }
    else
    {
      axis[this->ProjectionNormal] = 1.;
    }
  }
  else
  {
    // Create axis of rotation and angle of rotation
    vtkMath::Cross(vpn, v, axis);
    if (vtkMath::Normalize(axis) == 0.0)
    {
      return;
    }
  }

  // Radius vector (from mean center to cursor position)
  double rv[3] = { p2[0] - this->Centroid[0], p2[1] - this->Centroid[1],
    p2[2] - this->Centroid[2] };

  // Distance between center and cursor location
  double rs = vtkMath::Normalize(rv);

  // Spin direction
  double ax_cross_rv[3];
  vtkMath::Cross(axis, rv, ax_cross_rv);

  // Spin angle
  double theta = 360.0 * vtkMath::Dot(v, ax_cross_rv) / rs;

  // Manipulate the transform to reflect the rotation
  this->Transform->Identity();
  this->Transform->Translate(this->Centroid[0], this->Centroid[1], this->Centroid[2]);
  this->Transform->RotateWXYZ(theta, axis);
  this->Transform->Translate(-this->Centroid[0], -this->Centroid[1], -this->Centroid[2]);

  // Set the handle points
  double newCtr[3];
  double ctr[3];
  for (int i = 0; i < this->NumberOfHandles; ++i)
  {
    this->HandleGeometry[i]->GetCenter(ctr);
    this->Transform->TransformPoint(ctr, newCtr);
    this->HandleGeometry[i]->SetCenter(newCtr);
    this->HandleGeometry[i]->Update();
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::CreateDefaultProperties()
{
  this->HandleProperty = vtkProperty::New();
  this->HandleProperty->SetColor(1, 1, 1);

  this->SelectedHandleProperty = vtkProperty::New();
  this->SelectedHandleProperty->SetColor(1, 0, 0);

  this->LineProperty = vtkProperty::New();
  this->LineProperty->SetRepresentationToWireframe();
  this->LineProperty->SetAmbient(1.0);
  this->LineProperty->SetColor(1.0, 1.0, 0.0);
  this->LineProperty->SetLineWidth(2.0);

  this->SelectedLineProperty = vtkProperty::New();
  this->SelectedLineProperty->SetRepresentationToWireframe();
  this->SelectedLineProperty->SetAmbient(1.0);
  this->SelectedLineProperty->SetAmbientColor(0.0, 1.0, 0.0);
  this->SelectedLineProperty->SetLineWidth(2.0);
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::SetProjectionPosition(double position)
{
  this->ProjectionPosition = position;
  if (this->ProjectToPlane)
  {
    this->ProjectPointsToPlane();
  }
  this->BuildRepresentation();
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::SetPlaneSource(vtkPlaneSource* plane)
{
  if (this->PlaneSource == plane)
  {
    return;
  }
  this->PlaneSource = plane;
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::Initialize()
{
  int i;
  for (i = 0; i < this->NumberOfHandles; ++i)
  {
    this->HandlePicker->DeletePickList(this->Handle[i]);
    this->HandleGeometry[i]->Delete();
    this->Handle[i]->Delete();
  }

  this->NumberOfHandles = 0;

  delete[] this->Handle;
  delete[] this->HandleGeometry;
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::SizeHandles()
{
  if (this->NumberOfHandles > 0)
  {
    double radius = this->SizeHandlesInPixels(1.5, this->HandleGeometry[0]->GetCenter());
    for (int i = 0; i < this->NumberOfHandles; ++i)
    {
      this->HandleGeometry[i]->SetRadius(radius);
    }
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::CalculateCentroid()
{
  this->Centroid[0] = 0.0;
  this->Centroid[1] = 0.0;
  this->Centroid[2] = 0.0;

  double ctr[3];
  for (int i = 0; i < this->NumberOfHandles; ++i)
  {
    this->HandleGeometry[i]->GetCenter(ctr);
    this->Centroid[0] += ctr[0];
    this->Centroid[1] += ctr[1];
    this->Centroid[2] += ctr[2];
  }

  this->Centroid[0] /= this->NumberOfHandles;
  this->Centroid[1] /= this->NumberOfHandles;
  this->Centroid[2] /= this->NumberOfHandles;
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::ReleaseGraphicsResources(vtkWindow* win)
{
  this->LineActor->ReleaseGraphicsResources(win);
  for (int cc = 0; cc < this->NumberOfHandles; cc++)
  {
    this->Handle[cc]->ReleaseGraphicsResources(win);
  }
}

//------------------------------------------------------------------------------
int vtkMySkeletonRepresentation::RenderOpaqueGeometry(vtkViewport* win)
{
  this->BuildRepresentation();

  int count = 0;
  count += this->LineActor->RenderOpaqueGeometry(win);
  for (int cc = 0; cc < this->NumberOfHandles; cc++)
  {
    count += this->Handle[cc]->RenderOpaqueGeometry(win);
  }
  return count;
}

//------------------------------------------------------------------------------
int vtkMySkeletonRepresentation::RenderTranslucentPolygonalGeometry(vtkViewport* win)
{
  int count = 0;
  count += this->LineActor->RenderTranslucentPolygonalGeometry(win);
  for (int cc = 0; cc < this->NumberOfHandles; cc++)
  {
    count += this->Handle[cc]->RenderTranslucentPolygonalGeometry(win);
  }
  return count;
}

//------------------------------------------------------------------------------
int vtkMySkeletonRepresentation::RenderOverlay(vtkViewport* win)
{
  int count = 0;
  count += this->LineActor->RenderOverlay(win);
  for (int cc = 0; cc < this->NumberOfHandles; cc++)
  {
    count += this->Handle[cc]->RenderOverlay(win);
  }
  return count;
}

//------------------------------------------------------------------------------
vtkTypeBool vtkMySkeletonRepresentation::HasTranslucentPolygonalGeometry()
{
  this->BuildRepresentation();
  int count = 0;
  count |= this->LineActor->HasTranslucentPolygonalGeometry();
  for (int cc = 0; cc < this->NumberOfHandles; cc++)
  {
    count |= this->Handle[cc]->HasTranslucentPolygonalGeometry();
  }
  return count;
}

//------------------------------------------------------------------------------
int vtkMySkeletonRepresentation::ComputeInteractionState(int X, int Y, int vtkNotUsed(modify))
{
  this->InteractionState = vtkMySkeletonRepresentation::Outside;
  if (!this->Renderer || !this->Renderer->IsInViewport(X, Y))
  {
    return this->InteractionState;
  }

  // Try and pick a handle first. This allows the picking of the handle even
  // if it is "behind" the poly line.
  int handlePicked = 0;

  vtkAssemblyPath* path = this->GetAssemblyPath(X, Y, 0., this->HandlePicker);

  // always get pick position
  this->HandlePicker->GetPickPosition(this->LastPickPosition);

  if (path != nullptr)
  {
    this->ValidPick = 1;
    this->InteractionState = vtkMySkeletonRepresentation::OnHandle;
    this->SetCurrentHandleIndex(this->GetHandleIndex(path->GetFirstNode()->GetViewProp()));
    handlePicked = 1;
    this->FirstSelected = (this->CurrentHandleIndex == 0);
  }
  else
  {
    this->SetCurrentHandleIndex(-1);
  }

  if (!handlePicked)
  {
    path = this->GetAssemblyPath(X, Y, 0., this->LinePicker);

    if (path != nullptr)
    {
      this->ValidPick = 1;
      this->LinePicker->GetPickPosition(this->LastPickPosition);
      this->HighlightLine(1);
      this->InteractionState = vtkMySkeletonRepresentation::OnLine;
    }
    else
    {
      this->HighlightLine(0);
    }
  }
  else
  {
    this->HighlightLine(0);
  }

  return this->InteractionState;
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::StartWidgetInteraction(double e[2])
{
  // Store the start position
  this->StartEventPosition[0] = e[0];
  this->StartEventPosition[1] = e[1];
  this->StartEventPosition[2] = 0.0;

  // Store the start position
  this->LastEventPosition[0] = e[0];
  this->LastEventPosition[1] = e[1];
  this->LastEventPosition[2] = 0.0;

  this->ComputeInteractionState(static_cast<int>(e[0]), static_cast<int>(e[1]), 0);
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::WidgetInteraction(double e[2])
{
  // Convert events to appropriate coordinate systems
  vtkCamera* camera = this->Renderer->GetActiveCamera();
  if (!camera)
  {
    return;
  }
  double focalPoint[4], pickPoint[4], prevPickPoint[4];
  double z, vpn[3];

  // Compute the two points defining the motion vector
  vtkInteractorObserver::ComputeWorldToDisplay(this->Renderer, this->LastPickPosition[0],
    this->LastPickPosition[1], this->LastPickPosition[2], focalPoint);
  z = focalPoint[2];
  vtkInteractorObserver::ComputeDisplayToWorld(
    this->Renderer, this->LastEventPosition[0], this->LastEventPosition[1], z, prevPickPoint);
  vtkInteractorObserver::ComputeDisplayToWorld(this->Renderer, e[0], e[1], z, pickPoint);

  // Process the motion
  if (this->InteractionState == vtkMySkeletonRepresentation::Moving && this->CurrentHandleIndex != -1)
  {
    this->MovePoint(prevPickPoint, pickPoint);
  }

  if (this->ProjectToPlane)
  {
    this->ProjectPointsToPlane();
  }

  this->BuildRepresentation();

  vtkXMLPolyDataWriter* writer = vtkXMLPolyDataWriter::New();
  std::string cache_path_skelMod2 = std::string(CACHE_PATH) + std::string("/skeletonMod.vtp");
  writer->SetFileName ( cache_path_skelMod2.c_str() );
  writer->SetInputData ( this->Skeleton->GetSkeleton() );
  writer->Write();
  writer->Delete();

  // Store the position
  this->LastEventPosition[0] = e[0];
  this->LastEventPosition[1] = e[1];
  this->LastEventPosition[2] = 0.0;
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::EndWidgetInteraction(double[2])
{
  this->HighlightLine(0);
  this->InteractionState = vtkMySkeletonRepresentation::Outside;
}

//------------------------------------------------------------------------------
double* vtkMySkeletonRepresentation::GetBounds()
{
  this->BuildRepresentation();

  vtkBoundingBox bbox;
  bbox.AddBounds(this->LineActor->GetBounds());
  for (int cc = 0; cc < this->NumberOfHandles; cc++)
  {
    bbox.AddBounds(this->HandleGeometry[cc]->GetOutput()->GetBounds());
  }
  bbox.GetBounds(this->Bounds);
  return this->Bounds;
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::SetLineColor(double r, double g, double b)
{
  this->GetLineProperty()->SetColor(r, g, b);
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::GetPolyData(vtkPolyData* pd)
{
  this->Skeleton->Update();
  pd->ShallowCopy(this->Skeleton->GetOutput());
}

//------------------------------------------------------------------------------
vtkDoubleArray* vtkMySkeletonRepresentation::GetHandlePositions()
{
  return vtkArrayDownCast<vtkDoubleArray>(this->Skeleton->GetPoints()->GetData());
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::BuildRepresentation()
{
  this->ValidPick = 1;
  // TODO: Avoid unnecessary rebuilds.
  // Handles have changed position, re-compute the points
  vtkPoints* points = this->Skeleton->GetPoints();
  if (points->GetNumberOfPoints() != this->NumberOfHandles)
  {
    points->SetNumberOfPoints(this->NumberOfHandles);
  }

  vtkBoundingBox bbox;
  for (int i = 0; i < this->NumberOfHandles; ++i)
  {
    double pt[3];
    this->HandleGeometry[i]->GetCenter(pt);
    points->SetPoint(i, pt);
    bbox.AddPoint(pt);
  }
  // this->Skeleton->SetClosed(this->Closed);
  this->Skeleton->Modified();
  points->Modified();

  double bounds[6];
  bbox.GetBounds(bounds);
  this->InitialLength = sqrt((bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) +
    (bounds[3] - bounds[2]) * (bounds[3] - bounds[2]) +
    (bounds[5] - bounds[4]) * (bounds[5] - bounds[4]));
  this->SizeHandles();

  // differentiate our polyline from the default
  SetLineColor(0,1,0);
  this->Modified();
}

//------------------------------------------------------------------------------
void vtkMySkeletonRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  if (this->HandleProperty)
  {
    os << indent << "Handle Property: " << this->HandleProperty << "\n";
  }
  else
  {
    os << indent << "Handle Property: (none)\n";
  }
  if (this->SelectedHandleProperty)
  {
    os << indent << "Selected Handle Property: " << this->SelectedHandleProperty << "\n";
  }
  else
  {
    os << indent << "Selected Handle Property: (none)\n";
  }
  if (this->LineProperty)
  {
    os << indent << "Line Property: " << this->LineProperty << "\n";
  }
  else
  {
    os << indent << "Line Property: (none)\n";
  }
  if (this->SelectedLineProperty)
  {
    os << indent << "Selected Line Property: " << this->SelectedLineProperty << "\n";
  }
  else
  {
    os << indent << "Selected Line Property: (none)\n";
  }

  os << indent << "Project To Plane: " << (this->ProjectToPlane ? "On" : "Off") << "\n";
  os << indent << "Projection Normal: " << this->ProjectionNormal << "\n";
  os << indent << "Projection Position: " << this->ProjectionPosition << "\n";
  os << indent << "Number Of Handles: " << this->NumberOfHandles << "\n";
  os << indent << "Closed: " << (this->Closed ? "On" : "Off") << "\n";
  os << indent << "InteractionState: " << this->InteractionState << endl;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
vtkStandardNewMacro(vtkMySkeletonRepresentation::HandleSource);

//------------------------------------------------------------------------------
vtkMySkeletonRepresentation::HandleSource::HandleSource()
  : UseSphere(true)
  , Radius(0.5)
{
  this->Center[0] = 0.0;
  this->Center[1] = 0.0;
  this->Center[2] = 0.0;

  this->Direction[0] = 1.0;
  this->Direction[1] = 0.0;
  this->Direction[2] = 0.0;

  this->SetNumberOfInputPorts(0);
}

//------------------------------------------------------------------------------
int vtkMySkeletonRepresentation::HandleSource::RequestData(
  vtkInformation*, vtkInformationVector**, vtkInformationVector* outputVector)
{
  auto output = vtkPolyData::GetData(outputVector);
  if (this->UseSphere)
  {
    vtkNew<vtkSphereSource> sphere;
    sphere->SetRadius(this->Radius);
    sphere->SetCenter(this->Center);
    sphere->SetThetaResolution(16);
    sphere->SetPhiResolution(8);
    sphere->Update();
    output->ShallowCopy(sphere->GetOutput(0));
  }
  else
  {
    vtkNew<vtkConeSource> cone;
    cone->SetRadius(this->Radius);
    cone->SetCenter(this->Center);
    cone->SetHeight(2.8 * this->Radius);
    cone->SetResolution(16);
    cone->SetDirection(this->Direction);
    cone->Update();
    output->ShallowCopy(cone->GetOutput(0));
  }
  return 1;
}

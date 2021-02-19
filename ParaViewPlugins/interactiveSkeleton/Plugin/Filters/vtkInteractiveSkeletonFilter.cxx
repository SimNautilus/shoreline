/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkInteractiveSkeletonFilter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include <ostream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <sstream>

#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkPolyLine.h"
#include "vtkBoundingBox.h"
#include "vtkDataObject.h"
#include <vtkCenterOfMass.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include "vtkXMLPolyDataReader.h"
#include "vtkObjectFactory.h"

// #include <Eigen/Core>
// #include <igl/embree/bone_heat.h>

#include "vtkInteractiveSkeletonFilter.h"



vtkStandardNewMacro(vtkInteractiveSkeletonFilter);

////////////////////////////////////////////////////////////////////////////////
vtkInteractiveSkeletonFilter::vtkInteractiveSkeletonFilter(){
  SetNumberOfInputPorts(2);
  SetNumberOfOutputPorts(1);
  this->Points = vtkPoints::New();
  this->HandlePositions = vtkDoubleArray::New();
}

////////////////////////////////////////////////////////////////////////////////
vtkInteractiveSkeletonFilter::~vtkInteractiveSkeletonFilter(){
  if(this->HandlePositions)
    this->HandlePositions->Delete();
  if(this->Points)
    this->Points->Delete();
}

////////////////////////////////////////////////////////////////////////////////
// Creates call CGAL algorithms to generate a skeleton of the input mesh
int vtkInteractiveSkeletonFilter::RequestData(
  vtkInformation *,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  //  Get the input and output data objects.
  //  Get the info objects
  vtkInformation *inInfoMsh = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfoSkl = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  //  Get the input mesh
  vtkPolyData *mshData = vtkPolyData::SafeDownCast(inInfoMsh->Get(
    vtkDataObject::DATA_OBJECT()));
  vtkPolyData *sklData = vtkPolyData::SafeDownCast(inInfoSkl->Get(
    vtkDataObject::DATA_OBJECT()));


  //////////////////////////////////////////////////////////////////////////////
  // Read the representation of the widget from file
  vtkXMLPolyDataReader* reader = vtkXMLPolyDataReader::New();
  reader->SetFileName ("/Users/coreylee/Git/cache/skeletonMod.vtp");
  reader->Update();

  vtkPolyData* sklWidget = reader->GetOutput();
  this->Points->DeepCopy(sklWidget->GetPoints());
  this->Modified();
  reader->Delete();

  const vtkIdType nVerts = mshData->GetNumberOfPoints();
  const vtkIdType nFaces = mshData->GetNumberOfCells();

  const vtkIdType nSkelVerts = sklData->GetNumberOfPoints();
  const vtkIdType nSkelEdges = sklData->GetNumberOfCells();

  /***************************************************************
  * create libigl data structure to start computing bone weights *
  ***************************************************************/
  // define libigl data structures for storing geometric data
  // Eigen::MatrixXd V(nVerts, 3);
  // Eigen::MatrixXi F(nFaces, 3);
  // Eigen::MatrixXd C(nSkelVerts, 3);
  // Eigen::VectorXi P(nSkelVerts);
  // Eigen::MatrixXi BE(nSkelEdges, 2);
  // Eigen::MatrixXd W; // skinning weights
  //
  // // Extract points
  // for (vtkIdType i = 0; i < nVerts; ++i) {
  //   double coords[3];
  //   mshData->GetPoint(i, coords);
  //   V(i,0) = coords[0];
  //   V(i,1) = coords[1];
  //   V(i,2) = coords[2];
  // }
  //
  // // Extract cells
  // for (vtkIdType i = 0; i<nFaces; ++i) {
  //   vtkCell* cell_ptr = mshData->GetCell(i);
  //   vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
  //   if(nb_vertices != 3) {
  //     vtkErrorMacro("The input mesh must be triangulated ");
  //   }
  //   for (vtkIdType k=0; k<nb_vertices; ++k)
  //     F(i,k) = cell_ptr->GetPointId(k);
  // }
  //
  // // create the skeleton in libigl terms
  // double pt[3];
  // for(int i = 0; i < nSkelVerts; ++i) {
  //   sklData->GetPoints()->GetPoint(i, pt);
  //   C(i, 0) = pt[0];
  //   C(i, 1) = pt[1];
  //   C(i, 2) = pt[2];
  //   P(i) = i;
  // }
  // for(int i = 0; i < nSkelEdges; ++i) {
  //   vtkSmartPointer<vtkIdList> pts = sklData->GetCell(i)->GetPointIds();
  //   BE(i,0) = pts->GetId(0);
  //   BE(i,1) = pts->GetId(1);
  // }
  //
  // // create skinning weights!
  // if(!igl::embree::bone_heat(V, F, C, P, BE, Eigen::MatrixXi(), W))
  // {
  //   vtkErrorMacro("Weight computation failed");
  // }

  /********************************************
  *       pass mesh through to output 0       *
  ********************************************/
  vtkPolyData *meshOut = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  meshOut->DeepCopy(mshData);

  //////////////////////////////////////////////////////////////////////////////
  // move mesh according to joint association
  // get the features array


  // apply the displacement field to the output
   vtkNew<vtkDoubleArray> const dispsOut;
   dispsOut->SetNumberOfComponents(3);
   dispsOut->SetName("displacements");
   dispsOut->Allocate(nVerts);

  if(meshOut->GetPointData()->HasArray("joint"))
  {
    vtkDataArray* jointData = meshOut->GetPointData()->GetArray("joint");
    double pt0[3] = {0,0,0};
    double pt1[3] = {0,0,0};
    double del[3] = {0,0,0};
    double mPt[3] = {0,0,0};

    for(int i = 0; i < nVerts; ++i) {
      // get the joint associated with this point
      int joint = jointData->GetTuple1(i);

      // get the displacement of that joint
      sklData->GetPoints()->GetPoint(joint, pt0);
      Points->GetPoint(joint, pt1);
      del[0] = pt1[0] - pt0[0];
      del[1] = pt1[1] - pt0[1];
      del[2] = pt1[2] - pt0[2];

      meshOut->GetPoints()->GetPoint(i, mPt);
      mPt[0] += del[0];
      mPt[1] += del[1];
      mPt[2] += del[2];
      // double pt[3] = {V(i,0), V(i,1), V(i,2)};
      meshOut->GetPoints()->SetPoint(i, mPt);

      if(meshOut->GetPointData()->HasArray("displacements")) {
        vtkDataArray* dispsIn = meshOut->GetPointData()->GetArray("displacements");
        double dIn[3];
        dispsIn->GetTuple(i, dIn);
        dispsOut->InsertNextTuple3(dIn[0] + del[0], dIn[1] + del[1], dIn[2] + del[2]);
      } else {
          dispsOut->InsertNextTuple3(del[0], del[1], del[2]);
      }

    }
   } else {
     vtkErrorMacro("vtk object does not have joint array");
   }

   meshOut->GetPointData()->AddArray(dispsOut);
  //////////////////////////////////////////////////////////////////////////////
  // HandlePositions = GetHandlePositions();
  // if(this->Points)
  // {
  //   if(this->Points->GetNumberOfPoints() == sklData->GetNumberOfPoints())
  //     sklData->SetPoints(this->Points);
  // }




  this->Modified();
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
void vtkInteractiveSkeletonFilter::SetPoint(vtkIdType id, double x, double y, double z)
{
  if (!this->Points)
  {
    return;
  }

  if (id >= this->Points->GetNumberOfPoints())
  {
    this->Points->Resize(id+1);
    // vtkErrorMacro(<< "point id " << id << " is larger than the number of points");
    // return;
  }

  this->Points->SetPoint(id, x, y, z);
  this->Modified();
  this->DebugOn();
  vtkDebugMacro(<<"hey we moved a point");
}

////////////////////////////////////////////////////////////////////////////////
void vtkInteractiveSkeletonFilter::SetPoints(vtkDoubleArray* pts)
{
  if(this->Points && pts){
    this->Points->SetData(pts);
    this->Modified();
  }
  this->DebugOn();
  vtkDebugMacro(<<"hey we moved some points");
}

////////////////////////////////////////////////////////////////////////////////
void vtkInteractiveSkeletonFilter::SetNumberOfPoints(int npts)
{
  if(this->Points)
  {
    this->Points->Resize(npts);
    this->Modified();
  }
}

////////////////////////////////////////////////////////////////////////////////
void vtkInteractiveSkeletonFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

////////////////////////////////////////////////////////////////////////////////
int vtkInteractiveSkeletonFilter::FillInputPortInformation(
  int, vtkInformation* info)
{
  // both input ports are vtkpolydata
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
int vtkInteractiveSkeletonFilter::FillOutputPortInformation(
  int, vtkInformation *info)
{
  // Always returns a vtkPolyData
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
int vtkInteractiveSkeletonFilter::RequestInformation(
  vtkInformation *,
  vtkInformationVector ** inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Sets the bounds of the output.
  outInfo->Set(vtkDataObject::BOUNDING_BOX(),
    inInfo->Get(vtkDataObject::BOUNDING_BOX()), 6);

  return 1;
}

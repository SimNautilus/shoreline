#include <ostream>
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>

#include "vtkInformationVector.h"
#include "scaleNormalFilter.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"

#include <igl/harmonic.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/doublearea.h>
#include <igl/barycenter.h>
#include <igl/per_corner_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

// Declare the plugin
vtkStandardNewMacro(scaleNormalFilter);
////////////////////////////////////////////////////////////////////////////////
// Constructor
// Fills the number of input and output objects.
// Initializes the members that need it.
scaleNormalFilter::scaleNormalFilter() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  this->DebugOn();
}

////////////////////////////////////////////////////////////////////////////////
// Gets the input
// Creates CGAL::Surface_mesh from vtkPolydata
// Calls the CGAL::isotropic_remeshing algorithm
// Fills the output vtkPolyData from the result.
int scaleNormalFilter::RequestData(
  vtkInformation *,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
  startTime = std::chrono::system_clock::now();
  //  Get the input and output data objects.
  //  Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //  Get the input
  vtkPolyData *polydata = vtkPolyData::SafeDownCast(inInfo->Get(
    vtkDataObject::DATA_OBJECT()));


  // get the features array
  vtkDataArray * featureData = nullptr;
  if(polydata->GetCellData()->HasArray("features"))
    featureData = polydata->GetCellData()->GetArray("features");
  if(featureData == nullptr) {
    vtkErrorMacro("vtk object does not have features array");
  }
  /********************************************
   * Create a igl surface mesh from the input mesh *
   ********************************************/
  //  Get nb of points and cells
  vtkIdType nVerts = polydata->GetNumberOfPoints();
  vtkIdType nFaces = polydata->GetNumberOfCells();


  // define libigl data structures for storing geometric data
  Eigen::MatrixXd V(nVerts, 3);
  Eigen::MatrixXi F(nFaces, 3);
  Eigen::MatrixXi featureFace(nFaces, 1); // feature identifiers by face

  // Extract points
  for (vtkIdType i = 0; i < nVerts; ++i) {
    double coords[3];
    polydata->GetPoint(i, coords);
    V(i,0) = coords[0];
    V(i,1) = coords[1];
    V(i,2) = coords[2];
  }

  // Extract cells
  for (vtkIdType i = 0; i<nFaces; ++i) {
    vtkCell* cell_ptr = polydata->GetCell(i);
    vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
    if(nb_vertices != 3) {
      vtkErrorMacro("The input mesh must be triangulated ");
    }
    featureFace(i,0) = *(featureData->GetTuple(i));
    for (vtkIdType k=0; k<nb_vertices; ++k)
      F(i,k) = cell_ptr->GetPointId(k);
  }

  /**********************************
   * Compute normals *
   **********************************/
  Eigen::MatrixXd N(nVerts, 3);
  Eigen::MatrixXd N_faces(nFaces, 3);
  // no real reason to use this particular normal algorithm. there are others
  igl::per_vertex_normals(V, F, N);
  igl::per_face_normals(V,F,N_faces);

  // allocate space for nodal normal vectors
  vtkNew<vtkDoubleArray> const normalsArray;
  normalsArray->SetNumberOfComponents(3);
  normalsArray->SetName("Normals");
  normalsArray->Allocate(nVerts);


  /**********************************************
  * Convert face regions to boundary conditions *
  **********************************************/
  // collect from list, regions to scale
  std::string srString = regionsToScale;
  std::string frString = regionsToFix;
  std::stringstream srStream(srString);
  std::stringstream frStream(frString);

  // get those lists into a vector of ints
  std::vector<int> srVec, frVec;
  int xx = -1;
  while(srStream >> xx) {
    srVec.push_back(xx);
  }
  xx = -1;
  while(frStream >> xx) {
    frVec.push_back(xx);
  }

  enum Move {unset = -1, free = 0, fix = 1, move = 2};
  // start by mapping features to boundary conditions
  std::vector<Move> feat2BC(featureFace.maxCoeff() + 1, Move::free);
  for(auto r : srVec)
    feat2BC[r] = Move::move;

  for(auto r : frVec)
    feat2BC[r] = Move::fix;

  // feat2BC[regionToScale] = Move::move;
  // feat2BC[regionToFix[0]] = Move::fix;
  // feat2BC[regionToFix[1]] = Move::fix;

  // buzz through vertices on every face and assign the vertices to the behavior
  // of the incident face. Prioritize move > fix > free
  Eigen::VectorXi featureVert(nVerts); // feature identifiers by face
  featureVert.fill(Move::unset);
  for(int f = 0; f < nFaces; ++f) {
    for(int v = 0; v < 3; ++v) {
      if(feat2BC[featureFace(f)] == Move::move) {
       N(F(f,v),0) = N_faces(f,0);
       N(F(f,v),1) = N_faces(f,1);
       N(F(f,v),2) = N_faces(f,2);
      }
      if(featureVert(F(f,v)) != Move::move) {
        if(featureVert(F(f,v)) != Move::fix) {
          featureVert(F(f,v)) = feat2BC[featureFace(f)];
        }
      }
    }
  }

  Eigen::VectorXi b(featureVert.count()); // boundary vertex indices
  Eigen::MatrixXd D_bc(b.size(), 3); // deformation field boundary conditions
  int itr = 0;
  // loop over every vertex and get the indices of the boundary
  for(int i = 0; i < featureVert.size(); ++i) {
    if(featureVert(i) > 0) {
      b(itr) = i;
      if(featureVert(i) == Move::fix) {
        D_bc(itr, 0) = 0.0;
        D_bc(itr, 1) = 0.0;
        D_bc(itr, 2) = 0.0;
      }
      if(featureVert(i) == Move::move) {
        D_bc(itr, 0) = N(i,0) * scaleValue;
        D_bc(itr, 1) = N(i,1) * scaleValue;
        D_bc(itr, 2) = N(i,2) * scaleValue;
      }
      itr++;
    }
  }
  /*****************************
  * Apply deformation *
  *****************************/
  Eigen::MatrixXd D; // deformation field
  igl::harmonic(V,F,b,D_bc,2,D); // compute deformation field
  V = V+D; // apply deformation

  /**********************************
   * Pass the igl data to the output *
   **********************************/
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkNew<vtkPoints> const vtkVerts;
  vtkNew<vtkCellArray> const vtkFaces;
  vtkVerts->Allocate(nVerts);
  vtkFaces->Allocate(nFaces);

  /**********************************
   * Compute new normals and apply  *
   **********************************/
  igl::per_corner_normals(V, F, 20, N);

  for(int i = 0; i < nVerts; ++i) {
    vtkVerts->InsertNextPoint(V(i,0), V(i,1), V(i,2));
    normalsArray->InsertNextTuple3(N(i,0), N(i,1), N(i,2));
  }
  // apply faces to the output
  for(int i = 0; i < nFaces; ++i)
  {
    vtkNew<vtkIdList> cell;
    cell->InsertNextId(F(i,0));
    cell->InsertNextId(F(i,1));
    cell->InsertNextId(F(i,2));

    vtkFaces->InsertNextCell(cell);
  }

  output->SetPoints(vtkVerts);
  // output->GetPointData()->SetScalars(curvature);
  // output->GetPointData()->SetNormals(normalsArray);
  output->GetCellData()->AddArray(featureData);
  output->SetPolys(vtkFaces);

  // apply the displacement field to the output
  vtkNew<vtkDoubleArray> const dispsOut;
  dispsOut->SetNumberOfComponents(3);
  dispsOut->SetName("displacements");
  dispsOut->Allocate(nVerts);
  if(polydata->GetPointData()->HasArray("displacements")) {
    vtkDataArray* dispsIn = polydata->GetPointData()->GetArray("displacements");
    for(int i = 0; i < nVerts; ++i) {
      double dIn[3];
      dispsIn->GetTuple(i, dIn);
      dispsOut->InsertNextTuple3(dIn[0] + D(i,0), dIn[1] + D(i,1), dIn[2] + D(i,2));
    }
  } else {
    for(int i = 0; i < nVerts; ++i) {
      dispsOut->InsertNextTuple3(D(i,0), D(i,1), D(i,2));
    }
  }
  output->GetPointData()->AddArray(dispsOut);

  output->Squeeze();


  endTime = std::chrono::system_clock::now();
  vtkDebugMacro(<< "scaleNormal: " <<std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms");
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
void scaleNormalFilter::PrintSelf(std::ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os,indent);
  os << "regionToScale : "    << regionsToScale   << std::endl;
  os << "regionsToFix  : " << regionsToFix << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
int scaleNormalFilter::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
int scaleNormalFilter::FillOutputPortInformation(int, vtkInformation *info)
{
  // Always returns a vtkPolyData
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
int scaleNormalFilter::RequestInformation(
  vtkInformation *,
  vtkInformationVector ** inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Sets the bounds of the output.
  outInfo->Set(vtkDataObject::BOUNDING_BOX(),
    inInfo->Get(vtkDataObject::BOUNDING_BOX()), 6);

  vtkPolyData *input= vtkPolyData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));

    return 1;
}

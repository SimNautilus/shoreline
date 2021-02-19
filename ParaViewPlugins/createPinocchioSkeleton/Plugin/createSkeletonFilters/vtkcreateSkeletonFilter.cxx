/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkcreateSkeletonFilter.cxx

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

#include "pinocchioApi.h"

#include "vtkcreateSkeletonFilter.h"

#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkcreateSkeletonFilter);

////////////////////////////////////////////////////////////////////////////////
vtkcreateSkeletonFilter::vtkcreateSkeletonFilter(){
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);
}

////////////////////////////////////////////////////////////////////////////////
vtkcreateSkeletonFilter::~vtkcreateSkeletonFilter(){}

////////////////////////////////////////////////////////////////////////////////
// Gets the input
// Creates CGAL::Surface_mesh from vtkPolydata
// Calls the CGAL::isotropic_remeshing algorithm
// Fills the output vtkPolyData from the result.
int vtkcreateSkeletonFilter::RequestData(
  vtkInformation *,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  this->DebugOn();

  //  Get the input and output data objects.
  //  Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
  //  Get the input
  vtkPolyData *polydata = vtkPolyData::SafeDownCast(inInfo->Get(
    vtkDataObject::DATA_OBJECT()));

  //  Get nb of points and cells
  const vtkIdType nVerts = polydata->GetNumberOfPoints();
  const vtkIdType nFaces = polydata->GetNumberOfCells();

  // get the features array
  vtkDataArray * featureData = nullptr;
  if(polydata->GetCellData()->HasArray("features"))
    featureData = polydata->GetCellData()->GetArray("features");
  if(featureData == nullptr) {
    vtkErrorMacro("vtk object does not have features array");
  }

  // list of features on each face
  std::vector<int> featureFace(nFaces);

  for (vtkIdType i = 0; i < nFaces; ++i)
    featureFace[i] = *(featureData->GetTuple(i));

  const int nJoints = boneCount+1;
  /**********************************
  * Collect list of vertices to rig *
  **********************************/

  std::map<int, bool> rigFeature;
  // collect from list, regions to rig
  std::string rigString = regionsToRig;
  std::stringstream rigStream(rigString);

  // get those lists into a vector of ints
  if(rigString == "all") { // rig the whole object
    for(auto f : featureFace) rigFeature[f] = true;
  } else { // read list of features
    std::vector<int> rigVec;
    int xx = -1;
    while(rigStream >> xx)
      rigVec.push_back(xx);

    // start by mapping features for rigging
    for(auto r : rigVec) rigFeature[r] = true;
  }

  /**********************************
  * Convert mesh data to Pinocchio  *
  **********************************/

  // create mesh of features to pass along to Pinocchio
  std::vector<double> verts;
  std::vector<int> halfEdges;

  // create map from vertices on the mesh to verts getting rigged
  std::map<int, int> vtk2Pinocchio;

  // buzz through faces
  int pID = 0;
  for(int f = 0; f < nFaces; ++f) {
    // check if the face is on a feature to rig
    if(rigFeature.find(featureFace[f]) != rigFeature.end()) {
      // get a pointer to the face
      vtkCell* cell_ptr = polydata->GetCell(f);
      // loop over points
      for(int i = 0; i < 3; ++i) {
        int vtkVertID = cell_ptr->GetPointId(i);
        // add the vertex to the list, otherwise just add the half edge
        if(vtk2Pinocchio.find(vtkVertID) == vtk2Pinocchio.end()) {
          vtk2Pinocchio[vtkVertID] = pID;
          pID++;
          double coords[3];
          polydata->GetPoint(vtkVertID, coords);
          verts.insert(verts.end(), coords, coords+3);
        }
        // add the half edge
        halfEdges.push_back(vtk2Pinocchio[vtkVertID]);
      }
    }
  }

  std::string debug = "verts passed to pinocchio: " + to_string(verts.size()/3) + ", edges passed: " + to_string(halfEdges.size());
  vtkDebugMacro( << debug);

  // define a mesh for the Pinocchio rigger system
  std::stringstream pinocchioDebug;
  Debugging::setOutStream(pinocchioDebug);
  Mesh pMesh = Mesh(verts, halfEdges);

  vtkDebugMacro( << "PINOCCHIO: " << pinocchioDebug.str());
  pinocchioDebug.str("");
  pinocchioDebug.clear();
  debug = "verts read by pinocchio: " + to_string(pMesh.vertices.size()) + ", edges read: " + to_string(pMesh.edges.size());
  vtkDebugMacro( << debug);

  // write the obj file for references
  // pMesh.writeObj("/Users/coreylee/Git/Nelson/ImmersiveSims/ApplicationCode/paraviewPlugins/immersive/createSkeleton/build/out.obj");
  /************************************
  * Generate the initial skeleton rig *
  ************************************/
  // get the bounds for the rig
  double bnds[6] = {0,0,0,0,0,0};
  polydata->ComputeBounds();
  polydata->GetBounds(bnds);

  double xMean = 0.5 * (bnds[0] + bnds[1]);
  double yMean = 0.5 * (bnds[2] + bnds[3]);
  double zMean = 0.5 * (bnds[4] + bnds[5]);
  switch(skeletonDirection) {
    case 0:
    bnds[2] = yMean; bnds[3] = yMean; bnds[4] = zMean; bnds[5] = zMean; break;
    case 1:
    bnds[0] = xMean; bnds[1] = xMean; bnds[4] = zMean; bnds[5] = zMean; break;
    case 2:
    bnds[0] = xMean; bnds[1] = xMean; bnds[2] = yMean; bnds[3] = yMean; break;
  }

  vtkSmartPointer<vtkPoints> skelPoints = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkPolyLine> skelLine = vtkSmartPointer<vtkPolyLine>::New();
  skelLine->GetPointIds()->SetNumberOfIds(nJoints);

  // create a chain of bones
  for(int b = 0; b < nJoints; ++b) {
    double x = bnds[0] + double(b)/(boneCount) * (bnds[1] - bnds[0]);
    double y = bnds[2] + double(b)/(boneCount) * (bnds[3] - bnds[2]);
    double z = bnds[4] + double(b)/(boneCount) * (bnds[5] - bnds[4]);
    double p[3] = {x,y,z};
    skelPoints->InsertNextPoint(p);
    skelLine->GetPointIds()->SetId(b,b);
  }

  // Create a cell array to store the lines in and add the lines to it
  vtkSmartPointer<vtkCellArray> skeletonCells =
    vtkSmartPointer<vtkCellArray>::New();
  skeletonCells->InsertNextCell(skelLine);

  // Create a polydata to store everything in
  vtkSmartPointer<vtkPolyData> skelOut = vtkPolyData::SafeDownCast(
     outInfo1->Get(vtkDataObject::DATA_OBJECT()));

  // Add the points to the dataset
  skelOut->SetPoints(skelPoints);

  // Add the lines to the dataset
  skelOut->SetLines(skeletonCells);

  // convert the vtkpolydata to a pinocchio skeleton
  Skeleton pSkel;
  string jointName, lastJointName;
  for(size_t i = 0; i < skelPoints->GetNumberOfPoints(); ++i) {
    double pt[3] = {0,0,0};
    skelPoints->GetPoint(i, pt);
    jointName = "j" + to_string(i);
    if(i == 0) {
      pSkel.makeJoint(jointName, Vector3(pt[0],pt[1],pt[2]));
    } else {
      pSkel.makeJoint(jointName, Vector3(pt[0],pt[1],pt[2]), lastJointName);
    }
    lastJointName = jointName;
  }
  pSkel.initCompressed();
  vtkDebugMacro( << "PINOCCHIO: " << pinocchioDebug.str());
  pinocchioDebug.str("");
  pinocchioDebug.clear();
  // pSkel.setFoot("j0");
  /***********************
  * Run Pinocchio rigger *
  ***********************/
  PinocchioOutput pOut = autorig(pSkel, pMesh);
  vtkDebugMacro( << "PINOCCHIO: " << pinocchioDebug.str());

  // allocate vtk arrays for joint weights
  // std::vector<> (nJoints);
  vtkNew<vtkDoubleArray> const boneWeights;
  boneWeights->SetNumberOfComponents(boneCount);
  boneWeights->Allocate(nVerts);
  boneWeights->SetName("boneWeights");
  boneWeights->Fill(0.0);

  if(pOut.embedding.size() == 0) {
    vtkErrorMacro("problem with embedding");
  } else {
    for(auto vert : vtk2Pinocchio) {
      int vInd = vert.first;
      int pInd = vert.second;
      Vector<double, -1> v = pOut.attachment->getWeights(pInd);
      std::vector<double> ptWts(boneCount);
      for(int j = 0; j < boneCount; ++j) {
        ptWts[j] = v[j];
        // double d = floor(0.5 + v[j] * 10000.) / 10000.;
        // boneWeights->SetComponent(vInd, j, v[j]);
      }
      boneWeights->InsertTuple(vInd, ptWts.data());
    }
  }

  // apply refined skeleton to the output skeleton
  polydata->GetBounds(bnds);
  for(int i = 0; i < pOut.embedding.size(); ++i) {
    double x[3] = {
      pOut.embedding[i][0] * (bnds[1] - bnds[0]) + bnds[0],
      pOut.embedding[i][1] * (bnds[3] - bnds[2]) + bnds[2],
      pOut.embedding[i][2] * (bnds[5] - bnds[4]) + bnds[4]
    };
    skelOut->GetPoints()->SetPoint(i, x);
  }
  vtkDebugMacro(<< "number of elements in embedding: " << pOut.embedding.size());
  // don't forget to cleanup Pinocchio
  delete pOut.attachment;

  /********************************************
  * Immediately pass mesh through to output 0 *
  ********************************************/
  vtkPolyData *meshOut = vtkPolyData::SafeDownCast(
    outInfo0->Get(vtkDataObject::DATA_OBJECT()));
  meshOut->CopyStructure(polydata);

  // add weights data to the mesh as a field
  meshOut->GetPointData()->AddArray(boneWeights);

  /*********************************************
  * Copy just the input vertices into output 1 *
  *********************************************/
  // skelOut->Squeeze();
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
void vtkcreateSkeletonFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

////////////////////////////////////////////////////////////////////////////////
int vtkcreateSkeletonFilter::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
int vtkcreateSkeletonFilter::FillOutputPortInformation(int, vtkInformation *info)
{
  // Always returns a vtkPolyData
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
int vtkcreateSkeletonFilter::RequestInformation(
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

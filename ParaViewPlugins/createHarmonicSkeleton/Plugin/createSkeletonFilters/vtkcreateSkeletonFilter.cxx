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
#include "vtkLine.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkPolyLine.h"
#include "vtkBoundingBox.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkObjectFactory.h"

// CGAL includes
#include <boost/multiprecision/gmp.hpp>
#include <Eigen/Core>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>

#include "vtkcreateSkeletonFilter.h"

// bool runIGLHeatWeights(
//   const Eigen::MatrixXd V,
//   const Eigen::MatrixXi F,
//   const Eigen::MatrixXd C,
//   const Eigen::VectorXi P,
//   const Eigen::MatrixXi BE,
//   Eigen::MatrixXd& W);

vtkStandardNewMacro(vtkcreateSkeletonFilter);

// CGAL typedefs
typedef CGAL::Simple_cartesian<double>                      K;
typedef CGAL::Surface_mesh<K::Point_3>                      SM;
typedef boost::property_map<SM, CGAL::vertex_point_t>::type VPMap;
typedef boost::graph_traits<SM>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<SM>::edge_descriptor            edge_descriptor;
typedef boost::graph_traits<SM>::face_descriptor            face_descriptor;
typedef boost::graph_traits<SM>::halfedge_descriptor        halfedge_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<SM>       Skeletonization;
typedef Skeletonization::Skeleton                           Skeleton;
typedef Skeleton::vertex_descriptor                         Skeleton_vertex;
typedef Skeleton::edge_descriptor                           Skeleton_edge;

////////////////////////////////////////////////////////////////////////////////
vtkcreateSkeletonFilter::vtkcreateSkeletonFilter(){
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(2);
}

////////////////////////////////////////////////////////////////////////////////
vtkcreateSkeletonFilter::~vtkcreateSkeletonFilter(){}

////////////////////////////////////////////////////////////////////////////////
// Creates call CGAL algorithms to generate a skeleton of the input mesh
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

  /********************************************
  * Create a CGAL::SurfaceMesh (SM) from the input mesh *
  ********************************************/
  SM sm;
  VPMap vpmap = get(CGAL::vertex_point, sm);

  // Extract points
  std::vector<vertex_descriptor> vertex_map(nVerts);
  for (vtkIdType i=0; i<nVerts; ++i) {
    double coords[3];
    polydata->GetPoint(i, coords);
    vertex_descriptor v = add_vertex(sm);
    put(vpmap, v, K::Point_3(coords[0], coords[1], coords[2]));
    vertex_map[i] = v;
  }

  // Extract cells
  for (vtkIdType i = 0; i<nFaces; ++i) {
    vtkCell* cell_ptr = polydata->GetCell(i);
    vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
    std::vector<vertex_descriptor> vr(nb_vertices);
    for (vtkIdType k=0; k<nb_vertices; ++k)
      vr[k] = vertex_map[cell_ptr->GetPointId(k)];
    CGAL::Euler::add_face(vr, sm);
  }

  // cleanup any isolated vertices
  std::vector<vertex_descriptor> isolated_vertices;
  for(SM::vertex_iterator vit = sm.vertices_begin();
      vit != sm.vertices_end();
      ++vit)
  {
    if(sm.is_isolated(*vit))
      isolated_vertices.push_back(*vit);
  }

  for (std::size_t i=0; i < isolated_vertices.size(); ++i)
    sm.remove_vertex(isolated_vertices[i]);

  /*****************************
  * Create skeleton *
  *****************************/
  Skeleton skeleton;
  CGAL::extract_mean_curvature_flow_skeleton(sm, skeleton);
  int nSkelVerts = boost::num_vertices(skeleton);
  int nSkelEdges = boost::num_edges(skeleton);
  std::cout << "Number of verts in the skeleton: " << nSkelVerts << "\n";
  std::cout << "Number of edges in the skeleton: " << nSkelEdges << "\n";
  // split skeleton into polylines
  //split_graph_into_polylines

  /********************************************
  * Immediately pass mesh through to output 0 *
  ********************************************/
  vtkPolyData *meshOut = vtkPolyData::SafeDownCast(
    outInfo0->Get(vtkDataObject::DATA_OBJECT()));
  meshOut->DeepCopy(polydata);

  // add weights data to the mesh as a field
  // meshOut->GetPointData()->AddArray(boneWeights);
  // setup joint index array
  vtkNew<vtkDoubleArray> const jointArray;
  jointArray->SetNumberOfComponents(1);
  jointArray->SetName("joint");
  jointArray->Allocate(meshOut->GetPoints()->GetNumberOfPoints());
  /*********************************************
  * fill in the vtk skeleton*
  *********************************************/
  vtkPolyData *skelOut = vtkPolyData::SafeDownCast(
    outInfo1->Get(vtkDataObject::DATA_OBJECT()));
  vtkNew<vtkPoints>    const vtkSkelVerts;
  vtkNew<vtkCellArray> const vtkSkelEdges;
  vtkSkelVerts->Allocate(nSkelVerts);
  vtkSkelEdges->Allocate(nSkelEdges);

  std::vector<vtkIdType> Vids(num_vertices(sm));
  vtkIdType inum = 0;
  // read vertices from skeleton
  for(Skeleton_vertex v : CGAL::make_range(vertices(skeleton))) {
    vtkSkelVerts->InsertNextPoint(
      CGAL::to_double(skeleton[v].point.x()),
      CGAL::to_double(skeleton[v].point.y()),
      CGAL::to_double(skeleton[v].point.z()));
    Vids[v] = inum++;
    // loop over all the nodes in the mesh associated with this skeleton joint
    for(vertex_descriptor vd : skeleton[v].vertices) {
      jointArray->InsertTuple1(vd, Vids[v]);
    }
  }

  meshOut->GetPointData()->AddArray(jointArray);

  // read edges from skeleton
  for(Skeleton_edge e : CGAL::make_range(edges(skeleton))) {
    vtkSmartPointer<vtkLine> line =
      vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0,Vids[source(e, skeleton)]);
    line->GetPointIds()->SetId(1,Vids[target(e, skeleton)]);
    vtkSkelEdges->InsertNextCell(line);
  }
  // cast points into output skeleton
  skelOut->SetPoints(vtkSkelVerts);
  skelOut->SetLines(vtkSkelEdges);
  skelOut->Squeeze();

  // write skeleton to file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  std::string cache_path = std::string(CACHE_PATH) + std::string( "/skeleton.vtp");
  writer->SetFileName ( cache_path.c_str() );
  writer->SetInputData ( skelOut );
  writer->Write();
  
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
int vtkcreateSkeletonFilter::FillOutputPortInformation(int,vtkInformation *info)
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

    return 1;
}

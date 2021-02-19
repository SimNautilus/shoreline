#include <ostream>
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>
#include "vtkInformationVector.h"
#include "exportDisplacementFilter.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"

// Declare the plugin
vtkStandardNewMacro(exportDisplacementFilter);
////////////////////////////////////////////////////////////////////////////////
// Constructor
// Fills the number of input and output objects.
// Initializes the members that need it.
exportDisplacementFilter::exportDisplacementFilter() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
  this->DebugOn();
}

////////////////////////////////////////////////////////////////////////////////
// Gets the input
// Creates CGAL::Surface_mesh from vtkPolydata
// Calls the CGAL::isotropic_remeshing algorithm
// Fills the output vtkPolyData from the result.
int exportDisplacementFilter::RequestData(
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

  /**********************************
  * Pass the mesh through to the output *
  **********************************/
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->CopyStructure(polydata);
  // get the features array
  if(polydata->GetCellData()->HasArray("features")) {
    vtkDataArray* featureData = polydata->GetCellData()->GetArray("features");
    output->GetCellData()->AddArray(featureData);
  }

  // get a displacement field
  if(polydata->GetPointData()->HasArray("displacements")) {
    vtkDataArray* dispData = polydata->GetPointData()->GetArray("displacements");
    output->GetPointData()->AddArray(dispData);
  } else {
    vtkIdType nVerts = output->GetNumberOfPoints();
    vtkNew<vtkDoubleArray> const disps;
    disps->SetNumberOfComponents(3);
    disps->SetName("displacements");
    disps->Allocate(nVerts);
    disps->Fill(0.0);
    output->GetPointData()->AddArray(disps);
  }
  /***********************************
  * Write displacement field to file *
  ***********************************/
  if(output->GetPointData()->HasArray("displacements")) {
    std::string fpath = filePath;
    if(fpath != "") {
      // hopefully the user entered a valid directory path!
      // this will only work on Unix!
      std::string filename = fpath + ((fpath.back() == '/') ? "" : "/")
                           + "displacements.csv";
      std::ofstream file;
      file.open(filename);

      // write displacements to file
      vtkDataArray* disps = output->GetPointData()->GetArray("displacements");
      vtkIdType nVerts = output->GetNumberOfPoints();
      for(vtkIdType i = 0; i < nVerts; ++i) {
        double d[3];
        disps->GetTuple(i, d);
        file << std::setprecision(16) << d[0] << " " << d[1] << " " << d[2] << "\n";
      }

      file.close();
    }
  }
  endTime = std::chrono::system_clock::now();
  vtkDebugMacro(<< "exportDisplacement: " <<std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << " ms");

  return 1;
}

////////////////////////////////////////////////////////////////////////////////
void exportDisplacementFilter::PrintSelf(std::ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os,indent);
}

////////////////////////////////////////////////////////////////////////////////
int exportDisplacementFilter::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
int exportDisplacementFilter::FillOutputPortInformation(int, vtkInformation *info)
{
  // Always returns a vtkPolyData
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
int exportDisplacementFilter::RequestInformation(
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

/* @file MeshCatalyst.h
 * @brief Header file for the MeshCatalyst class
 * @details MeshCatalyst provides routines for sending data over Catalyst
 *          to a Paraview client
 */

 #pragma once

// Shoreline Incldues
#include "MeshCourier.h"

#ifdef BUILD_WITH_PARAVIEW
#include <vtkCPDataDescription.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
vtkCPProcessor* Processor = NULL; // static data
#endif

class MeshCatalyst : public MeshCourier {
public:
  //////////////////////////////////////////////////////////////////////////////
#ifdef BUILD_WITH_PARAVIEW
  // Catalyst Functionality
  /* @brief initialize the catalyst runtime environment */
  static void CatalystInit() {// int numScripts, char* scripts[]
    if(Processor == NULL)
    {
      Processor = vtkCPProcessor::New(); Processor->Initialize();
    }
    std::string script = "../pythonScripts/allinputsgridwriter.py";
    vtkCPPythonScriptPipeline* pipeline = vtkCPPythonScriptPipeline::New();
    pipeline->Initialize(script.c_str());
    Processor->AddPipeline(pipeline);
    pipeline->Delete();
    // scripts are passed in as command line arguments
    // for(int i=0;i<numScripts;i++) {
    //   vtkCPPythonScriptPipeline* pipeline = vtkCPPythonScriptPipeline::New();
    //   pipeline->Initialize(scripts[i]);
    //   Processor->AddPipeline(pipeline);
    //   pipeline->Delete();
    // }
  }

  /* @brief finalize the catalyst runtime environment */
  static void CatalystFinalize() {
    if(Processor) {
      Processor->Delete(); Processor = NULL;
    }
  }

  /* @brief The coprocessing step
     @details This function takes an STLMesh and passes it to catalyst
   */
  template <typename T>
  static void CatalystCoProcess( int timeStep, double time,
    std::map<T, std::array<double,3>>& verts, std::vector<std::array<T,3>>& ien)
  {
    vtkCPDataDescription* datadescription = vtkCPDataDescription::New();
    datadescription->AddInput("input");
    datadescription->SetTimeData(time, timeStep);
    if(Processor->RequestDataDescription(datadescription) != 0) {
      // Catalyst needs to output data
      // create the vtkPolyData object to pass to Catalyst
      vtkSmartPointer<vtkPolyData> polydata = STLMesh2VTKPolyData(verts, ien);
      datadescription->GetInputDescriptionByName("input")->SetGrid(polydata);
      polydata->Delete();

      // pass data to coprocessor
      Processor->CoProcess(datadescription);
    }
    datadescription->Delete();
  }
#endif
}

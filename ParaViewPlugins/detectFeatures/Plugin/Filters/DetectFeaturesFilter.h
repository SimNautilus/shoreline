#ifndef DetectFeaturesFilter_h
#define DetectFeaturesFilter_h
// Gives access to macros for communication with the UI
#include "vtkFiltersCoreModule.h"
#include "vtkGeometryFilter.h"

// Inherit from the desired filter
class DetectFeaturesFilter : public vtkGeometryFilter
{
public:
  // VTK requirements
  static DetectFeaturesFilter* New();
  vtkTypeMacro(DetectFeaturesFilter, vtkGeometryFilter);
  // Prints the values of the specific data
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Communicate with the UI
  vtkSetMacro(AngleInDegrees, double);
  vtkGetMacro(AngleInDegrees, double);
  vtkSetMacro(AngleInfo, double);
  vtkGetMacro(AngleInfo, double);

  // Pipeline functions:
  // Performs the isotropic remeshing algorithm and fills the output object here.
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *)override;
  // Specifies the type of the input objects
  int FillInputPortInformation(int, vtkInformation *info)override;
  // Specifies the type of the output object.
  int FillOutputPortInformation(int, vtkInformation *info)override;

protected:
  DetectFeaturesFilter();
  ~DetectFeaturesFilter(){}

  // Computes the bbox's diagonal length to set the default target edge length.
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
  // Data set by the UI and used by the algorithm
  double AngleInDegrees = 90;
  double AngleInfo = 90;

  // needed but not implemented
  DetectFeaturesFilter(const DetectFeaturesFilter&);
  void operator=(const DetectFeaturesFilter&);
};
#endif

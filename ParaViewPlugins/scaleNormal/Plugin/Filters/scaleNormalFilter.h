#ifndef scaleNormalFilter_h
#define scaleNormalFilter_h

  #include "scaleNormalFiltersModule.h" // for export macro

// Gives access to macros for communication with the UI
#include "vtkFiltersCoreModule.h"
#include "vtkGeometryFilter.h"

// Inherit from the desired filter
class scaleNormalFilter : public vtkGeometryFilter
{
public:
  // VTK requirements
  static scaleNormalFilter* New();
  vtkTypeMacro(scaleNormalFilter, vtkGeometryFilter);
  // Prints the values of the specific data
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Communicate with the UI

  vtkSetStringMacro(regionsToScale);
  vtkGetStringMacro(regionsToScale);

  vtkSetStringMacro(regionsToFix);
  vtkGetStringMacro(regionsToFix);

  vtkSetMacro(scaleValue, double);
  vtkGetMacro(scaleValue, double);

  // Pipeline functions:
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *)override;
  // Specifies the type of the input objects
  int FillInputPortInformation(int, vtkInformation *info)override;
  // Specifies the type of the output object.
  int FillOutputPortInformation(int, vtkInformation *info)override;

protected:
  scaleNormalFilter();
  ~scaleNormalFilter(){}

  // Computes the bbox's diagonal length to set the default target edge length.
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
  // Data set by the UI and used by the algorithm
  char*  regionsToScale = NULL;
  char*  regionsToFix = NULL;
  double scaleValue = 0;
  // needed but not implemented
  scaleNormalFilter(const scaleNormalFilter&);
  void operator=(const scaleNormalFilter&);
};
#endif

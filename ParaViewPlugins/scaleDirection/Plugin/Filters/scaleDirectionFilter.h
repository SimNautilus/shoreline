#ifndef scaleDirectionFilter_h
#define scaleDirectionFilter_h

#include "scaleDirectionFiltersModule.h" // for export macro

// Gives access to macros for communication with the UI
#include "vtkFiltersCoreModule.h"
#include "vtkGeometryFilter.h"

// Inherit from the desired filter
class scaleDirectionFilter : public vtkGeometryFilter
{
public:
  // VTK requirements
  static scaleDirectionFilter* New();
  vtkTypeMacro(scaleDirectionFilter, vtkGeometryFilter);
  // Prints the values of the specific data
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Communicate with the UI

  vtkSetStringMacro(regionsToScale);
  vtkGetStringMacro(regionsToScale);

  vtkSetStringMacro(regionsToFix);
  vtkGetStringMacro(regionsToFix);

  vtkSetVector3Macro(scaleVector, double);
  vtkGetVector3Macro(scaleVector, double);

  vtkSetMacro(HarmonicOrder, int);
  vtkGetMacro(HarmonicOrder, int);

  // Pipeline functions:
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *)override;
  // Specifies the type of the input objects
  int FillInputPortInformation(int, vtkInformation *info)override;
  // Specifies the type of the output object.
  int FillOutputPortInformation(int, vtkInformation *info)override;

protected:
  scaleDirectionFilter();
  ~scaleDirectionFilter(){}

  // Computes the bbox's diagonal length to set the default target edge length.
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
  // Data set by the UI and used by the algorithm
  char*  regionsToScale = NULL;
  char*  regionsToFix = NULL;
  double scaleVector[3];
  int HarmonicOrder = 2;

  // needed but not implemented
  scaleDirectionFilter(const scaleDirectionFilter&);
  void operator=(const scaleDirectionFilter&);
};
#endif

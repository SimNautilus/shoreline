#ifndef exportDisplacementFilter_h
#define exportDisplacementFilter_h
// Gives access to macros for communication with the UI
#include "vtkFiltersCoreModule.h"
#include "vtkGeometryFilter.h"

// Inherit from the desired filter
class exportDisplacementFilter : public vtkGeometryFilter
{
public:
  // VTK requirements
  static exportDisplacementFilter* New();
  vtkTypeMacro(exportDisplacementFilter, vtkGeometryFilter);
  // Prints the values of the specific data
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Communicate with the UI
  vtkSetStringMacro(filePath);
  vtkGetStringMacro(filePath);

  // Pipeline functions:
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *)override;
  // Specifies the type of the input objects
  int FillInputPortInformation(int, vtkInformation *info)override;
  // Specifies the type of the output object.
  int FillOutputPortInformation(int, vtkInformation *info)override;

protected:
  exportDisplacementFilter();
  ~exportDisplacementFilter(){}

  // Computes the bbox's diagonal length to set the default target edge length.
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
  // Data set by the UI and used by the algorithm
  char* filePath = NULL;

  // needed but not implemented
  exportDisplacementFilter(const exportDisplacementFilter&);
  void operator=(const exportDisplacementFilter&);
};
#endif

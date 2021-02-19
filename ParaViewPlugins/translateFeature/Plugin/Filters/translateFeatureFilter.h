#ifndef translateFeatureFilter_h
#define translateFeatureFilter_h

#include "translateFeatureFiltersModule.h" // for export macro

// Gives access to macros for communication with the UI
#include "vtkFiltersCoreModule.h"
#include "vtkGeometryFilter.h"

// Inherit from the desired filter
class translateFeatureFilter : public vtkGeometryFilter
{
public:
  // VTK requirements
  static translateFeatureFilter* New();
  vtkTypeMacro(translateFeatureFilter, vtkGeometryFilter);
  // Prints the values of the specific data
  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Communicate with the UI

  // vtkSetMacro(regionToMove, int);
  // vtkGetMacro(regionToMove, int);

  vtkSetMacro(HarmonicOrder, int);
  vtkGetMacro(HarmonicOrder, int);

  vtkSetStringMacro(regionsToTranslate);
  vtkGetStringMacro(regionsToTranslate);

  vtkSetStringMacro(regionsToFix);
  vtkGetStringMacro(regionsToFix);

  // vtkSetMacro(regionToFix, int);
  // vtkGetMacro(regionToFix, int);
  vtkSetVector3Macro(moveVector, double);
  vtkGetVector3Macro(moveVector, double);
  //
  // vtkSetVector2Macro(regionToFix, int);
  // vtkGetVector2Macro(regionToFix, int);
  // Pipeline functions:
  // Performs the isotropic remeshing algorithm and fills the output object here.
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *)override;
  // Specifies the type of the input objects
  int FillInputPortInformation(int, vtkInformation *info)override;
  // Specifies the type of the output object.
  int FillOutputPortInformation(int, vtkInformation *info)override;

protected:
  translateFeatureFilter();
  ~translateFeatureFilter(){}

  // Computes the bbox's diagonal length to set the default target edge length.
  int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
  // Data set by the UI and used by the algorithm
  double Length;
  double LengthInfo;
  int MainIterations;


  char*  regionsToTranslate = NULL;
  char*  regionsToFix = NULL;
  // int regionToMove;
  // int regionToFix[2];
  double moveVector[3];
  int HarmonicOrder = 2;

  // needed but not implemented
  translateFeatureFilter(const translateFeatureFilter&);
  void operator=(const translateFeatureFilter&);
};
#endif

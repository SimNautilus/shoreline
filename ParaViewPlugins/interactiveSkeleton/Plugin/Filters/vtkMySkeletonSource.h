/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMySkeletonSource.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkMySkeletonSource
 * @brief   create a skeleton
 */

#ifndef vtkMySkeletonSource_h
#define vtkMySkeletonSource_h

#include "interactiveSkeletonFiltersModule.h" // for export macro
#include "vtkPolyDataAlgorithm.h"

class vtkPolyData;
class vtkCellArray;
class vtkPoints;

class INTERACTIVESKELETONFILTERS_EXPORT vtkMySkeletonSource :
  public vtkPolyDataAlgorithm
{
public:
  static vtkMySkeletonSource* New();
  vtkTypeMacro(vtkMySkeletonSource, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkPoints* GetPoints(){return this->Skeleton->GetPoints();}
  void GetPoints(vtkPoints* points);
  void GetLines(vtkCellArray* lines);

  int GetNumberOfPoints();
  int GetNumberOfLines();

  vtkGetObjectMacro(Skeleton, vtkPolyData);
  void SetSkeleton(vtkPolyData* skl);
protected:
  vtkMySkeletonSource();
  ~vtkMySkeletonSource() override = default;

  int RequestData(
    vtkInformation*,
    vtkInformationVector**,
    vtkInformationVector*
  ) override;

  vtkSmartPointer<vtkPolyData> Skeleton;

private:
  vtkMySkeletonSource(const vtkMySkeletonSource&) = delete;
  void operator=(const vtkMySkeletonSource&) = delete;
};

#endif

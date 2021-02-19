/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMySkeletonSource.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMySkeletonSource.h"

#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <cmath>

vtkStandardNewMacro(vtkMySkeletonSource);

//------------------------------------------------------------------------------
vtkMySkeletonSource::vtkMySkeletonSource()
{
  this->SetNumberOfInputPorts(0);
}

//------------------------------------------------------------------------------
int vtkMySkeletonSource::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  // get the info object
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the output
  vtkPolyData* output =
    vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->CopyStructure(this->Skeleton);

  return 1;
}

//------------------------------------------------------------------------------
void vtkMySkeletonSource::SetSkeleton(vtkPolyData* skl)
{
  if(!this->Skeleton)
  {
    this->Skeleton = vtkPolyData::New();
  }
  if(skl)
  {
    this->Skeleton->CopyStructure(skl);
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonSource::GetPoints(vtkPoints* points)
{
  this->Update();
  if(this->Skeleton)
  {
    this->Skeleton->GetPoints()->ShallowCopy(points);
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonSource::GetLines(vtkCellArray* lines)
{
  this->Update();
  if(this->Skeleton)
  {
    this->Skeleton->GetLines()->ShallowCopy(lines);
  }
}
//------------------------------------------------------------------------------
int vtkMySkeletonSource::GetNumberOfPoints()
{
  this->Update();
  return this->Skeleton->GetNumberOfPoints();
}

//------------------------------------------------------------------------------
int vtkMySkeletonSource::GetNumberOfLines()
{
  this->Update();
  return this->Skeleton->GetNumberOfLines();
}

//------------------------------------------------------------------------------
void vtkMySkeletonSource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

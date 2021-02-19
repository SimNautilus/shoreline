/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMySkeletonWidget.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkMySkeletonWidget
 * @brief   widget for vtkMySkeletonRepresentation.
 *
 * vtkMySkeletonWidget is the vtkAbstractWidget subclass for
 * vtkMySkeletonRepresentation which manages the interactions with
 * vtkMySkeletonRepresentation. This is based on vtkMySkeletonWidget.
 *
 * This widget allows the creation of a polyline interactively by adding or removing points
 * based on mouse position and a modifier key.
 *
 * - ctrl+click inserts a new point on the selected line
 * - shift+click deletes the selected point
 * - alt+click adds a new point anywhere depending on last selected point.
 *   If the first point is selected, the new point is added at the beginning,
 *   else it is added at the end.
 *
 * @sa
 * vtkMySkeletonRepresentation, vtkMySkeletonWidget
 */

#ifndef vtkMySkeletonWidget_h
#define vtkMySkeletonWidget_h

#include "vtkAbstractWidget.h"
#include "interactiveSkeletonFiltersModule.h" // for export macro

class vtkMySkeletonRepresentation;

class INTERACTIVESKELETONFILTERS_EXPORT vtkMySkeletonWidget : public vtkAbstractWidget
{
public:
  static vtkMySkeletonWidget* New();
  vtkTypeMacro(vtkMySkeletonWidget, vtkAbstractWidget);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Specify an instance of vtkMyWidgetRepresentation used to represent this
   * widget in the scene. Note that the representation is a subclass of
   * vtkProp so it can be added to the renderer independent of the widget.
   */
  void SetRepresentation(vtkMySkeletonRepresentation* r)
  {
    this->Superclass::SetWidgetRepresentation(reinterpret_cast<vtkWidgetRepresentation*>(r));
  }

  /**
   * Create the default widget representation if one is not set. By default,
   * this is an instance of the vtkMySkeletonRepresentation class.
   */
  void CreateDefaultRepresentation() override;

  /**
   * Override superclasses' SetEnabled() method because the line
   * widget must enable its internal handle widgets.
   */
  void SetEnabled(int enabling) override;

protected:
  vtkMySkeletonWidget();
  ~vtkMySkeletonWidget() override;

  int WidgetState;
  enum _WidgetState
  {
    Start = 0,
    Active
  };

  // These methods handle events
  static void SelectAction(vtkAbstractWidget*);
  static void EndSelectAction(vtkAbstractWidget*);
  static void TranslateAction(vtkAbstractWidget*);
  static void ScaleAction(vtkAbstractWidget*);
  static void MoveAction(vtkAbstractWidget*);

  vtkCallbackCommand* KeyEventCallbackCommand;
  static void ProcessKeyEvents(vtkObject*, unsigned long, void*, void*);

private:
  vtkMySkeletonWidget(const vtkMySkeletonWidget&) = delete;
  void operator=(const vtkMySkeletonWidget&) = delete;
};

#endif

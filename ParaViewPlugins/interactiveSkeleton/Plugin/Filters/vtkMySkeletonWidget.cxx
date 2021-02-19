/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMySkeletonWidget.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMySkeletonWidget.h"

#include "vtkCallbackCommand.h"
#include "vtkCommand.h"
#include "vtkEvent.h"
#include "vtkObjectFactory.h"
#include "vtkMySkeletonRepresentation.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkWidgetCallbackMapper.h"
#include "vtkWidgetEvent.h"
#include "vtkWidgetEventTranslator.h"

vtkStandardNewMacro(vtkMySkeletonWidget);
//------------------------------------------------------------------------------
vtkMySkeletonWidget::vtkMySkeletonWidget()
{
  this->WidgetState = vtkMySkeletonWidget::Start;
  this->ManagesCursor = 1;

  // Define widget events
  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonPressEvent, vtkWidgetEvent::Select,
    this, vtkMySkeletonWidget::SelectAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonReleaseEvent,
    vtkWidgetEvent::EndSelect, this, vtkMySkeletonWidget::EndSelectAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::MiddleButtonPressEvent,
    vtkWidgetEvent::Translate, this, vtkMySkeletonWidget::TranslateAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::MiddleButtonReleaseEvent,
    vtkWidgetEvent::EndTranslate, this, vtkMySkeletonWidget::EndSelectAction);
  this->CallbackMapper->SetCallbackMethod(
    vtkCommand::RightButtonPressEvent, vtkWidgetEvent::Scale, this, vtkMySkeletonWidget::ScaleAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::RightButtonReleaseEvent,
    vtkWidgetEvent::EndScale, this, vtkMySkeletonWidget::EndSelectAction);
  this->CallbackMapper->SetCallbackMethod(
    vtkCommand::MouseMoveEvent, vtkWidgetEvent::Move, this, vtkMySkeletonWidget::MoveAction);

  this->KeyEventCallbackCommand = vtkCallbackCommand::New();
  this->KeyEventCallbackCommand->SetClientData(this);
  this->KeyEventCallbackCommand->SetCallback(vtkMySkeletonWidget::ProcessKeyEvents);
}

//------------------------------------------------------------------------------
vtkMySkeletonWidget::~vtkMySkeletonWidget()
{
  this->KeyEventCallbackCommand->Delete();
}

//------------------------------------------------------------------------------
void vtkMySkeletonWidget::SetEnabled(int enabling)
{
  int enabled = this->Enabled;

  // We do this step first because it sets the CurrentRenderer
  this->Superclass::SetEnabled(enabling);

  // We defer enabling the handles until the selection process begins
  if (enabling && !enabled)
  {
    if (this->Parent)
    {
      this->Parent->AddObserver(
        vtkCommand::KeyPressEvent, this->KeyEventCallbackCommand, this->Priority);
      this->Parent->AddObserver(
        vtkCommand::KeyReleaseEvent, this->KeyEventCallbackCommand, this->Priority);
    }
    else
    {
      this->Interactor->AddObserver(
        vtkCommand::KeyPressEvent, this->KeyEventCallbackCommand, this->Priority);
      this->Interactor->AddObserver(
        vtkCommand::KeyReleaseEvent, this->KeyEventCallbackCommand, this->Priority);
    }
  }
  else if (!enabling && enabled)
  {
    if (this->Parent)
    {
      this->Parent->RemoveObserver(this->KeyEventCallbackCommand);
    }
    else
    {
      this->Interactor->RemoveObserver(this->KeyEventCallbackCommand);
    }
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonWidget::SelectAction(vtkAbstractWidget* w)
{
  // We are in a static method, cast to ourself
  vtkMySkeletonWidget* self = vtkMySkeletonWidget::SafeDownCast(w);

  // Get the event position
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  if (!self->CurrentRenderer || !self->CurrentRenderer->IsInViewport(X, Y))
  {
    self->WidgetState = vtkMySkeletonWidget::Start;
    return;
  }

  // Begin the widget interaction which has the side effect of setting the
  // interaction state.
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  self->WidgetRep->StartWidgetInteraction(e);
  int interactionState = self->WidgetRep->GetInteractionState();
  if (interactionState == vtkMySkeletonRepresentation::Outside && !self->Interactor->GetAltKey())
  {
    return;
  }
  // We are definitely selected
  self->WidgetState = vtkMySkeletonWidget::Active;
  self->GrabFocus(self->EventCallbackCommand);

  // if (self->Interactor->GetAltKey())
  // {
  //   // push point.
  //   reinterpret_cast<vtkMySkeletonRepresentation*>(self->WidgetRep)
  //     ->SetInteractionState(vtkMySkeletonRepresentation::Pushing);
  // }
  // else if (interactionState == vtkMySkeletonRepresentation::OnLine &&
  //   self->Interactor->GetControlKey())
  // {
  //   // insert point.
  //   reinterpret_cast<vtkMySkeletonRepresentation*>(self->WidgetRep)
  //     ->SetInteractionState(vtkMySkeletonRepresentation::Inserting);
  // }
  // else if (interactionState == vtkMySkeletonRepresentation::OnHandle &&
  //   self->Interactor->GetShiftKey())
  // {
  //   // remove point.
  //   reinterpret_cast<vtkMySkeletonRepresentation*>(self->WidgetRep)
  //     ->SetInteractionState(vtkMySkeletonRepresentation::Erasing);
  // }

  reinterpret_cast<vtkMySkeletonRepresentation*>(self->WidgetRep)
    ->SetInteractionState(vtkMySkeletonRepresentation::Moving);

  // start the interaction
  self->EventCallbackCommand->SetAbortFlag(1);
  self->StartInteraction();
  self->InvokeEvent(vtkCommand::StartInteractionEvent, nullptr);
  self->Render();
}

//------------------------------------------------------------------------------
void vtkMySkeletonWidget::TranslateAction(vtkAbstractWidget* w)
{
  // Not sure this should be any different than SelectAction
  vtkMySkeletonWidget::SelectAction(w);
}

//------------------------------------------------------------------------------
void vtkMySkeletonWidget::ScaleAction(vtkAbstractWidget* w)
{
  // We are in a static method, cast to ourself
  vtkMySkeletonWidget* self = reinterpret_cast<vtkMySkeletonWidget*>(w);

  // Get the event position
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  if (!self->CurrentRenderer || !self->CurrentRenderer->IsInViewport(X, Y))
  {
    self->WidgetState = vtkMySkeletonWidget::Start;
    return;
  }

  // Begin the widget interaction which has the side effect of setting the
  // interaction state.
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  self->WidgetRep->StartWidgetInteraction(e);
  int interactionState = self->WidgetRep->GetInteractionState();
  if (interactionState == vtkMySkeletonRepresentation::Outside)
  {
    return;
  }

  // We are definitely selected
  self->WidgetState = vtkMySkeletonWidget::Active;
  self->GrabFocus(self->EventCallbackCommand);
  // Scale
  reinterpret_cast<vtkMySkeletonRepresentation*>(self->WidgetRep)
    ->SetInteractionState(vtkMySkeletonRepresentation::Scaling);

  // start the interaction
  self->EventCallbackCommand->SetAbortFlag(1);
  self->StartInteraction();
  self->InvokeEvent(vtkCommand::StartInteractionEvent, nullptr);
  self->Render();
}

//------------------------------------------------------------------------------
void vtkMySkeletonWidget::MoveAction(vtkAbstractWidget* w)
{
  vtkMySkeletonWidget* self = reinterpret_cast<vtkMySkeletonWidget*>(w);

  // See whether we're active
  if (self->WidgetState == vtkMySkeletonWidget::Start)
  {
    return;
  }

  // compute some info we need for all cases
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, adjust the representation
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  self->WidgetRep->WidgetInteraction(e);

  // moving something
  self->EventCallbackCommand->SetAbortFlag(1);
  self->InvokeEvent(vtkCommand::InteractionEvent, nullptr);
  self->Render();
}

//------------------------------------------------------------------------------
void vtkMySkeletonWidget::EndSelectAction(vtkAbstractWidget* w)
{
  vtkMySkeletonWidget* self = reinterpret_cast<vtkMySkeletonWidget*>(w);
  if (self->WidgetState == vtkMySkeletonWidget::Start)
  {
    return;
  }

  // compute some info we need for all cases
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, adjust the representation
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);

  self->WidgetRep->EndWidgetInteraction(e);

  // EndWidgetInteraction for this widget can modify/add/remove points
  // Make sure the representation is updated
  self->InvokeEvent(vtkCommand::InteractionEvent, nullptr);

  // Return state to not active
  self->WidgetState = vtkMySkeletonWidget::Start;
  reinterpret_cast<vtkMySkeletonRepresentation*>(self->WidgetRep)
    ->SetInteractionState(vtkMySkeletonRepresentation::Outside);
  self->ReleaseFocus();

  self->EventCallbackCommand->SetAbortFlag(1);
  self->EndInteraction();
  self->InvokeEvent(vtkCommand::EndInteractionEvent, nullptr);
  self->Render();
}

//------------------------------------------------------------------------------
void vtkMySkeletonWidget::ProcessKeyEvents(vtkObject*, unsigned long event, void* clientdata, void*)
{
  vtkMySkeletonWidget* self = static_cast<vtkMySkeletonWidget*>(clientdata);
  vtkRenderWindowInteractor* iren = self->GetInteractor();
  vtkMySkeletonRepresentation* rep = vtkMySkeletonRepresentation::SafeDownCast(self->WidgetRep);
  switch (event)
  {
    case vtkCommand::KeyPressEvent:
      switch (iren->GetKeyCode())
      {
        case 'x':
        case 'X':
          rep->SetXTranslationAxisOn();
          break;
        case 'y':
        case 'Y':
          rep->SetYTranslationAxisOn();
          break;
        case 'z':
        case 'Z':
          rep->SetZTranslationAxisOn();
          break;
        default:
          break;
      }
      break;
    case vtkCommand::KeyReleaseEvent:
      switch (iren->GetKeyCode())
      {
        case 'x':
        case 'X':
        case 'y':
        case 'Y':
        case 'z':
        case 'Z':
          rep->SetTranslationAxisOff();
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonWidget::CreateDefaultRepresentation()
{
  if (!this->WidgetRep)
  {
    this->WidgetRep = vtkMySkeletonRepresentation::New();
  }
}

//------------------------------------------------------------------------------
void vtkMySkeletonWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

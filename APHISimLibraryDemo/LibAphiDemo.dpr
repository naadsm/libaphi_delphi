program LibAphiDemo;

(*
LibAphiDemo.dpr
----------------
Begin: 2008/12/11
Last revision: $Date: 2010-05-20 17:32:54 $ $Author: areeves $
Version: $Revision: 1.9 $
Project: APHI Delphi Library for Simulation Modeling: Demo application
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2008 - 2010 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
*)

{$INCLUDE Defs.inc}


(*
The Delphi interface overwrites user comments in the
"uses" section.  Remember to copy this section back from a good version!!
(Also remember that comments in braces in the "uses" section have a
special purpose, so don't mess with them.)
*)

uses
  //*************************************************************************************
  // This block can be used as a template for any new simulation modeling application.
  //*************************************************************************************
  // Built-in Delphi units
  //----------------------
  Windows,
  Forms,
  SysUtils,
  ComObj,
  Controls,

  // APHI General Purpose Delphi Library
  //------------------------------------
  AppLog,
  ARMath,
  ARMathAdvanced,
  BasicGIS,
  CmdLine,
  ControlUtils,
  CStringList,
  CsvParser,
  DebugWindow,
  FunctionPointers,
  I88n,
  ImageResources,
  IniHandler,
  MyDelphiArrayUtils,
  myDialogs,
  MyGraphicsUtils,
  MyStrUtils,
  Points,
  RegExpDefs,
  RegExpr,
  Resources,
  RoundToXReplacement_3c,
  SqlClasses,
  StringSuperList,
  UnicodeDev,
  WindowsUtils,
  ZipFunctions,

  // APHI Delphi Library for Simulation Modeling
  //--------------------------------------------
  // NOTE: an application-specific unit called FunctionEnums.pas must exist
  // to define the functions used in a model (see below).  This unit can be
  // generated from the file called FunctionEnums.pas.template, which is
  // included with the LibAPHI source code.
  AphiRng,
  ChartFunction,
  RelFunction,
  ProbDensityFunctions,
  FunctionDictionary,
  Models,
  SimInput,
  {$IFDEF DATABASE_ENABLED}
  ModelDatabase,
  {$ENDIF}
  
  // Graphical components in the APHI General Purpose Delphi Library
  //----------------------------------------------------------------
  FrameAcceptCancel in 'general_purpose_gui\FrameAcceptCancel.pas' {FrameAcceptCancel: TFrame},
  FrameChartBase in 'general_purpose_gui\FrameChartBase.pas' {FrameChartBase: TFrame},
  DialogLongMessage in 'general_purpose_gui\DialogLongMessage.pas' {DialogLongMessage},

  // Graphical components in the APHI Delphi Library for Simulation Modeling
  //------------------------------------------------------------------------
  FrameArrayHistogram in 'libaphi_delphi_gui\FrameArrayHistogram.pas' {FrameArrayHistogram: TFrame},

  // PDF/relational function editing, part of the APHI Delphi Library for Simulation Modeling
  //-----------------------------------------------------------------------------------------
  FrameFunctionParams2 in 'libaphi_delphi_gui\function_editor\FrameFunctionParams2.pas' {FrameFunctionParams2: TFrame},
  FrameChartPointsEditor in 'libaphi_delphi_gui\function_editor\FrameChartPointsEditor.pas' {FrameChartPointsEditor: TFrame},
  FormFunctionEditor in 'libaphi_delphi_gui\function_editor\FormFunctionEditor.pas' {FormFunctionEditor},
  FramePointEditorGrid in 'libaphi_delphi_gui\function_editor\FramePointEditorGrid.pas' {FramePointEditorGrid: TFrame},
  FrameFunctionEditor in 'libaphi_delphi_gui\function_editor\FrameFunctionEditor.pas' {FrameFunctionEditor: TFrame},
  //*************************************************************************************


  // Application-specific units
  //---------------------------  
  FunctionEnums in 'FunctionEnums.pas', // See the note above
  DemoModel in 'DemoModel.pas',
  DemoSimInput in 'DemoSimInput.pas',

  FormMain in 'FormMain.pas' {FormMain},
  FormAboutAphiDemo in 'FormAboutAphiDemo.pas' {FormAboutAphiDemo}
;

{$R *.res}

begin // MAIN
  {$IFDEF DEBUG}
    setDEBUGGING( true );
  {$ELSE}
    setDEBUGGING( false );
  {$ENDIF}
  
  Application.Initialize();
  Application.CreateForm( TFormMain, frmMain );
  Application.Run();
end.

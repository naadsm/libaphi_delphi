program AphiRngTest;

(*
AphiRngTest.dpr
---------------
Begin: 2008/11/10
Last revision: $Date: 2009-02-23 20:38:02 $ $Author: areeves $
Version: $Revision: 1.3 $
Project: APHI Delphi Library for Simulation Modeling: Demo for random number generation
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2008 - 2009 Animal Population Health Institute, Colorado State University

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
  // Standard Delphi units
  //----------------------
  Forms,
  Classes,
  Contnrs,
  SysUtils,

  // APHI General Purpose Delphi Library
  //------------------------------------
  MyStrUtils,
  DebugWindow {DebugWindow},

  // APHI Delphi Library for Simulation Modeling
  //--------------------------------------------
  AphiRng,
  ChartFunction,
  ProbDensityFunctions,

  // Application-specific units
  //---------------------------
  FunctionEnums in 'FunctionEnums.pas',  
  FormMain in 'FormMain.pas' {FormMain}
;

  
{$R *.res}

begin
  setDEBUGGING( true );

  Application.Initialize;
  Application.CreateForm( TFormMain, frmMain );
  Application.Run;
end.


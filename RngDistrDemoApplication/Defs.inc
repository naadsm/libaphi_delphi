(*
Defs.inc
--------
Begin: 2006/03/01
Last revision: $Date: 2008-12-11 21:27:59 $ $Author: areeves $
Version: $Revision: 1.2 $
Project: APHI Delphi Library for Simulation Modeling: Demo for random number generation
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2006 - 2008 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
*)


// The following symbol must be defined for the PDF classes to have access to an SQL database.
//{$DEFINE DATABASE_ENABLED}

// The following symbol must be defined for the PDF classes to generate calculations, instead
// of simply storing parameters. (For the test application to work, it should ALWAYS be defined.)
{$DEFINE USE_GSL}

// Uncomment these two lines to compile as a console application.
//{$APPTYPE CONSOLE}
//{$DEFINE CONSOLEAPP}
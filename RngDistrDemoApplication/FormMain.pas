unit FormMain;

(*
FormMain.pas/dfm
----------------
Begin: 2008/11/11
Last revision: $Date: 2009-04-06 17:40:47 $ $Author: dschwick $
Version: $Revision: 1.2 $
Project: APHI Delphi Library for Simulation Modeling: Demo for random number generation
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2008 - 2009 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
*)

interface

  uses
    // Standard Delphi units
    Windows,
    Messages,
    SysUtils,
    Variants,
    Classes,
    Graphics,
    Controls,
    Forms,
    Dialogs,
    StdCtrls, 
    ExtCtrls
  ;

  type TFormMain = class( TForm )
      lblDescr: TLabel;
      lblUrl: TLabel;
      btnTest: TButton;

      Panel1: TPanel;
      Panel2: TPanel;
      Panel3: TPanel;
      Panel4: TPanel;
      Panel5: TPanel;
      Panel6: TPanel;

      procedure btnTestClick(Sender: TObject);
      procedure lblUrlClick(Sender: TObject);

    protected
      procedure hideBevels( ctrl: TWinControl );
      
      procedure testPiecewisePdf();

    public
      constructor create( AOwner: TComponent ); override;
      destructor destroy(); override;

    end
  ;

  var
    frmMain: TFormMain;

implementation

{$R *.dfm}

  uses
    // Standard Delphi units
    ShellAPI,
  
    // APHI General Purpose Delphi Library
    DebugWindow,
    MyStrUtils,
    Points,

    // APHI Library for Simulation Modeling
    AphiRng,
    ChartFunction,
    ProbDensityFunctions
  ;

  constructor TFormMain.create( AOwner: TComponent );
    begin
      inherited create( AOwner );

      hideBevels( self );

      lblDescr.Caption :=
        'This application tests the SPRNG-based'
        + ' random number generator'
        + ' and then generates a set of random variates from'
        + ' each of the probability density function supported'
        + ' in the APHI Library for Simulation Modeling.'
        + endl + endl
        + 'Click the button below, and random variates will be displayed'
        + ' in the debugging window.'
        + endl + endl
        + 'For more information, please see the web page at the URL displayed below.'
      ;

      if( libAphiRngLoaded() ) then
        begin
          dbcout2( 'Random number generator library is ready.' );
          dbcout2( 'RNG was initialized with seed value of ' + intToStr( rngSeed() ) );
        end
      else
        begin
          dbcout2( 'Random number generator did not load.' );
          btnTest.Enabled := false;
        end
      ;
    end
  ;


  destructor TFormMain.destroy();
    begin
      inherited destroy();
    end
  ;


  procedure TFormMain.hideBevels( ctrl: TWinControl );
    var
      i: integer;
    begin
      for i := 0 to ctrl.ControlCount - 1 do
        begin
          if( ctrl.Controls[i] is TPanel ) then
            begin
              ( ctrl.Controls[i] as TPanel ).BevelOuter := bvNone;
              hideBevels( ctrl.Controls[i] as TPanel );
            end
          ;
        end
      ;
    end
  ;


  procedure TFormMain.lblUrlClick(Sender: TObject);
    begin
      lblUrl.Font.Color := clRed;
      lblUrl.Repaint();

      ShellExecute(
        Application.Handle,
        PChar( 'open' ),
        PChar( lblUrl.caption ),
        PChar( 0 ),
        nil,
        SW_NORMAL
      );

      lblUrl.Font.Color := clPurple;
    end
  ;


  procedure TFormMain.testPiecewisePdf();
    var
      pdf: TPdf;
      points: RPointArray;
    begin
      setLength( points, 6 );
      points[0].x := 0;  points[0].y := 0;
      points[1].x := 10; points[1].y := 0.038468;
      points[2].x := 20; points[2].y := 0.00097;
      points[3].x := 30; points[3].y := 0.00097;
      points[4].x := 45; points[4].y := 0.039567;
      points[5].x := 60; points[5].y := 0;

      pdf := TPdfPiecewise.create( points, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );
    end
  ;


  procedure TFormMain.btnTestClick(Sender: TObject);
    var
      i, j: integer;
      pdf: TPdf;
    begin
      dbClear();

      // Some simple tests
      //------------------
      dbcout2( 'Random doubles between 0 and 1:' );
      for i := 0 to 9 do
        dbcout2( rngRand() )
      ;

      dbcout2( '' );

      dbcout2( 'Random integers between 1 and 10: ' );
      for i := 0 to 9 do
        dbcout2( rngRandInt( 1,11 ) ) // Note the use of maxPlusOne!
      ;

      dbcout2( '' );

      // Testing the RNG seed
      //---------------------
      dbcout2( 'Random doubles between 0 and 1 after setting the RNG seed to 527:' );
      for j := 0 to 2 do
        begin
          setRngSeed( 527 );
          dbcout2( 'Seed reset to ' + intToStr( rngSeed() ) );
          for i := 0 to 9 do
            dbcout2( rngRand() )
          ;
          dbcout2( '' );
        end
      ;

      setRngSeed();
      dbcout2( 'Seed reset to ' + intToStr( rngSeed() ) );
      dbcout2( '' );

      // Testing the probability density functions
      //------------------------------------------
      dbcout2( 'Random variates from probability density functions:' );
      dbcout2( '' );

      pdf := TPdfBeta.create( 13, 3, 2, 50, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfBetaPert.create( 1, 7, 10, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfExponential.create( 3, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfGamma.create( 5, 5, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfGaussian.create( 10, 2, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfLogistic.create( 30, 3, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfLogLogistic.create( 3, 10, 10, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfLognormal.create( 10, 2, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfPareto.create( 2, 2, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfPearson5.create( 10, 30, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      testPiecewisePdf();

      pdf := TPdfPoint.create( 2.5, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfTriangular.create( 1, 7, 10, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfUniform.create( 10, 20, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      pdf := TPdfWeibull.create( 5, 5, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      dbcout2( 'Attempting to generate variates from an invalid PDF:' );
      pdf := TPdfTriangular.create( 7, 1, 10, UUnknown );
      pdf.randTest();
      freeAndNil( pdf );
      dbcout2( '' );

      dbcout2( 'Done!' );
    end
  ;


end.
 
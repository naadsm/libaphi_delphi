unit FormMain;

(*
FormMain.pas/dfm
----------------
Begin: 2008/12/11
Last revision: $Date: 2010-05-20 17:32:54 $ $Author: areeves $
Version: $Revision: 1.8 $
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
    ExtCtrls,
    ComCtrls,
    ToolWin,
    Grids,
    ActnMan,
    ActnCtrls,
    ActnMenus,
    ActnList,
    XPStyleActnCtrls,
    ImgList,

    // APHI Delphi controls
    REEdit,

    // APHI General Purpose Delphi Library
    FrameChartBase,

    // APHI Delphi Library for Simulation Modeling
    FunctionDictionary,
    FrameFunctionEditor,
    FrameArrayHistogram,

    // Application-specific units
    DemoSimInput,
    DemoModel
  ;

  type TFormMain = class(TForm)
      pnlPdfs: TPanel;
      smcPdf1: TFrameFunctionEditor;
      smcPdf2: TFrameFunctionEditor;
      imgPdf1: TImage;
      lblPdf1: TLabel;
      imgPdf2: TImage;
      lblPdf2: TLabel;
      mnuMain: TActionMainMenuBar;
      pnlFooter: TPanel;
      pnlButtons: TPanel;
      btnRun: TButton;
      btnTest: TButton;
      ActionManager1: TActionManager;
      acnFile: TAction;
      acnQuit: TAction;
      ImageList1: TImageList;
      acnRun: TAction;
      lblIterations: TLabel;
      reIterations: TREEdit;
      pnlOutput: TPanel;
      fraOutputHistogram: TFrameArrayHistogram;
      acnEdit: TAction;
      acnSaveHistogram: TAction;
      acnCopyHistogram: TAction;
      acnPrintHistogram: TAction;
      dlgPrint: TPrintDialog;
      dlgSaveWMF: TSaveDialog;
      acnHelp: TAction;
      acnAbout: TAction;

      // Construction/initialization
      procedure FormCreate(Sender: TObject);

      // Menu items
      procedure acnRunExecute(Sender: TObject);
      procedure acnSaveHistogramExecute(Sender: TObject);
      procedure acnPrintHistogramExecute(Sender: TObject);
      procedure acnCopyHistogramExecute(Sender: TObject);
      procedure acnQuitExecute(Sender: TObject);
      procedure acnAboutExecute(Sender: TObject);

      // Data entry
      procedure processTextEntry( Sender: TObject );

      // Other GUI event handlers
      procedure btnTestClick(Sender: TObject);
      procedure reIterationsChange(Sender: TObject);

    protected
      // Data structures for the model
      _simInput: TDemoSimInput;

      // Properties used for running the model
      _modelIsRunning: boolean;

      // Construction/initialization/destruction
      procedure translateUI();

      // Functions for setting up and displaying the PDF/REL editors
      procedure updateMasterDisplay();
      procedure prepEditors();
      procedure clearLists();
      procedure setEditorDictionary();
      procedure showCharts();

      // Functions for running the model
      procedure clearOutputChart();
      procedure runModel();

      // Other GUI display functions
      procedure setControlsEnabled();

    public
      // Construction/initialization/destruction
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
    StrUtils,
    
    // APHI General Purpose Delphi Library
    I88n,
    MyStrUtils,
    DebugWindow,
    ControlUtils,
    DialogLongMessage,
    RegExpDefs,
    MyDelphiArrayUtils,
    MyDialogs,
    
    // APHI Delphi Library for Simulation Modeling
    ChartFunction,
    AphiRng,
    
    // Application-specific units
    FunctionEnums,
    FormAboutAphiDemo
  ;


//-----------------------------------------------------------------------------
// Construction/initialization/destruction
//-----------------------------------------------------------------------------
  constructor TFormMain.create( AOwner: TComponent );
    begin
      inherited create( AOwner );
      translateUI();

      // Set up the GUI
      //---------------
      {$IFDEF DEBUG}
        btnTest.Visible := true;
      {$ELSE}
        btnTest.Visible := false;
      {$ENDIF}

      pnlFooter.bevelOuter := bvNone;
      pnlButtons.bevelOuter := bvNone;

      reIterations.InputExpression := RE_INTEGER_INPUT;

      fraOutputHistogram.setTitle( tr( 'Model output' ) );


      // Configure the PDF editors
      //---------------------------
      // Tell the PDF editor who's the boss.
      smcPdf1.setForm( self );

      // Specify the types of functions that the editor will have.
      smcPdf1.chartType := CTPdf;
      smcPdf1.allowPdfTypesAll();
      smcPdf1.xUnits := UUnitless;

      // This is purely for aesthetic reasons.
      smcPdf1.showDisabledLabel := false;

      // When the value of a PDF editor changes, the display of the main
      // form may need to be updated to reflect the change.  These function
      // pointers allow the editor to update the main form when necessary.
      smcPdf1.ptrUpdateAppDisplay := updateMasterDisplay;
      smcPdf1.ptrUpdateAppDisplayForChartChange := setControlsEnabled;

      // Do it all again for the other function editor.
      smcPdf2.setForm( self );
      smcPdf2.chartType := CTPdf;
      smcPdf2.allowPdfTypesAll();
      smcPdf2.xUnits := UUnitless;
      smcPdf2.showDisabledLabel := false;
      smcPdf2.ptrUpdateAppDisplay := updateMasterDisplay;
      smcPdf2.ptrUpdateAppDisplayForChartChange := setControlsEnabled;


      // Set up data structures and properties
      //--------------------------------------
      // simInput is where most of the work happens.  The
      // rest of the application is extraneous and is provided
      // only for the benefit of the user.
      _simInput := TDemoSimInput.create();
      _modelIsRunning := false;


      // Show the form!
      //---------------
      updateMasterDisplay();
    end
  ;

  
  procedure TFormMain.FormCreate(Sender: TObject);
    begin
      Assert(not Scaled, 'You should set Scaled property of Form to False!');

      if( Screen.PixelsPerInch <> 96 ) then
        ScaleBy( Screen.PixelsPerInch, 96 )
      ;
    end
  ;


  procedure TFormMain.translateUI();
    begin
      // Do nothing (yet).
    end
  ;


  destructor TFormMain.destroy();
    begin
      freeAndNil( _simInput );

      inherited destroy();
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Functions for setting up the PDF/REL editors
//-----------------------------------------------------------------------------
  procedure TFormMain.updateMasterDisplay();
    begin
      // If a chart changes, the existing output is no longer valid.
      // Make it go away.
      clearOutputChart();

      // The "model" (an instance of TDemoSimInput) stores a list of all of the PDFs that it uses.
      // The PDF editors have responsibility for different PDFs in the list.
      // This function sorts out which editor will have responsibility for which PDF.
      prepEditors();

      // More sophisticated applications (like NAADSM) might also have several "models", which are
      // also stored in a list.  If this is the case, then the editors need to be aware of which
      // model the functions they are responsible for.  In this simple application, we don't
      // need to worry about any of this: there is only one model, so there is no need for a list.


      // Show the charts
      //----------------
      showCharts();

      if( 0 < _simInput.simIterations ) then
        reIterations.Text := intToStr( _simInput.simIterations )
      else
        reIterations.Text := ''
      ;

      // Set up other GUI elements
      //--------------------------
      setControlsEnabled();
    end
  ;


  procedure TFormMain.clearLists();
    begin
      smcPdf1.ClearList();
      smcPdf2.ClearList();
    end
  ;


  procedure TFormMain.setEditorDictionary();
    begin
      smcPdf1.setFunctionDict( _simInput.functionDictionary );
      smcPdf2.setFunctionDict( _simInput.functionDictionary );
    end
  ;


  procedure TFormMain.prepEditors();
    var
      it: TFunctionDictionaryIterator;
    begin
      clearLists();

      it := TFunctionDictionaryIterator.create( _simInput.functionDictionary );

      repeat
        if( nil <> it.value() ) then
          begin
            if ( not it.value().removed ) then
              begin
                case ( it.value().fn.dbField ) of
                  integer( DEMOPdf1 ):
                    begin
                       //it.value().refCounter := 0;
                       smcPdf1.appendFunction( it.value().fn );
                    end
                  ;

                  integer( DEMOPdf2 ):
                    begin
                      //it.value().RefCounter := 0;
                      smcPdf2.appendFunction( it.value().fn );
                    end
                  ;
                end; // case statement
              end
            ;
          end
        ;

        it.incr();
      until ( nil = it.value() );

      freeAndNil( it );

      setEditorDictionary();
    end
  ;


  procedure TFormMain.showCharts();
    var
      m: TDemoModel;
    begin
      m := _simInput.model;
      
      smcPdf1.showChart( m, m.pdf1, DEMOPdf1 );
      smcPdf2.showChart( m, m.pdf2, DEMOPdf2 );
    end
  ;
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Other data entry
//-----------------------------------------------------------------------------
  procedure TFormMain.processTextEntry( Sender: TObject );
  	begin
      _simInput.simIterations := myStrToInt( reIterations.text, -1 );
    end
  ;
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Running the model
//-----------------------------------------------------------------------------
  procedure TFormMain.clearOutputChart();
    begin
      fraOutputHistogram.clear();

      acnCopyHistogram.Enabled := false;
      acnPrintHistogram.Enabled := false;
      acnSaveHistogram.Enabled := false;
    end
  ;


  procedure TFormMain.runModel();
    var
      dlg: TDialogLongMessage;
      errMsg: string;
    begin
      errMsg := '';
      if( not( _simInput.validate( @errMsg ) ) ) then
        begin
          dlg := TDialogLongMessage.create(
            self,
            tr( 'Problems found' ),
            tr( 'Problems were found with this scenario.  These problems must be corrected before this scenario can be run:' ),
            errMsg
          );
          dlg.showModal();
          dlg.free();
        end
      else
        begin
          screen.Cursor := crHourglass;

          clearOutputChart();

          _modelIsRunning := true;

          setControlsEnabled();

          _simInput.run();

          // Display output
          fraOutputHistogram.populate( _simInput.output );

          acnCopyHistogram.Enabled := true;
          acnPrintHistogram.Enabled := true;
          acnSaveHistogram.Enabled := true;

          _modelIsRunning := false;
          setControlsEnabled();
          screen.Cursor := crDefault;
        end
      ;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
      // Other GUI display functions
//-----------------------------------------------------------------------------
  procedure TFormMain.setControlsEnabled();
    var
      val: boolean;
    begin
      val :=
        ( nil <> smcPdf1.chart )
      and
        ( nil <> smcPdf2.chart )
      and
        not( _modelIsRunning )
      ;

      btnTest.Enabled := val;
      acnRun.Enabled := val;
      smcPdf1.enabled := val;
      smcPdf2.enabled := val;

      acnQuit.Enabled := not( _modelIsRunning );

      application.ProcessMessages();
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Menu items and GUI events
//-----------------------------------------------------------------------------
  procedure TFormMain.acnQuitExecute(Sender: TObject);
    begin
      // Exit the application
      self.Close();
    end
  ;


  procedure TFormMain.acnRunExecute(Sender: TObject);
    begin
      runModel();
    end
  ;


  procedure TFormMain.acnSaveHistogramExecute(Sender: TObject);
    var
      success: boolean;
    begin
      if( dlgSaveWMF.Execute ) then
        begin
          success := fraOutputHistogram.saveChartToFile( dlgSaveWMF.FileName );

          if( success ) then
            begin
              if( sender = acnSaveHistogram ) then
                msgOK(
                  ansiReplaceStr( tr( 'This chart has been successfully saved to xyz.' ), 'xyz', abbrevPath( dlgSaveWMF.FileName ) ),
                  '',
                  IMGSuccess
                )
              ;
            end
          else
            msgOK(
              ansiReplaceStr( tr( 'This chart could not be saved to xyz.' ), 'xyz', abbrevPath( dlgSaveWMF.FileName ) ),
              '',
              IMGWarning
            )
          ;
        end
      ;

      fraOutputHistogram.Repaint();
      repaint();
    end
  ;


  procedure TFormMain.acnPrintHistogramExecute(Sender: TObject);
    var
      success: boolean;
    begin
      if( sender = acnPrintHistogram ) then
        begin
          if( not( dlgPrint.Execute() ) ) then exit;
        end
      ;

      success := fraOutputHistogram.printChart();

      if( success ) then
        begin
          if( sender = acnPrintHistogram ) then
            msgOK(
              tr( 'The chart has been sent to the selected printer.' ),
              '',
              IMGSuccess
            )
          ;
        end
      else
        msgOK(
          tr( 'The chart could not be printed.' ),
          '',
          IMGWarning
        )
      ;

      fraOutputHistogram.Repaint();
      repaint();
    end
  ;


  procedure TFormMain.acnCopyHistogramExecute(Sender: TObject);
    var
      success: boolean;
    begin
      success := fraOutputHistogram.copyChartToClipboard();

      if( success ) then
        begin
          if( sender = acnCopyHistogram ) then
            msgOK(
              tr( 'This chart has been successfully copied to the clipboard.' ),
              '',
              IMGSuccess
            )
          ;
        end
      else
        msgOK(
          tr( 'This chart could not be copied to the clipboard.' ),
          '',
          IMGWarning
        )
      ;

      fraOutputHistogram.Repaint();
      repaint();
    end
  ;


  procedure TFormMain.acnAboutExecute(Sender: TObject);
    var
      frm: TFormAboutAphiDemo;
    begin
      frm := TFormAboutAphiDemo.create( self );
      frm.showModal();
      freeAndNil( frm );
    end
  ;


  procedure TFormMain.reIterationsChange(Sender: TObject);
    begin
      clearOutputChart();
    end
  ;


  procedure TFormMain.btnTestClick(Sender: TObject);
    var
      i: integer;
    begin
      for i := 0 to 4999 do
        dbcout2( _simInput.model.pdf1.rand() )
      ;
    end
  ;
//-----------------------------------------------------------------------------



end.

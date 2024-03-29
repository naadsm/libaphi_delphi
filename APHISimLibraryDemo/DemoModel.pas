unit DemoModel;

(*
DemoModel.pas
-------------
Begin: 2008/12/11
Last revision: $Date: 2009-11-10 01:23:23 $ $Author: areeves $
Version: $Revision: 1.5 $
Project: APHI Delphi Library for Simulation Modeling: Demo application
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2008 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
*)

interface

  uses
    // APHI Delphi Library for Simulation Modeling
    Models, // Defines the based classes
    ChartFunction,
    ProbDensityFunctions,
    FunctionDictionary,
    
    // Application-specific units
    FunctionEnums
  ;

  type TDemoModel = class( TModelWithFunctions )
    protected
      // Functions/properties for model setup
      //-------------------------------------
      _pdf1Name: string;
      _pdf2Name: string;
      
      // Variables used while running the model
      //---------------------------------------
      // (None are needed in this simple model)

      // Outputs generated by the model
      //-------------------------------
      _result: integer;

      // Subroutines for running the model
      //----------------------------------
      // (None are needed in this simple model)
     
      // Protected functions overridden from TModel
      //-------------------------------------------
      function getUpdated(): boolean; override;
      procedure setUpdated( val: boolean ); override;      
      
      // Property getters/setters for probability density functions
      //-----------------------------------------------------------
      procedure setPdf1Name( val: string );
      procedure setPdf2Name( val: string );

      function getPdf1Name(): string;
      function getPdf2Name(): string;
      
      function getPdf1(): TPdf;
      function getPdf2(): TPdf;
      
    public
      // Construction/initialization/destruction
      //----------------------------------------
      constructor create( sim: TObject );
      destructor destroy(); override; 

      // Public functions overridden from TModel and TModelWithFunctions
      //----------------------------------------------------------------
      function validate( err: PString = nil ): boolean; override;
      procedure debug(); override;
      //function populateDatabase( db: TSMDatabase; const forceInsert: boolean = false ): integer; reintroduce;
      procedure setChart( const whichChart: TSMChart; fn: TChartFunction; addlInfo: integer = -1 ); override;
      function chart( const whichChart: TSMChart; addlInfo: integer = -1 ): TChartFunction; override;
      procedure removeChart( const chartName: string ); override;
      function hasChartName( const chartName: string; const whichChart: TSMChart ): boolean; override;
      function functionsAreValid(): boolean; override;
      
      // Functions for running the model
      //--------------------------------
      procedure run();
      
      // Properties
      //-----------      
      property pdf1Name: string read getPdf1Name write setPdf1Name;
      property pdf2Name: string read getPdf2Name write setPdf2Name;      

      property pdf1: TPdf read getPdf1;
      property pdf2: TPdf read getPdf2;

      property result: integer read _result; 
    end
  ;
  

implementation

  uses
    // Standard Delphi units
    SysUtils,
    Math,
    
    // APHI General Purpose Delphi Library
    I88n,
    MyStrUtils,
    DebugWindow,

    // APHI Delphi Library for Simulation Modeling
    AphiRng,

    // Application-specific units
    DemoSimInput
  ;


  const
    DBSHOWMSG: boolean = false; // Set to true to enable debugging messages for this unit.


//-----------------------------------------------------------------------------
// Construction/destruction
//-----------------------------------------------------------------------------
  constructor TDemoModel.create( sim: TObject );
    begin
      inherited create();

      _sim := sim;
      
      pdf1Name := '';
      pdf2Name := '';
    end
  ;
  
  
  destructor TDemoModel.destroy();
    begin
      // The function dictionary is freed elsewhere.
      // Functions are handled by the function dictionary:
      // don't free them here, but do decrement their counters.
      decrFnRefCounter( pdf1Name );
      decrFnRefCounter( pdf2Name );

      inherited destroy();    
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Protected functions overridden from TModel
//-----------------------------------------------------------------------------
  // FIX ME: What's up with this??
  function TDemoModel.getUpdated(): boolean;
    begin
      result := _updated;
    end
  ;


  procedure TDemoModel.setUpdated( val: boolean );
    begin
      (_sim as TDemoSimInput).updated := val;
      inherited setUpdated( val );
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Chart handling: public functions overridden from TModelWithFunctions
//-----------------------------------------------------------------------------
  procedure TDemoModel.setChart( const whichChart: TSMChart; fn: TChartFunction; addlInfo: integer = -1 );
    var
      newName: string;
    begin
      if( nil = fn ) then
        newName := ''
      else
        newName := fn.name
      ;

      case whichChart of
        DEMOPdf1: self.pdf1Name := newName;
        DEMOPdf2: self.pdf2Name := newName;
        else
          raise exception.create( 'Unrecognized whichChart (' + intToStr( integer( whichChart ) ) + ') in TProductionType.setChart' )
        ;
      end;

      setUpdated( true );
    end
  ;
  
  
  function TDemoModel.chart( const whichChart: TSMChart; addlInfo: integer = -1 ): TChartFunction;
    var
      ret_val: TChartFunction;
    begin
      ret_val := nil;

      if ( self.fnDictionary <> nil ) then
        begin
          case whichChart of
            DEMOPdf1:
              if ( self.fnDictionary.contains( self.pdf1Name ) ) then
                ret_val := self.fnDictionary.value( self.pdf1Name ).fn
              ;
            DEMOPdf2:
              if ( self.fnDictionary.contains( self.pdf2Name ) ) then
                ret_val := self.fnDictionary.value( self.pdf2Name ).fn
              ;
          end;
        end
      ;
        
      result := ret_val;
    end
  ;
 
 
  procedure TDemoModel.removeChart( const chartName: string );
    begin
      if( chartName = self.pdf1Name ) then self.pdf1Name := '';
      if( chartName = self.pdf2Name ) then self.pdf2Name := '';
      
      // The _updated flag will be set by the properties above, if necessary  
    end
  ;


  function TDemoModel.hasChartName( const chartName: string; const whichChart: TSMChart ): boolean;
    begin
      result := false;
      
      case whichChart of
        DEMOPdf1: result := ( chartName = self.pdf1Name );
        DEMOPdf2: result := ( chartName = self.pdf2Name );
      end;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Validation and debugging: public functions overridden from TModel
//-----------------------------------------------------------------------------
  function TDemoModel.validate( err: PString = nil ): boolean;
    var
      submsg: string;
    begin
      result := true;

      submsg := ''; 
      if( nil = pdf1 ) then
        begin
          if( nil <> err ) then err^ := err^ + '  ' + tr( 'PDF 1 is not set.' ) + endl;
          result := false;
        end
      else if( not( pdf1.validate( @submsg ) ) ) then
        begin
          if( nil <> err ) then err^ := err^ + '  ' + tr( 'PDF 1 is not valid:' ) + submsg + endl;
          result := false;
        end
      ;
      
      submsg := '';
      if( nil = pdf2 ) then
        begin
          if( nil <> err ) then err^ := err^ + '  ' + tr( 'PDF 2 is not set.' ) + endl;
          result := false;
        end
      else if( not( pdf2.validate( @submsg ) ) ) then
        begin
          if( nil <> err ) then err^ := err^ + '  ' + tr( 'PDF 2 is not valid:' ) + submsg + endl;
          result := false;
        end
      ;

      if( false = result ) then
        begin
          if( nil <> err ) then
            err^ := tr( 'Model is not valid:' ) + endl + err^ + endl;
        end
      ;
    end
  ;

  
  function TDemoModel.functionsAreValid(): boolean;
    begin
      result := true;
      
      if( fnDictionary.contains( _Pdf1Name ) ) then
        begin
          if( not( fnDictionary.value( _Pdf1Name ).fn is TPdf ) ) then
            begin
              setPdf1Name( '' );
              result := false;
            end
          ;
        end
      ;

      if( fnDictionary.contains( _Pdf2Name ) ) then
        begin
          if( not( fnDictionary.value( _Pdf2Name ).fn is TPdf ) ) then
            begin
              setPdf2Name( '' );
              result := false;
            end
          ;
        end
      ;
    end
  ;


  procedure TDemoModel.debug();
    begin
      dbcout( '-------------BEGIN DEMO MODEL', true );

      if( nil <> pdf1 ) then
        begin
          dbcout( 'PDF 1', true );
          pdf1.debug();
        end
      else
        dbcout( 'NO PDF 1 DEFINED', true )
      ;

      if( nil <> pdf2 ) then
        begin
          dbcout( endl + 'PDF 2', true );
          pdf2.debug();
        end
      else
        dbcout( endl + 'NO PDF 2 DEFINED', true )
      ;

      dbcout( '-------------END DEMO MODEL', true );    
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// PDF-related properties
//-----------------------------------------------------------------------------
  procedure TDemoModel.setPdf1Name( val: string );
    begin
      val := trim( val );

      // Decrement the reference counter for the old function...
      if( self.fnDictionary.contains( _pdf1Name ) ) then
        dbcout2( 'XX Ref counter = ' + intToStr( self.fnDictionary.value( _pdf1Name ).refCounter ) )
      else
        dbcout2( 'There is no function called "' + val + '"' )
      ;

      decrFnRefCounter( _pdf1Name );
      // ...and increment the reference counter for the new function.
      incrFnRefCounter( val );

      _pdf1Name := val;

      setUpdated( true );
    end
  ;


  procedure TDemoModel.setPdf2Name( val: string );
    begin
      val := trim( val );

      // Decrement the reference counter for the old function...
      decrFnRefCounter( _pdf2Name );
      // ...and increment the reference counter for the new function.
      incrFnRefCounter( val );

      _pdf2Name := val;

      setUpdated( true );
    end
  ;


  function TDemoModel.getPdf1Name(): string; begin result := _pdf1Name; end;
  function TDemoModel.getPdf2Name(): string; begin result := _pdf2Name; end;


  function TDemoModel.getPdf1(): TPdf;
    begin
      if( nil = fnDictionary ) then
        result := nil
      else
        begin
          if( fnDictionary.contains( _pdf1Name ) ) then
            begin
              if( fnDictionary.value( _pdf1Name ).fn is TPdf ) then
                result := fnDictionary.value( _pdf1Name ).fn as TPdf
              else
                begin
                  setPdf1Name( '' );
                  result := nil;
                end
              ;
            end
          else
            result := nil
          ;
        end
      ;
    end
  ;


  function TDemoModel.getPdf2(): TPdf;
    begin
      if( nil = fnDictionary ) then
        result := nil
      else
        begin
          if( fnDictionary.contains( _pdf2Name ) ) then
            begin
              if( fnDictionary.value( _pdf2Name ).fn is TPdf ) then
                result := fnDictionary.value( _pdf2Name ).fn as TPdf
              else
                begin
                  setPdf2Name( '' );
                  result := nil;
                end
              ;
            end
          else
            result := nil
          ;
        end
      ;
    end
  ;
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Running the model
//-----------------------------------------------------------------------------
  procedure TDemoModel.run();
    begin
      _result := pdf1.randNonNegInt() + pdf2.randNonNegInt();
    end
  ;
//-----------------------------------------------------------------------------


end.


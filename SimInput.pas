{*
SimInput.pas
------------
Begin: 2008/09/09
Last revision: $Date: 2010-10-27 19:39:25 $ $Author: rhupalo $
Version number: $Revision: 1.9 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2010 Animal Population Health Institute, Colorado State University
                                       
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
--------------------------------------------------
This unit defines the base class for the primary data structure
business object for accessing and managing all scenario input parameters
accessed while running a simulation.
}

 (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator
 *)


unit SimInput;

interface

  {$INCLUDE Defs.inc}

  uses
    // Standard Delphi units
    SysUtils,

    // APHI General Purpose Delphi Library
    FunctionPointers,

    // APHI Delphi Library for Simulation Modeling
    {$IFDEF DATABASE_ENABLED}
    ModelDatabase,
    {$ENDIF}
    FunctionDictionary
  ;
  
  type
  {*
    The base class for SMSimulationInput, the primary data structure
    business object for accessing and managing all scenario input parameters
  }
  TSimInput = class
    protected
      _updated: boolean;            /// indicates whether any chart functions were changed, so the database can be updated

      // Simulation inputs (properties)
      //--------------------------------
      _fnDictionary: TFunctionDictionary; /// Container for chart functions

      _useFixedRandomSeed: boolean; /// whether the random number generator is to be initialized with a user-defined seed
      _randomSeed: integer;         /// user-defined random seed number

      _simIterations: integer;      /// the number of simulation iterations that are to be completed
      
      // properties used while the simulation is running
      //-------------------------------------------------
      _currentIteration: integer;       /// current value of the iteration counter, a value between 0 and _simIterations
      _lastCompleteIteration: integer;  /// the last value of the iteration counter, indicating the number of iterations completed
      
      // Functions for internal use
      //---------------------------
      procedure initialize();

      // Properties
      //-----------
      function getUseFixedRandomSeed(): boolean;
      function getRandomSeed(): integer;

      procedure setUseFixedRandomSeed( const val: boolean );
      procedure setRandomSeed( const val: integer );

      function getSimIterations(): integer;
      procedure setSimIterations( const val: integer );
      
      function getUpdated(): boolean; virtual;
      procedure setUpdated( val: boolean ); virtual;
      
      // properties used while the simulation is running
      //-------------------------------------------------
      function getLastCompleteIteration(): integer; 
      
    public     
      constructor create(); overload;
      constructor create( const src: TSimInput ); overload;
      
      destructor destroy(); override;

      /// To be implemented by decendant class to validate chart and perhaps other input parameters
      function validate( msg: PString = nil ): boolean; virtual; abstract;

      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( const updateAction: TDBUpdateActionType ): boolean; virtual; abstract;
      // FIX ME: consider making removeDbFunction virtual but not abstract.
      function removeDbFunction( const fnID: integer ): boolean; virtual; abstract;
      {$ENDIF}

      // Input properties
      //------------------

      /// provides read-access to _fnDictionary
      property functionDictionary: TFunctionDictionary read _fnDictionary;
      /// provides access _useFixedRandomSeed
      property useFixedRandomSeed: boolean read getUseFixedRandomSeed write setUseFixedRandomSeed;
      /// provides access to _randomSeed
      property randomSeed: integer read getRandomSeed write setRandomSeed;
      /// provides access to _simIterations
      property simIterations: integer read getSimIterations write setSimIterations;
      /// provides access to _updated
      property updated: boolean read getUpdated write setUpdated;

      // Functions used when the model is running
      //-----------------------------------------
      /// To be used by desendant classes but not necessarily implemented yet
      procedure run( fn: TObjFnVoid0 = nil ); virtual; abstract;

      // Properties used with outputs
      //------------------------------
      /// read-only access to _lastCompleteIteration
      property lastCompleteIteration: integer read getLastCompleteIteration;
      /// read-only access to _currentIteration
      property currentIteration: integer read _currentIteration;
    end
  ;

implementation

  uses
    DebugWindow
  ;


  /// Used by constuctors to initialize private members to safe but not necessarily usefull values
  procedure TSimInput.initialize();
    begin
      _simIterations := -1;
      _updated := false;
      
      _useFixedRandomSeed := false;
      _randomSeed := 527;
    end
  ;

  /// Creates an initialized object
  constructor TSimInput.create();
    begin
      inherited create();
      initialize();

      _fnDictionary := TFunctionDictionary.create( self );
    end
  ;

  /// destroys the object and frees memory
  destructor TSimInput.destroy();
    begin
      freeAndNil( _fnDictionary );
      inherited destroy();
    end
  ;
  
  {*
    Creates a siminput object and copies values fom src
    @param src source simInput object to copy data from
  }
  constructor TSimInput.create( const src: TSimInput );
    begin
      inherited create();
      
      _simIterations := src._simIterations;
      
      _useFixedRandomSeed := src._useFixedRandomSeed;
      _randomSeed := src._randomSeed;
      
      _fnDictionary := TFunctionDictionary.create( src._fnDictionary, self, false );
      
      _updated := src._updated;
    end
  ;

  /// Get function for property useFixedRandomSeed, returning the value of _useFixedRandomSeed
  function TSimInput.getUseFixedRandomSeed(): boolean; begin result := _useFixedRandomSeed; end;
  /// Get function for property randomSeed, returning the value of _randomSeed
  function TSimInput.getRandomSeed(): integer; begin result := _randomSeed; end;
  /// Get function for property simIterations, returning the value of _simIterations
  function TSimInput.getSimIterations(): integer; begin result := _simIterations; end;

  /// Set function for property useFixedRandomSeed, setting _useFixedRandomSeed to val
  procedure TSimInput.setUseFixedRandomSeed( const val: boolean );
    begin
      if( val <> _useFixedRandomSeed ) then
        begin
          _useFixedRandomSeed := val;
          setUpdated( true );
        end
      ;
    end
  ;

  /// Set function for property randomSeed, setting _randomSeed to val
  procedure TSimInput.setRandomSeed( const val: integer );
    begin
      if( val <> _randomSeed ) then
        begin
          _randomSeed := val;
          setUpdated( true );
        end
      ;
    end
  ;

  /// Set function for property simIterations, setting _simIterations to val
  procedure TSimInput.setSimIterations( const val: integer );
    begin
      if( val <> _simIterations ) then
        begin
          _simIterations := val;
          setUpdated( true );
        end
      ;
    end
  ;

  {*
    Get function for property updated returning the value of _updated
    @comment If _updated is false the value of _fnDictionary.updated is
    returned as a safety measure.
  }
  function TSimInput.getUpdated(): boolean;
    begin
      result := _updated;

      if( false = result ) then
        result := _fnDictionary.updated
      ;
    end
  ;

  /// Set function for property updated, setting _updated to val
  procedure TSimInput.setUpdated( val: boolean );
    begin
      _updated := val;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// properties used while the simulation is running
//-----------------------------------------------------------------------------
  /// Get function for property lastCompleteIteration, returning the value of _lastCompleteIteration
	function TSimInput.getLastCompleteIteration(): integer; begin result := _lastCompleteIteration; end;
//-----------------------------------------------------------------------------


end.

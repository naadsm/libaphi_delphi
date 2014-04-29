{*
Models.pas
----------
Begin: 2005/06/06
Last revision: $Date: 2010-10-21 19:08:37 $ $Author: rhupalo $
Version number: $Revision: 1.14 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2005 - 2010 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.


This unit defines bases classes for project specific sm_model input parameter classes
}

  (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator
  *)


unit Models;

{$INCLUDE Defs.inc}

interface

  uses
    Contnrs,

    FunctionDictionary,

    {$IFDEF DATABASE_ENABLED}
    ModelDatabase,
    {$ENDIF}
    ProbDensityFunctions,
    ChartFunction,
    FunctionEnums
  ;


  type
  /// Base class for TModelWithFunctions and sm_model classes that do not include chart functions
  TModel = class
    protected
      _updated: boolean;  /// to indicate whether data has changed so the database can be updated
      _sim: TObject;      /// access to the model input parameters, see TSMSimulationInput

      function getSim(): TObject;
      procedure setSim( sim: TObject ); virtual;

      procedure setUpdated( val: boolean ); virtual;
      /// To be used by descendant classes to scan their member collections and indcate if any updates were made
      function getUpdated(): boolean; virtual; abstract;

      function getIsValid(): boolean;

    public
      constructor create(); overload;
      constructor create( const src: TModel ); overload;
      /// To be used by descendant classes to send messages to the debugging output window
      procedure debug(); virtual; abstract;

      /// To be used by descendant to validate user input of parameter values
      function validate( err: pstring = nil ): boolean; virtual; abstract;

      function ssXml(): string; virtual;

      {$IFDEF DATABASE_ENABLED}
      /// To be used by descendant classes usually to return the database ID number of the model.  This ID may be newly assigned. 
      function populateDatabase( const updateAction: TDBUpdateActionType ): integer; virtual; abstract;
      {$ENDIF}

      property sim: TObject read getSim write setSim; /// access to _sim
      property updated: boolean read getUpdated write setUpdated;  /// access to _updated
      property isValid: boolean read getIsValid;     /// read-only flag indicating if the input parameter values are valid
    end
  ;


  type
  /// Base class for sm_model classes that include chart functions in their parameter collections
  TModelWithFunctions = class( TModel )
    protected
      function getFnDictionary(): TFunctionDictionary;
      /// To be used by descendant classes to return a set of charts from their input parameter collections
      function getChartSet(): TChartSet; virtual; abstract;

    public
      constructor create(); overload;
      constructor create( const src: TModelWithFunctions ); overload;

      // Needed for chart editing
      //-------------------------
      /// To be used by descendant classes to remove a chart from their charts collection
      procedure removeChart( const chartName: string ); virtual; abstract;
      /// To be used by descendant classes to return a chart from their charts collection
      function chart( const whichChart: TSMChart; addlInfo: integer = -1 ): TChartFunction; virtual; abstract;
      /// To be used by descendant classes to associate a chart with a particular input parameter
      procedure setChart( const whichChart: TSMChart; fn: TChartFunction; addlInfo: integer = -1 ); virtual; abstract;

      // This function should not need to be overridden often.
      // See NAADSM ProductionType.pas for an example where it makes sense to do so.
      // This function also must be reimplemented if the same chart can be used in multiple fields within a model.
      procedure changeChart(
        const whichChart: TSMChart;
        const oldChartName: string;
        newChart: TChartFunction;
        addlInfo: integer = -1
      ); virtual;

      ///  Used by descendant classes to ensure that information read from the database produces the right kinds of functions.
      function functionsAreValid(): boolean; virtual; abstract;

      /// Used by descendant classes to determine whether a model uses the specified chart
      function hasChartName( const chartName: string; const whichChart: TSMChart ): boolean; virtual; abstract;

      // Used to keep track of how many times individual functions in the dictionary are used.
      procedure decrFnRefCounter( val: string );
      procedure incrFnRefCounter( val: string );
      function fnRefCounter( val: string ): integer;

      property fnDictionary: TFunctionDictionary read getFnDictionary; /// read-only access to the function dictionary of _sim
      property chartSet: TChartSet read getChartSet; /// returns the chart collection of the model
    end
  ;


  // FIX ME: Consider making TModelList a derivative of QIntegerObjectMap (using the
  // model IDs as the keys) to improve performance for random access to list elements.
  // Adding new elements may be a minor problem, since the permanent IDs are assigned
  // only upon insertion into the database.
  type
  /// Base class for container to manage a collection of models
  TModelList = class( TObjectList )
    protected
      _sim: TObject;  /// access to the model input parameters, see TSMSimulationInput

      function getFnDictionary(): TFunctionDictionary;

      // Typical list functions
      procedure setObject( index: integer; item: TModel );
      function getObject( index: integer ): TModel;
      
    public
      constructor create(); overload;
      constructor create( const src: TModelList ); overload;

      // This list owns its objects (see the contructors), so
      // objects are automatically freed when the list is destroyed.

      // Typical list functions
      function append( dm: TModel ): integer; virtual;
      property objects[ index: integer]: TModel read getObject write setObject; default;
      function at( i: integer ): TModel;

      // Needed for chart editing
      procedure removeChart( const chartName: string ); virtual;
      procedure changeChart(
        const whichChart: TSMChart;
        const oldChartName: string;
        newChart: TChartFunction;
        addlInfo: integer = -1
      ); virtual;

      // Validation and debugging
      function validate( err: PString = nil ): boolean; virtual;
      /// To be used by descendant classes to indicate whether all the models have valid chart functions
      function functionsAreValid(): boolean; virtual; abstract;
      procedure debug(); virtual;

      // Properties
      /// read-only access to the function dictionary of _sim
      property fnDictionary: TFunctionDictionary read getFnDictionary;
      /// read-only access to _sim
      property sim: TObject read _sim;
    end
  ;

  
  // FIX ME: iterator class has no error checking.
  type
  /// To iterate through a list of models
  TModelListIterator = class
    protected
      _list: TModelList;        /// model list to navigate
      _currentIndex: integer;   /// current position in _list

      function getCount(): word;
      function getIsEmpty(): boolean;

      function _toFirst(): TModel;
      function _toLast(): TModel;
      function _current(): TModel;

    public
      constructor create( list: TModelList );
      destructor destroy(); override;

      procedure incr();
      procedure decr();

      property count: word read getCount;  /// number of models in the list
      property isEmpty: boolean read getIsEmpty; /// whether the listis empty
      property currentIndex: integer read _currentIndex; /// current position in the list
    end
  ;
  

implementation

  uses
    SysUtils,

    DebugWindow,
    
    SimInput
  ;

  /// Creates an empty model, _sim is set to nil.
  constructor TModel.create();
    begin
      inherited create();

      _updated := false;
      _sim := nil;
    end
  ;

  {*
    Creates an empty model, copying the value of _updated from src
    @param src source TModel object
    @comment It is up to the descendant class to do more with src
  }
  constructor TModel.create( const src: TModel );
    begin
      inherited create();
      
      _sim := nil;
      _updated := src._updated;
    end
  ;

  {*
    Get function for property sim
    @return object reference to _sim
  }
  function TModel.getSim(): TObject;
    begin
      if( nil = _sim ) then
        dbcout( '_sim is nil in TModel.getSim()', true )
      ;
      result := _sim;
    end
  ;

  {*
    Set function for property sim
    @param sim an instance of TSMSimulationInput
  }
  procedure TModel.setSim( sim: TObject );
    begin
      if( nil = sim ) then
        dbcout( 'sim is nil in TModel.setSim()', true )
      ;
      _sim := sim;
    end
  ;

  {*
    Get function for property isValid
    @return true if all model parameter values are valid, else false
  }
  function TModel.getIsValid(): boolean; begin result := validate( nil ); end;

  {*
    Set function for property isUpdated
    @param val value to set _updated
  }
  procedure TModel.setUpdated( val: boolean ); begin _updated := val; end;


  {*
    Cautionary measure for a virtual method that should be over-ridden by descendant classes
    @return nothing
    @throws exception indicating error
  }
  function TModel.ssXml(): string;
    begin
      raise exception.create( '"Abstract" error: call made to TModel.ssXml()' );
    end
  ;


  /// Creates an empty model, _sim is set to nil.
  constructor TModelWithFunctions.create();
    begin
      inherited create();
    end
  ;

  {*
    Creates an empty model, copying the value of _updated from src
    @param src source TModelWithFunctions object
    @comment It is up to the descendant class to do more with src
  }
  constructor TModelWithFunctions.create( const src: TModelWithFunctions );
    begin
      inherited create( src );
    end
  ;
  

  {*
     Get function for property fnDictionary
     @return object reference to the function dictionary of _sim
  }
  function TModelWithFunctions.getFnDictionary(): TFunctionDictionary;
    begin
      if( nil <> _sim ) then
        result := (_sim as TSimInput).functionDictionary
      else
        begin
          dbcout( '_sim is nil in TModelWithFunctions.getFnDictionary()', true );
          result := nil;
        end
      ;
    end
  ;


  {*
    Changes the pdf or relational function assigned to whichChart from oldChartName to newChart
    @param whichChart an item from a set of TSMCharts defined in FunctionEnums
    @param oldChartName current chart name to be changed
    @param newChart a PDF and relational function object
    @param addlInfo not used in this instance, assign -1
    @comment This function should not need to be overridden often.
    See NAADSM ProductionType.pas for an example where it makes sense to do so.
    This function also must be reimplemented if the same chart can be used in multiple fields within a model.
  }
  procedure TModelWithFunctions.changeChart(
        const whichChart: TSMChart;
        const oldChartName: string;
        newChart: TChartFunction;
        addlInfo: integer = -1
      );                                                                           
    begin
      if( self.hasChartName( oldChartName, whichChart ) ) then
        self.setChart( whichChart, newChart, addlInfo )
      ;
    end
  ;


  {*
    Decrements the reference counter for a function
    @param val key for the function in the function dictionary
  }
  procedure TModelWithFunctions.decrFnRefCounter( val: string );
    begin
      val := trim( val );

      if( '' <> val ) then
        begin
          if( fnDictionary.contains( val ) ) then
            fnDictionary.value( val ).decrRefCounter()
          ;
        end
      ;
    end
  ;


  {*
    Increments the reference counter for a function
    @param val key for the function in the function dictionary
  }
  procedure TModelWithFunctions.incrFnRefCounter( val: string );
    begin
      val := trim( val );

      if( '' <> val ) then
        begin
          if( fnDictionary.contains( val ) ) then
            fnDictionary.value( val ).incrRefCounter()
          else
            raise exception.create( 'Missing function dictionary item (' + val + ') in TModelWithFunctions.incrFnRefCounter' )
          ;
        end
      ;
    end
  ;


  {*
    Provides the reference counter value of function having the key of val
    @param val key for the function in the function dictionary
    @return the reference count of the function
  }
  function TModelWithFunctions.fnRefCounter( val: string ): integer;
    begin
      val := trim( val );

      if( '' <> val ) then
        begin
          if( fnDictionary.contains( val ) ) then
            result := fnDictionary.value( val ).refCounter
          else
            raise exception.create( 'Missing function dictionary item (' + val + ') in TModelWithFunctions.fnRefCounter' )
          ;
        end
      else
        result := 0
      ;
    end
  ;

  /// Creates an empty model list
  constructor TModelList.create();
    begin
      inherited create( true );
    end
  ;

  {*
    Creates an empty model, copying the value of _updated from src
    @param src source TModelList object
    @comment It is up to the descendant class to do more with src
  }
  constructor TModelList.create( const src: TModelList );
    begin
      inherited create( true );
    end
  ;


  {*
    Cautionary measure for a virtual method that should be over-ridden by descendant classes
    @throws exception indicating error
  }
  procedure TModelList.removeChart(
        const chartName: string
      );
    begin
      raise exception.create( '"Abstract" error: call made to TModelList.removeChart()' );
    end
  ;


  {*
    Cautionary measure for a virtual method that should be over-ridden by descendant classes
    @throws exception indicating error
  }
  procedure TModelList.changeChart(
        const whichChart: TSMChart;
        const oldChartName: string;
        newChart: TChartFunction;
        addlInfo: integer = -1
      );
    begin
      raise exception.create( '"Abstract" error: call made to TModelList.changeChart()' );
    end
  ;


  {*
    Cautionary measure for a virtual method that should be over-ridden by descendant classes
    @throws exception indicating error
  }
  function TModelList.validate( err: PString = nil ): boolean;
    begin
      raise exception.create( '"Abstract" error: call made to TModelList.validate()' );
      result := false;
    end
  ;


  {*
     Get method for property fnDictionary
     @return object reference to the function dictionary of _sim
  }
  function TModelList.getFnDictionary(): TFunctionDictionary;
    begin
      if( nil <> _sim ) then
        result := (_sim as TSimInput).functionDictionary
      else
        begin
          dbcout( '_sim is nil in TModelList.getFnDictionary()', true );
          result := nil;
        end
      ;
    end
  ;


  {*
    Cautionary measure for a virtual method that should be over-ridden by descendant classes
    @throws exception indicating error
  }
  procedure TModelList.debug();
    begin
      raise exception.create( '"Abstract" error: call made to TModelList.debug()' );
    end
  ;


//-----------------------------------------------------------------------------
// List: Typical list functions
//-----------------------------------------------------------------------------
  {*
    Returns object reference to model at location i in the list
  }
  function TModelList.at( i: integer ): TModel;
    begin
      result := getObject( i );
    end
  ;


  {*
    Appends a model object to the list
    @param dm object to append
    @return index value of the object whch was added
  }
  function TModelList.append( dm: TModel ): integer;
    begin
      result := inherited Add( dm );
    end
  ;


  {*
     Adds a model object to the list at list position index
     @param index list position to put item
     @param item model to be added
  }
  procedure TModelList.setObject( index: integer; item: TModel );
    begin
      inherited SetItem( index, item );
    end
  ;


  {*
     Returns an object reference to the model at position index
     @throws An exception is raised if value of index is out of bounds
  }
  function TModelList.getObject( index: integer ): TModel;
    begin
      if( ( index < 0 ) or ( index > self.Count - 1 ) ) then
        begin
          raise exception.Create( 'Index out of bounds (' + intToStr( index ) + ') in TModelList with ' + intToStr( self.Count ) + ' items.' );
          result := nil;
        end
      else
        result := inherited GetItem( index ) as TModel
      ;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Iterator
//-----------------------------------------------------------------------------

  {*
    Creates a list iterator for list
    @param list model list that will be navigated
  }
  constructor TModelListIterator.create( list: TModelList );
    begin
      _list := list;
      _currentIndex := 0;
    end
  ;


  /// Frees object from memory
  destructor TModelListIterator.destroy();
    begin
      // Nothing to do in here
      inherited destroy();
    end
  ;


  /// Returns an object reference to the first model in the list or nil if the list is empty
  function TModelListIterator._toFirst(): TModel;
    begin
      if( 0 = _list.Count ) then
        begin
          _currentIndex := -1;
          result := nil;
        end
      else
        begin
          _currentIndex := 0;
          result := _list.at(0) as TModel;
        end
      ;
    end
  ;

  /// Returns an object reference to the last model in the list or nil if the list is empty
  function TModelListIterator._toLast(): TModel;
    begin
      if( 0 = _list.Count ) then
        begin
          _currentIndex := -1;
          result := nil;
        end
      else
        begin
          _currentIndex := _list.count-1;
          result := _list.at( _currentIndex ) as TModel;
        end
      ;
    end
  ;

  /// Returns an object reference to the current model in the list or nil if the list is empty 
  function TModelListIterator._current(): TModel;
    begin
      if( _currentIndex < 0 ) or ( _currentIndex > _list.Count -1 ) then
        result := nil
      else
        result := _list.at( _currentIndex ) as TModel
      ;
    end
  ;

  /// Increments _currentIndex
  procedure TModelListIterator.incr();
    begin
      inc( _currentIndex );
    end
  ;

  /// Decrements _currentIndex
  procedure TModelListIterator.decr();
    begin
      dec( _currentIndex );
    end
  ;

  /// Returns the number of models in the list
  function TModelListIterator.getCount(): word;
    begin
      result := _list.count;
    end
  ;

  /// True if list count is 0 else false
  function TModelListIterator.getIsEmpty(): boolean;
    begin
      result := (_list.count = 0 );
    end
  ;
//-----------------------------------------------------------------------------

end.

{*
FunctionDictionary.pas
----------------------
Begin: 2006/01/04
Last revision: $Date: 2012-10-23 22:24:38 $ $Author: areeves $
Version number: $Revision: 1.31 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2006 - 2012 Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

The units in this class deal with managing the chart functions associated with the input parameters.
}

 (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator
 *)


unit FunctionDictionary;

{$INCLUDE Defs.inc}

interface

  uses
    Classes,

    QStringMaps,

    {$IFDEF DATABASE_ENABLED}
    SqlClasses,
    {$ENDIF}

    ChartFunction
  ;

	type
  /// Reference to a PDF or  REL chart function object, and its status and useage
  TFunctionDictionaryItem = class
  	protected
      _refCounter: integer;  /// reference counter of the number of input parameters using this function
      _removed: boolean;     /// whether the function has been removed (no longer reference by any input parameters)
      _new: boolean;         /// whether the function is new, not yet in the database
      _modified: boolean;    /// whether the function has been modified, flag indicating DB needs to be updated

      function getRefCounter(): integer;
      procedure setRefCounter( val: integer );

      procedure setRemoved( val: boolean );
      procedure setNew( val: boolean );
      procedure setModified( val: boolean );
      function getRemoved(): boolean;
      function getNew(): boolean;
      function getModified(): boolean;

      function getUpdated(): boolean;
  	public
    	_fn: TChartFunction;  /// a chart function object

    	constructor create( val: TChartFunction; const refCounter: longword = 0 ); overload;
      constructor create( const src: TFunctionDictionaryItem ); overload;
      destructor destroy(); override;

      procedure incrRefCounter();
      procedure decrRefCounter();
      procedure resetRefCounter();

      procedure debug();

      property fn: TChartFunction read _fn;                          /// chart function
                                                                     /// number of input parameters referencing _fn
      property refCounter: integer read getRefCounter; // write setRefCounter;  // AR 10/19/09 It should never be necessary to set the refCounter.  This is handled whenever charts are changed.
      property removed: boolean read getRemoved write setRemoved;    /// indicates no input parameters reference this item so it can be removed from the database
      property new: boolean read getNew write setNew;                /// indicates that is is a new chart function so it is added to the database
      property modified: boolean read getModified write setModified; /// indicates _fn has been modified so the database is updated

      property updated: boolean read getUpdated;                     /// indicates this item has been removed, modified, or is new
  	end
  ;


  type
  /// Container for chart functions, indexed by chart name
  TFunctionDictionary = class( TQStringObjectMap )
    protected
      _sim: TObject; /// access to the model input parameters, see TSMSimulationInput

      function getUpdated(): boolean;

    public
      constructor create( sim: TObject ); overload;

      {$IFDEF DATABASE_ENABLED}
      constructor create( db: TSqlDatabase; sim: TObject ); overload;
      procedure fillFromDB( db: TSqlDatabase );
      {$ENDIF}
      constructor create( const src: TFunctionDictionary; sim: TObject; preserveReferenceCounts: boolean ); overload;
      destructor destroy(); override;

      procedure assign( Source: TObject ); override;

      function value( const key: string ): TFunctionDictionaryItem; reintroduce;
      function itemAtIndex( idx: integer ): TFunctionDictionaryItem; reintroduce;

      procedure addFunction( const fn: TChartFunction ); overload;
      procedure addFunction( const fn: TChartFunction; const fnName: string ); overload;

      // Returns the name of the function that is an identical match for fn, if there is one, otherwise returns the name of fn.
      function checkAndInsert( fn: TChartFunction ): string;

      procedure incrRefCounter( const fnName: string );
      procedure decrRefCounter( const fnName: string );

      // Returns false if the named function is not in the dictionary or is not updated.
      function functionExistsAndIsUpdated( const fnName: string ): boolean;

      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase ): boolean;
      {$ENDIF}

      procedure debug();
      procedure showRefCounts();

      property updated: boolean read getUpdated; /// indicates if any dictionary items have been updated
    end
  ;


  type
  /// Used to iterate through a TFunctionDictionary object
  TFunctionDictionaryIterator = class( TQStringObjectMapIterator )
    public
      function value(): TFunctionDictionaryItem;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function createFunctionFromDB( db: TSqlDatabase; chartID: integer ): TChartFunction;
  {$ENDIF}

implementation

  uses
    SysUtils,
    StrUtils,
    Variants,
    
    DebugWindow,
    MyStrUtils,

    ProbDensityFunctions,
    RelFunction,
    SimInput
    
    {$IFDEF DATABASE_ENABLED}
    , FunctionEnums
    {$ENDIF}
  ;

  const
    DBSHOWMSG: boolean = false; /// Set to true to enable debugging messages for this unit

  {$IFDEF DATABASE_ENABLED}

  {*
     Helper function called by TFunctionDictionary.fillFromDB to created a chart function from the database
     @param db object accessing the NAADSM parameters database
     @param chartID primary key value for table inChart
     @return a PDF or Relational chart for data associated with chartID
  }
  function createFunctionFromDB( db: TSqlDatabase; chartID: integer ): TChartFunction;
  	var
    	q: string;
      res: TSqlResult;
      row: TSqlRow;
    	pdft: TPdfType;
      pdf: TChartFunction;
      rel: TChartFunction;
  	begin
    	q := 'SELECT '
      	+ '`isPdf`, '
        + '`chartType`, '
        + '`mean`, `stddev`, '
        + '`min`, `mode`, `max`, '
        + '`alpha`, `alpha2`, `beta`, '
        + '`location`, `scale`, `shape`, '
        + '`n`, `p`, '
        + '`m`, `d`, '
        + '`dMin`, `dMax`, '
        + '`theta`, `a`, '
        + '`s`, '
        + '`xAxisUnits`, '
        + '`yAxisUnits`, '
        + '`chartID`, '
        + '`chartName`, '
        + '`fieldName`, '
        + ' `notes`'
        + 'FROM `inChart` WHERE `chartID` = ' + intToStr( chartID )
      ;

      res := TSqlResult.create( q, db );

      row := res.fetchArrayFirst();
      
      if( row.field( 'isPdf' ) ) then
        begin
          pdft := pdfType( row.field('chartType') );

          case pdft of
            // Continuous types
            //-----------------
            PdfBeta: pdf := TPdfBeta.create( row.field('alpha'), row.field('alpha2'), row.field('min'), row.field('max'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfBetaPERT: pdf := TPdfBetaPERT.create( row.field('min'), row.field('mode'), row.field('max'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfExponential: pdf := TPdfExponential.create( row.field('mean'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfGamma: pdf := TPdfGamma.create( row.field('alpha'), row.field('beta'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfGaussian: pdf := TPdfGaussian.create( row.field('mean'), row.field('stddev'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfHistogram: pdf := TPdfHistogram.create( db, chartID, chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfInverseGaussian: pdf := TPdfInverseGaussian.create( row.field('mean'), row.field('shape'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfLogistic: pdf := TPdfLogistic.create( row.field('location'),row.field('scale'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfLogLogistic: pdf := TPdfLogLogistic.create( row.field('location'), row.field('scale'), row.field('shape') , chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfLognormal: pdf := TPdfLognormal.create( row.field('mean'), row.field('stddev'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfPareto: pdf := TPdfPareto.create( row.field('theta'), row.field('a'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfPearson5: pdf := TPdfPearson5.create( row.field('alpha'), row.field('beta'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfPiecewise: pdf := TPdfPiecewise.create( db, chartID, chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfPoint: pdf := TPdfPoint.create( row.field('mode'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfTriangular: pdf := TPdfTriangular.create( row.field('min'), row.field('mode'), row.field('max'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfUniform: pdf := TPdfUniform.create( row.field('min'), row.field('max'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfWeibull: pdf := TPdfWeibull.create( row.field('alpha'), row.field('beta'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );

            // Discrete types
            //---------------
            PdfBinomial: pdf := TPdfBinomial.create( row.field('n'), row.field('p'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfDiscreteUniform: pdf := TPdfDiscreteUniform.create( row.field('dMin'), row.field('dMax'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfHypergeometric: pdf := TPdfHypergeometric.create( row.field('n'), row.field('d'), row.field('m'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfNegativeBinomial: pdf := TPdfNegativeBinomial.create( row.field('s'), row.field('p'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
            PdfPoisson: pdf := TPdfPoisson.create( row.field('mean'), chartUnitTypeFromDatabase( row.field('xAxisUnits') ) );
          else
            begin
              raise exception.create( 'Unrecognized PDF type ("' + row.field('chartType') + '") in functionFromDB' );
              pdf := nil;
            end
          end;

          if( nil <> pdf ) then
            begin
              pdf.name := trim( row.field( 'chartName' ) );
              pdf.id := row.field( 'chartID' );
              pdf.dbField := ord( strToSMChart( row.field( 'fieldName' ) ) );
              if( null <> row.field( 'notes' ) ) then
                pdf.notes := row.field( 'notes' )
              else
                pdf.notes := ''
              ;
            end
          ;

          result := pdf;
        end
      else // relational chart
      	begin
          rel := TRelFunction.create( db, chartID, chartUnitTypeFromDatabase( row.field( 'xAxisUnits') ), chartUnitTypeFromDatabase( row.field( 'yAxisUnits' ) ) );

          rel.name := trim( row.field( 'chartName' ) );
          rel.id := row.field( 'chartID' );
          rel.dbField := ord( strToSmChart( row.field( 'fieldName' ) ) );
          if( null <> row.field( 'notes' ) ) then
            rel.notes := row.field( 'notes' )
          else
            rel.notes := ''
          ;
          result := rel;
        end
      ;

      freeAndNil( res );
    end
  ;
  {$ENDIF}

//-----------------------------------------------------------------------------
// TFunctionDictionaryItem
//-----------------------------------------------------------------------------

  {*
    Constructor to create a function dictionary item for val and an input reference counter
    @param val a PDF or REL type of chart function
    @param refCounter reference counter value indicating how many input parameters reference val
  }
	constructor TFunctionDictionaryItem.create( val: TChartFunction; const refCounter: longword = 0 );
  	begin
      inherited create();

   		_fn := val;

      _refCounter := refCounter;

      _removed := false;
      _new := false;
      _modified := false;
    end
  ;

  {*
    Constructor to create a function dictionary item from src
    @param src the source FunctionDictionaryItem to copy
  }
  constructor TFunctionDictionaryItem.create( const src: TFunctionDictionaryItem );
    begin
      inherited create();

      _fn := src._fn.createCopy();

      _refCounter := src._refCounter;

      _removed := src._removed;
      _new := src._new;
      _modified := src._modified;
    end
  ;


  /// Frees the dictionary item from memory
  destructor TFunctionDictionaryItem.destroy();
  	begin
      _fn.Free(); // FIX ME 5/23: Is this right??
      inherited destroy();
    end
  ;


  /// outputs to the debug output window information regarding the function's name and reference count
  procedure TFunctionDictionaryItem.debug();
  	begin
      dbcout( endl + '--- TFunctionDictionaryItem.debug', true );

      if( nil = _fn ) then
        dbcout( '*** fn is nil!  This should not be possible!', true )
      else
        begin
          dbcout( 'Function named "' + _fn.name + '" is used ' + intToStr( refCounter ) + ' times.', true );
          _fn.debug();
        end
      ;

      dbcout( '--- Done TFunctionDictionaryItem.debug' + endl, true );
    end
  ;

  {*
    Indicates that this chart function is no longer referenced by any input parameter
    @param val sets _removed, and if true, then the reference count is set to 0
  }
  procedure TFunctionDictionaryItem.setRemoved( val: boolean );
  	begin
    	_removed := val;
  		if( _removed ) then _refCounter := 0;
    end
  ;

  /// get function for property refCounter returning the value of _refCounter
	function TFunctionDictionaryItem.getRefCounter(): integer; begin result := _refCounter; end;
  /// set function for property refCounteer setting the value of _refCounter from val
  procedure TFunctionDictionaryItem.setRefCounter( val: integer ); begin _refCounter := val; end;

  /// increments the reference counter by 1
  procedure TFunctionDictionaryItem.incrRefCounter();
    begin
      inc( _refCounter );
    end
  ;

  /// decrements reference counter by 1 and raises exception if the value falls below 0
  procedure TFunctionDictionaryItem.decrRefCounter();
    begin
      dec( _refCounter );

      if( 0 > _refCounter ) then
        raise exception.Create( '_refCounter has dropped below 0 in TFunctionDictionaryItem.decrRefCounter(): "' + self.fn.name + '"' )
      ;
    end
  ;

  /// resets the reference counter to zero
  procedure TFunctionDictionaryItem.resetRefCounter(); begin _refCounter := 0; end;

  /// set function for property new setting the value of _new from val
  procedure TFunctionDictionaryItem.setNew( val: boolean ); begin _new := val; end;
  /// set function for property modified setting the value of _modified from val
  procedure TFunctionDictionaryItem.setModified( val: boolean ); begin _modified := val; end;

  /// get function for property removed returning the value of _removed
  function TFunctionDictionaryItem.getRemoved(): boolean; begin Result := _removed; end;
  /// get function for property new returning the value of _new
  function TFunctionDictionaryItem.getNew(): boolean; begin Result := _new; end;
  /// get function for property modified returning the value of _modified
  function TFunctionDictionaryItem.getModified(): boolean; begin Result := _modified; end;

  {*
    Indicates whether the item is new, modified, or is to be removed
    @return true if the function dictionary is new or modified or is to be removed, else false
  }
  function TFunctionDictionaryItem.getUpdated(): boolean;
    begin
      result :=
        _new
      or
        _removed
      or
        _modified
      ;
    end
  ;
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// TFunctionDictionary
//-----------------------------------------------------------------------------

  {$IFDEF DATABASE_ENABLED}

  {*
    Creates a function dictionary and fills it from chart functions from the database
    @param db object accessing the NAADSM parameters database
    @param sim  access to the model input parameters, see TSMSimulationInput
  }
  constructor TFunctionDictionary.create( db: TSqlDatabase; sim: TObject );
  	begin
    	inherited create();
      _sim := sim;

      fillFromDB( db );
    end
  ;


  {*
    Helper function used by constructor to fill the dictionary with chart functions from the database
    @param db object accessing the NAADSM parameters database
  }
  procedure TFunctionDictionary.fillFromDB( db: TSqlDatabase );
  	var
      res: TSqlResult;
      row: TSqlRow;
      cht: TChartFunction;
      item: TFunctionDictionaryItem;
    begin
      res := TSqlResult.create( 'SELECT `chartID` FROM `inChart` ORDER BY `chartName`', db );

      row := res.fetchArrayFirst();
      while( row <> nil ) do
      	begin
        	cht := createFunctionFromDB( db, row.field(0) );
          item := TFunctionDictionaryItem.create( cht );
          self.insert( cht.name, item );
          row := res.fetchArrayNext();
        end
      ;

      freeAndNil( res );
    end
  ;
  {$ENDIF}

  {*
    Creates an empty function dictionary
    @param sim  access to the model input parameters, used to set _sim
  }
  constructor TFunctionDictionary.create( sim: TObject );
    begin
      inherited create();
      _sim := sim;
    end
  ;


  {*
    Creates a function dictionary by copying the items in src
    @param src source of dictionary item chart functions
    @param sim  access to the model input parameters, used to set _sim
    @param preserveReferenceCounts if true the copied reference count value of each item is maintained, else each is set to 0
  }
  constructor TFunctionDictionary.create( const src: TFunctionDictionary; sim: TObject; preserveReferenceCounts: boolean );
    var
      it: TFunctionDictionaryIterator;
    begin
      inherited create( src );
      _sim := sim;

      if( not preserveReferenceCounts ) then
        begin
          it := TFunctionDictionaryIterator.create( self );

          while( nil <> it.value() ) do
            begin
              it.value()._refCounter := 0;
              it.incr();
            end
          ;

          it.free();
        end
      ;
    end
  ;

  /// Frees the dictionary and dictionary items from memory
  destructor TFunctionDictionary.destroy();
    begin
      deleteValues();
      inherited destroy();
    end
  ;


  {*
    Returns the function dictionary item for key
    @param key the name of the chart function associated with this dictionary item
    @return the function dictionary item (and chart function) named key
  }
  function TFunctionDictionary.value( const key: string ): TFunctionDictionaryItem;
    begin
      result := inherited value(key) as TFunctionDictionaryItem;
    end
  ;


  {*
    Returns the function dictionary item for position idx
    @param idx the index position of the chart function in the dictionary
    @return the function dictionary item (and chart function) at idx
  }
  function TFunctionDictionary.itemAtIndex( idx: integer ): TFunctionDictionaryItem;
    begin
      result := inherited itemAtIndex( idx ) as TFunctionDictionaryItem;
    end
  ;


  {*
     Adds a dictionary item to self for fn
     @param fn chart function to be added to the dictionary
     @comment the key for the new item will be chart function's name
  }
  procedure TFunctionDictionary.addFunction( const fn: TChartFunction );
    var
      item: TFunctionDictionaryItem;
    begin
      item := TFunctionDictionaryItem.create( fn, 0 );
      self.insert( fn.name, item );
    end
  ;


  {*
     Adds a dictionary item to self for fn but naming it fnName rather than fn.Name
     @param fn chart function to be added to the dictionary
     @param fnName the name for the new item
  }
  procedure TFunctionDictionary.addFunction( const fn: TChartFunction; const fnName: string );
    var
      item: TFunctionDictionaryItem;
    begin
      fn.name := fnName;
      item := TFunctionDictionaryItem.create( fn, 0 );
      self.insert( fn.name, item );
    end
  ;


  {*
     Increments the reference counter of the dictionary item named fnName
     @param fnName name of the function dictionary item
     @comment an exception is raised if fnName is not in the dictionary
  }
  procedure TFunctionDictionary.incrRefCounter( const fnName: string );
    begin
      if( '' <> fnName ) then
        begin
          if( self.contains( fnName ) ) then
            self.value( fnName ).incrRefCounter()
          else
            raise exception.create( 'Missing function dictionary item (' + fnName + ') in TFunctionDictionary.incrFnRefCounter' )
          ;
        end
      ;
    end
  ;


  {*
     Decrements the reference counter of the dictionary item named fnName
     @param fnName name of the function dictionary item
     @comment an exception is NOT raised if fnName is not in the dictionary
  }
  procedure TFunctionDictionary.decrRefCounter( const fnName: string );
    begin
      if( '' <> fnName ) then
        begin
          if( self.contains( fnName ) ) then
            self.value( fnName ).decrRefCounter()
          ;
        end
      ;
    end
  ;


  {*
    Checks to see if fn has the same chart type, chart name, and axes units as an item in the list.
    If so, replace fn with the one in the list. If not add an item for fn to the dictionary.

    @param fn a chart function that may or may not already be in the list
    @return end result chart function name
    @comment WARNING: parameter fn may be freed by this function.  This behavior is by design:
    if fn is an exact duplicate of a function that already exists in the dictionary,
    the existing instance should be used instead, and the duplicate instance should
    be thrown out.

    If access to the TChartFunction selected by checkAndInsert() is needed, use
    code like this:

    fnName := dictionary.checkAndInsert( someChartFunction );
    someChartFunction := dictionary.value( fnName ).fn;
  }
  function TFunctionDictionary.checkAndInsert( fn: TChartFunction ): string;
    var
      i: integer;
      existing: TChartFunction;
      matchFound: boolean;
      dictItem: TFunctionDictionaryItem;

      // FIX ME: should this go into MyStrUtils?  See also ChartFunction.compare()
      function stripParens( str: string ): string;
        var
          strLen: integer;
          parenPos: integer;
        begin
          result := trim( str );

          strLen := length( str );
          if( ')' = charAt( result, strLen - 1 ) ) then
            begin
              parenPos := lastDelimiter( '(', result );
              if( 0 < parenPos ) then
                result := trim( leftStr( result, parenPos - 1 ) )
              ;
            end
          ;
        end
      ;
    begin
      // Is there already a function in the same field with a similar name and identical values?
      // If so, replace this function with the one already in the list.
      // If not, add this function to the list.
      existing := nil;
      matchFound := false;

      for i := 0 to self.count - 1 do
        begin
          existing := self.itemAtIndex(i).fn;

          //if( existing.dbField = fn.dbField ) then
            //begin
              if( lower( stripParens( fn.name ) ) = lower( stripParens( existing.name ) ) ) then
                begin
                  if( fn.compare( existing ) ) then
                    begin
                      matchFound := true;
                      break;
                    end
                  ;
                end
              ;
            //end
          //;
        end
      ;

      if( not( matchFound ) ) then
        begin
          dictItem := TFunctionDictionaryItem.create( fn, 0 );
          dictItem.new := true;
          self.insert( fn.name, dictItem );
          result := fn.name;
        end
      else // A match was found, use the existing function
        begin
          if( nil = existing ) then
            raise exception.Create( 'Existing function is nil in TFunctionDictionary.checkAndInsert' )
          else
            begin
              fn.Free();
              result := existing.name;
            end
          ;
        end
      ;
    end
  ;


  {*
     Populates self from Source, assumed to be of type TQStringObjectMap
     @param Source object that is cast to TQStringObjectMap, contributing function dictionary items
     @comment Fix me: Is this used to get around Delphi unit circular reference error?
  }
  procedure TFunctionDictionary.assign( Source: TObject );
    var
      src: TQStringObjectMap;
      it: TFunctionDictionaryIterator;
      newObj: TFunctionDictionaryItem;
    begin
      dbcout( '+++ Begin TFunctionDictionary.assign', DBSHOWMSG );
      src := source as TQStringObjectMap;
      it := TFunctionDictionaryIterator.create( src );

      repeat
        if( nil <> it.value() ) then
          begin
            newObj := TFunctionDictionaryItem.Create( it.value() );
            self.insert( it.key(), newObj );
          end
        ;

        it.incr();
      until ( nil = it.value() );

      it.Free();

      dbcout( '--- Done TFunctionDictionary.assign', DBSHOWMSG );
    end
  ;


  /// outputs to the debug output window information regarding each function's name and reference count
  procedure TFunctionDictionary.debug();
    var
      it: TFunctionDictionaryIterator;
    begin
      dbcout( endl + '+++ Beginning TFunctionDictionary.debug()...', true );
      it := TFunctionDictionaryIterator.create( self );

      repeat
        if( nil <> it.value() ) then
          it.value().debug()
        ;

        it.incr();
      until ( nil = it.value() );

      it.Free();

      dbcout( '--- Done debugging TFunctionDictionary.' + endl, true );
    end
  ;

  /// Shows each functions name and reference count
  procedure TFunctionDictionary.showRefCounts();
    var
      it: TFunctionDictionaryIterator;
      s: string;
    begin
      dbcout( '+++ Beginning TFunctionDictionary reference counts...', true );
       it := TFunctionDictionaryIterator.create( self );

      repeat
        if( nil <> it.value() ) then
          begin
            s := it.value().fn.name + ': used ' + intToStr( it.value().refCounter ) + ' times.';
            dbcout( s, true );
          end
        ;

        it.incr();
      until ( nil = it.value() );

      it.Free();

      dbcout( '--- Done with TFunctionDictionary reference counts.' + endl, true );
    end
  ;

  {*
    Checks on the status of a chart function named fnName
    @param fnName name of a chart function
    @return True if the function exists and is flagged as updated. False can
    mean either an item named fnName was not found or it was found but not flagged as updated ...
  }
  function TFunctionDictionary.functionExistsAndIsUpdated( const fnName: string ): boolean;
    begin
      if( self.contains( fnName ) ) then
        result := self.value( fnName ).updated
      else
        result := false
      ;
    end
  ;

  /// Returns true if any item in the fucntion dictionary is flagged as updated else false
  function TFunctionDictionary.getUpdated(): boolean;
    var
      it: TFunctionDictionaryIterator;
    begin
      result := false;

      it := TFunctionDictionaryIterator.create( self );

      repeat
        if( nil <> it.value() ) then
          begin
            if( it.value().updated ) then
              begin
                result := true;
                break;
              end
            ;
          end
        ;

        it.incr();
      until ( nil = it.value() );

      it.Free();
    end
  ;


  {$IFDEF DATABASE_ENABLED}

  {*
    Updates the database from the dictionary with any chart function additions, removals, or modifications
    @param db object accessing the NAADSM parameters database
    @return true if all updates suceeded else false
  }
  function TFunctionDictionary.populateDatabase( db: TSqlDatabase ): boolean;
    var
      it: TFunctionDictionaryIterator;
      fnID: integer;
      success: boolean;
      tmp: integer;
      fn: TChartFunction;
    begin
      result := true; // until shown otherwise
      success := true; // until shown otherwise

      it := TFunctionDictionaryIterator.create( self );

      repeat
        if( nil <> it.value() ) then
          begin
            if( not it.value().removed ) then
              begin
                fn := it.value().fn;
                tmp := fn.populateDatabase( db, not it.value().new );
                success := ( 0 < tmp );
                it.value().new := false;
                it.value().modified := false;
              end
            else
              begin
                //Remove this item from the database
                fnID := it.value().fn.id;

                if( 0 < fnID ) then
                  begin
                    if( nil <> _sim ) then
                      success := ( _sim as TSimInput ).removeDbFunction( fnID )
                    else
                      raise exception.Create( '_sim is nil in TFunctionDictionary.populateDatabase()' )
                    ;
                  end
                ;
              end;
          end;
        ;

        if( not success ) then
          result := false
        ;

        it.incr();
      until ( nil = it.value() );

      it.Free();
    end
  ;
  {$ENDIF}
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Dictionary iterator
//-----------------------------------------------------------------------------

  {*
     Used to iterate through a TFunctionDictionary object
     @return a TFunctionDictionaryItem object from the dictionary
  }
  function TFunctionDictionaryIterator.value(): TFunctionDictionaryItem;
    begin
      result := inherited value() as TFunctionDictionaryItem;
    end
  ;
//-----------------------------------------------------------------------------

end.

{*
ModelDatabase.pas
-----------------
Begin: 2009/02/10
Last revision: $Date: 2011-10-25 04:49:08 $ $Author: areeves $
Version: $Revision: 1.18 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2009 - 2010 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.


This unit defines the base class for the business object that manages the NAADSM
database file, including versioning, schema update, and synchronizing the
temporary and permanent database files.

}

  (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator
  *)

unit ModelDatabase;

interface

  uses
    // Standard Delphi units
    Forms,
    Contnrs,

    // APHI General-Purpose Delphi Library
    SqlClasses,

    // QClasses for Delphi
    QStringMaps,

    // APHI Delphi Library for Simulation Modeling
    ChartFunction
  ;


  type
  /// Reasons for schema update listed in order of severity
  TDBSchemaUpdateReason = (
    DBUpdateUnspecified,
    DBUpdateOK,           /// This is the current version of the DB schema
    DBUpdateMinorChanges, /// The reason for the database schema update has no real effect on the end user
    DBUpdateOutputsAdded, /// New outputs are generated, but existing outputs are still valid
    DBUpdateSpecChange    /// A fairly substantial change to the model specification occurred
    //DBUpdateOutputsDropped // Some existing outputs are no longer available: this would be a good reason to save a backup copy!
  );


  type
  /// Database version check result
  TDBCheckResult = (
    DBVersionUnspecified,
    DBVersionCurrent,      /// Version current
    DBVersionUpdated,      /// Version was updated to current version
    DBVersionObsolete,     /// Version is not current
    DBVersionUnrecognized  /// Version is not recognized
  );


  type
  /// Reasons for updating a database
  TVersionUpdateReason = (
    VERSOK,               /// Changes do not affect existing inputs or outputs (the database updater handles minor updates)
    VERSBug,              /// A bug was discovered in the previous version: outputs may be invalid
    VERSModelSpecChange,  /// The model specification has changed: inputs and/or outputs are invalid
      // (This shouldn't happen without at least a change in minor app version and input translation)
    VERSUnrecognized      /// Some other program mucked around with things.
  );


  type
  /// SQL actions
  TDBUpdateActionType = (
      MDBAForceInsert,
      MDBAForceUpdate,
      MDBAAuto
  );


  type
  /// Structure to hold the names of a table and table field, specifically for PDFs and RELs
  TChartField = class
    protected
      _table: string;   /// table name
      _field: string;   /// field name

    public
      constructor create( const table, field: string );
      destructor destroy; override;

      property table: string read _table;  /// read-only access to _table
      property field: string read _field;  /// read-only access to _field
    end
  ;


  type
  /// Container structure for a list of ChartField objects
  TChartFieldList = class( TObjectList )
    public
      constructor create();
      destructor destroy(); override;
      function at( const i: integer ): TChartField;
    end
  ;


  type
  /// Base class to interact with a NAADSM model SQL database
  TModelDatabase = class( TSqlDatabase )
    protected
      // File management
      //----------------
      _workingDBFileName: string;   /// fullpath name to temporary database
      _permanentDBFileName: string; /// fullpath name to the user's scenario database
      _workingDBHasChanged: boolean;/// flag indicating changes have been made to the data
      _dbSaved: boolean;            /// indicates whether a permanent database exists; true: File|Open, false: File|New
      _isReadOnly: boolean;         /// write-status of the database file, true if the file is read-only, else false
     
      // PDF/REL handling
      //-----------------
      _chartFields: TChartFieldList; /// for managing PDFs and RELs

      // Database schema version checking
      //---------------------------------
      _vNumber: string;  /// database version number
      _vID: integer;     /// database version ID
      
      // Other junk
      //-----------
      _frmMain: TForm;  /// model's main form object

      // File management
      //----------------
      function getWorkingDBFileName(): string;
      function getPermanentDBFileName(): string;
      function getWorkingDBHasChanged(): boolean;
      function getDBSaved(): boolean;
      procedure setDBSaved( val: boolean );
      procedure setPermanentDBFileName( const val: string );
      procedure setWorkingDBHasChanged( val: boolean ); virtual;
      function getLanguageCode(): string; virtual;
      function executeWithoutChange( q: string ): boolean; virtual;

      // Schema management
      //------------------
      procedure setSchemaVersion( const versNumber, versApp, versDate, versURL, versID: string );

      /// To be used by descendant classes when creating a new database from scratch.
      procedure makeDbTables(); virtual; abstract;

      /// To be used by descendant classes to retrieve the current database schema version
      function getCurrentDbSchemaVersion(): string; virtual; abstract;
      
      /// To be used by descendant classes to check that version number (_vNumber) and version ID (_vID) are appropriate and correspond to each other.
      function schemaIDOK(): boolean; virtual; abstract;

      { Is the version of the newly opened database obsolete (i.e., will this application not attempt to update it)? }
      function versionObsolete(): boolean; virtual;

      /// To be used by descendant classes to determin why was the database schema updated: it wasn't, minor changes, new outputs, or new specification?
      function setUpdateReason(): TDBSchemaUpdateReason; virtual; abstract;


      // Simulation management
      //----------------------
      /// To be used by descendant classes to populate input parameter PDFs and RELs from the database
      procedure fillChartFields(); virtual; abstract;
      /// To be used by descendant classes to check if the simulation has completed
      function getSimulationComplete(): boolean; virtual; abstract;
      function getContainsOutput(): boolean; virtual;

    public
      // Construction/initialization/destruction
      //----------------------------------------
      { WARNING: it is up to the interface to check for the existence of the permanent database file. }
      constructor create(
        const fileName: string;
        const dbAction: TSqlDBAction;
        errMsg: PChar = nil;
        parent: TForm = nil
      ); virtual;

      destructor destroy(); override;

      // Database functionality
      //-----------------------
      function execute( q: string ): boolean; override;

      // File management
      //----------------
      function save( newFileName: string = '' ): boolean; virtual;

      // New database creation
      //----------------------
      { *Should* return true on success, false on failure. }
      function processDDLCommands( const resFileContents: string ): boolean;

      // PDF/REL handling
      //-----------------
      {* Returns the ID number of an identical chart with the same field.  If none is found, returns -1 }
      // FIX ME: AR - I don't think this version works as it should.  See comments below.
      //function compareFunctions( fieldName:String; Chart:TChartFunction; var nFunctionsWithSameName: integer ): integer;
      function removeChartFunction( const id: integer ): boolean;

      // Schema management
      //------------------
      { Can be used to alert a user if a database schema update is necessary. }
      function checkVersion( var updateReason: TDBSchemaUpdateReason ): TDBCheckResult; virtual; // Should rarely be overridden.
      /// To be used by descendant classes to indicate whether schema versionID neeeds to be updated and the reason why
      function versionUpdateReason( versionID: pstring = nil ): TVersionUpdateReason; virtual; abstract;

      { Update database to the latest schema. Return true if schema update was actually required.
        Set updateSuccess to true/false if update was attempted and was successful/unsuccessful. }
      function updateSchema( var updateSuccess: boolean ): boolean; virtual;

      // Output management
      //------------------
      procedure setRngSettings( const useFixedRandomSeed: boolean;  const randomSeed: integer );
      function completedIterations(): integer;
      function quickInsert( dict: TQStringVariantMap; keyField: string = '' ): integer; override;
      procedure quickUpdate( indexValue: variant; fieldName: string; newValue: variant ); override;
      procedure initializeAllOutputRecords(); virtual;
      procedure recordStartTime( const versionNumber: string ); virtual;
      procedure incrementCompleteIterations(); virtual;
      procedure recordEndTime(); virtual;


      // File management properties
      //---------------------------
      property workingDBFileName: string read getWorkingDBFileName;  /// read-only access to _workingDBFileName
      property permanentDBFileName: string read getPermanentDBFileName write setPermanentDBFileName;  /// access to _permanentDBFileName
      property workingDBHasChanged: boolean read getWorkingDBHasChanged write setWorkingDBHasChanged; /// access to _workingDBHasChanged
      property dbSaved: boolean read getDBSaved write setDBSaved; /// access to _dbSaved
      property isReadOnly: boolean read _isReadOnly;              /// access to _isReadOnly
      property languageCode: string read getLanguageCode;         /// returns 'en'

      // Simulation management properties
      //---------------------------------
      /// returns true or false, implementation is carried out in descendant classes
      property simulationComplete: boolean read getSimulationComplete;
       /// true if at least one iteration has been completed, else false
      property containsOutput: boolean read getContainsOutput;        


      // Schema management properties
      //-----------------------------
      /// returns the version of the current database schema, implementation carried out in descendant classes
      property currentDbSchemaVersion: string read getCurrentDBSchemaVersion;
    end
  ;
  
  
  
implementation

  uses
    SysUtils,
    StrUtils,
    Variants,

    WindowsUtils,
    NetworkUtils,
    DebugWindow,
    MyStrUtils,
    CStringList,

    FunctionDictionary
  ;


  const
    DBSHOWMSG: boolean = false; /// Set to true to enable debugging messages for this unit


  {*
    Creates a ChartField object
    @param table database table name
    @param field name of a field in table
  }
  constructor TChartField.create( const table, field: string );
    begin
      inherited create();
      _table := table;
      _field := field;
    end
  ;

  /// destroys the object and frees memory
  destructor TChartField.destroy;
    begin
      inherited destroy();
    end
  ;

  /// Creates an empty list
  constructor TChartFieldList.create();
    begin
      inherited create( true );
    end
  ;


  /// destroys the object and frees memory
  destructor TChartFieldList.destroy();
    begin
      inherited destroy();
    end
  ;


  {*
    Returns the ChartField object at position i in the list
    @param i list index value
    @return ChartField object at position i
  }
  function TChartFieldList.at( const i: integer ): TChartField;
    begin
      result := self[i] as TChartField;
    end
  ;


//-----------------------------------------------------------------------------
// Construction/destruction
//-----------------------------------------------------------------------------

  {*
    Creates the model database business object and temporary database file
    @param fileName full path to permanent database file
    @param dbAction indicates if an existing permanent database is read from or if a new one is to be created
    @param errMsg currently not used
    @param parent sets _frmMain to parent
    @comment The abstract method fillChartFields() is called, but it is upto descendant classes to implement.
  }
  constructor TModelDatabase.create( const fileName: string; const dbAction: TSqlDBAction; errMsg: PChar = nil; parent: TForm = nil  );
    var
      copied: boolean;
      tmp: string;
    begin
      _frmMain := parent;

      // Network drives are unreliable.
      // Make sure that the working file is saved to a local drive.
      if( isNetworkDrive( currentDir() ) ) then
        tmp := tempFileName( tempDir() )
      else
        tmp := tempFileName( currentDir() )
      ;

      if( dbAction = DBCreate ) then
        begin
          _permanentDBFileName := fileName;
          _workingDBFileName := directory( tmp ) + '$$$' + shortFileName( tmp );
          workingDBHasChanged := true;
          _dbSaved := false;

          // Create the working database, not the permanent one.
          inherited create( DBMSAccess, _workingDBFileName, DBCreate );

          _isReadOnly := false;

          makeDBTables();

        end
      else if( dbAction = DBOpen ) then
        begin
          // Store some file names
          _permanentDBFileName := fileName;
          _workingDBFileName := directory( tmp ) + '$$$' + shortFileName( tmp );
          workingDBHasChanged := false;
          _dbSaved := true;

          _isReadOnly := fileIsReadOnly( fileName );

          // Make a temporary copy of the selected database
          copied := windowsCopyFile( _permanentDBFileName, _workingDBFileName );
          if( not copied ) then
            raise exception.create( 'TSMDatabase.create: Could not copy the selected scenario file.' )
          ;

          // Open the copy
          inherited create( DBMSAccess, _workingDBFileName, DBOpen );

          if( not _errorOnOpen ) then
            begin
              // Compact the database
              compact( true ); // it's pointless to raise an exception if compact fails.

              if( not(isOpen) ) then open();
            end
          else
            _isOpen := false
          ;
        end
      ;

      _chartFields := TChartFieldList.create();
      fillChartFields();
    end
  ;


  {*
     Closes connection to the temporary database, deletes this file, and frees itself from memory.
  }
  destructor TModelDatabase.destroy();
    begin
      // FIX ME: for some reason, after an exception occurs, the temporary file just won't go away.
      close();

      if( 0 <> length( self.workingDBFileName ) ) then
        begin
          dbcout( 'Deleting ' + self.workingDBFileName, DBSHOWMSG );

          setFileReadWrite( self.workingDBFileName );

          if( fileExists( self.workingDBFileName ) ) then
            begin
              if( not deleteFile( self.workingDBFileName ) ) then
                {$IFDEF DEBUG}
                  msgOK( ansiReplaceStr( tr( 'The temporary file xyz could not be deleted.' ), 'xyz', self.workingDBFileName ) )
                {$ENDIF}
              ;
            end
          ;
        end
      ;

      _chartFields.free();
      
      inherited destroy();
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// New database creation
//-----------------------------------------------------------------------------

  {*
    Executes a series of SQL Data Definition Language commands to create a new database
    @param resFileContents contents from a resource file holding DDL commands
    @comment the DDL files of schema versions are stored in the projects's sm_database folder
    @return true if all the commands were executed else false
  }
  function TModelDatabase.processDDLCommands( const resFileContents: string ): boolean;
    var
      str: string;
      startingList: TCStringList;
      ddlList: TCStringList;
      ptr: pchar;
      testStr: string;
      cmd: PChar;
  	begin
      result := true;

      startingList := TCStringList.create();
      startingList.text := resFileContents;

      // Remove SQL comments and empty lines from the starting list
      ptr := startingList.first();

      while ptr <> nil do
        begin
          testStr := trim( ptr );

          if( leftStr( testStr, 2 ) = '--' ) or( length( testStr ) = 0 ) then
            begin
              startingList.remove();
              ptr := startingList.current();
            end
          else
            ptr := startingList.next()
          ;

        end
      ;

      // Explode remaining strings at each semicolon to produce the DDL table and key creation commands
      str := startingList.text;
      ddlList := TCStringList.create( str, ';', true );

      //if( DBSMDATABASE ) then ddlList.debug( true );

      // Execute the commands!
      cmd := ddlList.first();

      while( cmd <> nil ) do
      	begin
        	str := trim( cmd );

          if( length( str ) > 0 ) then
            begin
              if( not( self.execute( str ) ) ) then
                begin
                  result := false;
                  break;
                end
              ;
            end
          ;

          cmd := ddlList.next();
        end
      ;

      // clean up
      startingList.Clear();
      ddlList.Clear();
      freeAndNil( startingList );
      freeAndNil( ddlList );
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Database modification
//-----------------------------------------------------------------------------

  {*
    Saves a copy of the working database as the current permanent database
    @param newFileName if specified, essentialy a "save as" naming the permanent
    database newFileName. If empty string then essentially a "save", the former
    permanent database will be overwritten.
    @return true if operation succeeded, else false
  }
  function TModelDatabase.save( newFileName: string = '' ): boolean;
    begin
      if( 0 = length( newFileName ) ) then
        newFileName := permanentDBFileName
      ;

      if( ( 0 = length( workingDBFileName )  ) or ( 0 = length( newFileName ) ) ) then
        begin
          raise exception.create( 'Missing file name in SMDatabase.save()' );
          result := false;
          exit;
        end
      ;

      dbcout( 'TSMDatabase.save()', DBSHOWMSG );
      
      waitForExecuteComplete(); // This seems to be worthless.
      close(); // Closing the database seems to buy enough time for the last operation to finish before a copy of the file is made.

      if( windowsCopyFile( self.workingDBFileName, newFileName ) ) then
        begin
          self.dbSaved := true;
          self.workingDBHasChanged := false;
          self.permanentDBFileName := newFileName;

          if( _isReadOnly ) then // The database must be closed and reopened to allow read/write privileges
            begin
              close();

              setFileReadWrite( newFileName );
              setFileReadWrite( self.workingDBFileName );
              _isReadOnly := false;

              open();
            end
          ;

          result := true;
        end
      else
        begin
          self.workingDBHasChanged := true;
          result := false;
        end
      ;

      open();
    end
  ;
  

  {*
     Executes SQL query q
     @param q SQL query command
     @return true if the command was executed, false if the command failed or the database is read-only
  }
  function TModelDatabase.execute( q: string ): boolean;
    begin
      // If the file/database is read-only, don't even attempt to make changes.
      
      if( _isReadOnly ) then
        result := false
      else
        begin
          workingDBHasChanged := true;
          result := inherited execute( q );
        end
      ;
    end
  ;


  function TModelDatabase.executeWithoutChange( q: string ): boolean;
    begin
      if( _isReadOnly ) then
        result := false
      else
        result := inherited execute( q )
      ;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// File management
//-----------------------------------------------------------------------------
  /// Get function for property workingDBFileName
  function TModelDatabase.getWorkingDBFileName(): string; begin result := _workingDBFileName end;

  /// Get function for property permanentDBFileName
  function TModelDatabase.getPermanentDBFileName(): string; begin result := _permanentDBFileName end;

  /// Get function for property workingDBHasChanged
  function TModelDatabase.getWorkingDBHasChanged(): boolean; begin result := _workingDBHasChanged end;

  /// Get function for property dbSaved
  function TModelDatabase.getDBSaved(): boolean; begin result := _dbSaved end;

  /// Set function for property dbSaved, setting _dbSaved to val
  procedure TModelDatabase.setDBSaved( val: boolean ); begin _dbSaved := val; end;

  /// Set function for property workingDBHasChanged, setting _workingDBHasChanged to val
  procedure TModelDatabase.setWorkingDBHasChanged( val: boolean );
    begin
      _workingDBHasChanged := val;
    end
  ;

  /// Set function for property permanentDBFileName, setting _permanentDBFileName to val
  procedure TModelDatabase.setPermanentDBFileName( const val: string ); begin _permanentDBFileName := val; end;

  /// Get function for property languageCode
  function TModelDatabase.getLanguageCode(): string; begin result := 'en'; end;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Schema management
//-----------------------------------------------------------------------------
  {*
    Updates fields of the database table DBSchemaVersion providing metadata about the schema
    @param versNumber version number
    @param versApp  not used yet and and now always set to 'NAADSMXXXX'
    @param versDate version date
    @param versURL version URL, this field not in the DB yet and therefore should be ''
    @param versID  version ID
    @comment see project's \sm_database\SMDatabase.pas
  }
  procedure TModelDatabase.setSchemaVersion( const versNumber, versApp, versDate, versURL, versID: string );
    var
      vDict: TQueryDictionary;
    begin

      vDict := TQueryDictionary.Create();

      execute( 'DELETE FROM `DBSchemaVersion`' );

      vDict['VersionNumber'] := sqlQuote( versNumber );
      vDict['VersionApplication'] := sqlQuote( versApp );
      vDict['VersionDate'] := sqlQuote( versDate );
      vDict['VersionInfoURL'] := sqlQuote( versURL );
      vDict['VersionID'] := versID;

      execute( writeQuery( 'DBSchemaVersion', QInsert, vdict ) );

      vDict.Free();
    end
  ;

  {*
    Is the version of the newly opened database obsolete (i.e., will this application not attempt to update it)?
    @return false only
  }
  function TModelDatabase.versionObsolete(): boolean;
    begin
      result := false;
    end
  ;


  {*
     Compares the database version information with what is the current schema version
     @param updateReason an output parameter, the value can be used by the caller to alert user if a schema update is necessary
     @return enumration describing the database version status - current, obsolete, etc.
  }
  function TModelDatabase.checkVersion( var updateReason: TDBSchemaUpdateReason ): TDBCheckResult;
  	var
      row: TSqlRow;
  	begin
      if( nil = _sqlResult ) then
        createSqlResult()
      ;

      _sqlResult.runQuery( 'SELECT `versionNumber`, `VersionID` FROM `DBSchemaVersion`' );

      if( not( _sqlResult.success ) ) then
      	result := DBVersionUnrecognized
      else if( 1 <> _sqlResult.numRows ) then
      	result := DBVersionUnrecognized
      else
      	begin
          row := _sqlResult.fetchArrayFirst();
          _vNumber := fixup( row.field('VersionNumber') );
          _vID := row.field( 'VersionID' );

          if( not( schemaIDOK() ) ) then
						result := DBVersionUnrecognized
          else
         		begin
              if( versionObsolete() ) then
              	result := DBVersionObsolete
              else if( currentDbSchemaVersion = _vNumber ) then
                begin
                  updateReason := DBUpdateOK;
                  result := DBVersionCurrent;
                end
              else
                begin
                	// This database can be updated.
                  updateReason := setUpdateReason();
                  result := DBVersionUpdated;
                end
              ;
            end
          ;
        end
      ;
    end
  ;


  {*
     Update database to the latest schema. 
     Set updateSuccess to true/false if update was attempted and was successful/unsuccessful.
     @return true if schema update was actually required but hardcoded to always return false.
  }
  function TModelDatabase.updateSchema( var updateSuccess: boolean ): boolean;
    begin
      updateSuccess := true;
      result := false;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// PDF/REL handling
//-----------------------------------------------------------------------------
  {* Returns the ID number of an identical chart with the same field.  If none is found, returns -1 }
  // FIX ME: I don't think this function is quite right.
  // Functions with the same name can occur regardless of the field.  This function doesn't seem to
  // check for that.
  //
  // It may be irrelvant, once functionDictionary has checkAndInsert().
  (*
  function TModelDatabase.compareFunctions( fieldName: String; Chart: TChartFunction; var nFunctionsWithSameName: integer ): Integer;
    var
      ret_val:Integer;
      q: string;
      db2: TSqlDatabase;
      res, res2: TSqlResult;
      row: TSqlRow;
      _chart: TChartFunction;
      I, Count, chartID:Integer;
      chartName: string;

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
      ret_val := -1;

      q := 'SELECT '
        + 'inChart.chartID, inChart.fieldName '
        //+ 'FROM inChart WHERE LEFT( inChart.fieldName, ' + intToStr( length( fieldName ) ) + ') = ' + SQLQuote( fieldName )
        + 'FROM inChart WHERE inChart.fieldName = ' + SQLQuote( fieldName )
      ;

      db2 := self as TSqlDatabase;
      res := TSqlResult.create( q, db2 );

      Count := res.numRows;

      // Are there other charts with a similar name?  If so, then the
      // name suffix needs to be changed (e.g. from "Distance distribution" to "Distance distribution (1)"
      if ( 0 < Count ) then
        begin
          chartName := stripParens( chart.name );
          res2 := TSqlResult.create(
            'SELECT inChart.chartID FROM inChart WHERE LEFT( inChart.chartName, ' + intToStr( length( chartName ) ) + ' ) = ' + sqlQuote( chartName ),
            db2
          );

          nFunctionsWithSameName := res2.numRows;
          res2.free();

          row := res.fetchArrayFirst();

          for I := 1 to Count do
            begin
              if ( row.field( 'chartID' ) ) then
                begin
                  chartID :=  StrToInt( row.field( 'chartID' ) );
                  _chart := createFunctionFromDB( self, chartID );

                  if ( Chart.compare( _chart ) ) then
                    begin
                      ret_val := chartID;
                      _chart.free();
                      break;
                    end;
                  _chart.free();
                end;

              row := res.fetchArrayNext();
            end
          ;
        end
      ;

      freeAndNil( res );

      result := ret_val;
    end
  ;
*)  

  {*
    Updates the database when a chart function is removed, deleting all input parameter references to it and the chart data itself.
    @param id identifier for the chart, the primary key value in table inChart
    @return 
  }
  function TModelDatabase.removeChartFunction( const id: integer ): boolean;
    var
      i: integer;
      q: string;
      cf: TChartField;
    begin
      result := true; // Until shown otherwise.

      // Update any/all fields that might use this chart.
      //-------------------------------------------------
      for i := 0 to _chartFields.Count - 1 do
        begin
          cf := _chartFields.at( i );

          if( self.fieldExists( cf.field, cf.table ) ) then
            begin
              q := 'UPDATE `' + cf.table + '` SET `' + cf.field + '` = NULL WHERE `' + cf.field + '` = ' + intToStr( id );
              result := result and execute( q );
            end
          ;
        end
      ;

      // Remove the chart points
      q := 'DELETE FROM `inChartDetail` WHERE `chartID` = ' + intToStr( id );
      result := result and execute( q );

      // Finally, delete the chart itself
      q := 'DELETE FROM `inChart` WHERE `chartID` = ' + intToStr( id );
      result := result and execute( q );
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Output management
//-----------------------------------------------------------------------------

  {*
    Updates the database with settings for the simulation random seed
    @param useFixedRandomSeed true to use random seed else false and the system will gnerate a seed
    @param randomSeed value to be used as the seed for the random number generator
  }
  procedure TModelDatabase.setRngSettings( const useFixedRandomSeed: boolean;  const randomSeed: integer );
    var
      q: string;
    begin
      q :=
        'UPDATE `inGeneral` SET '
        + '  `useFixedRandomSeed` = ' + usBoolToText( useFixedRandomSeed )
        + ', `randomSeed` = ' + intToStr( randomSeed )
      ;

      execute( q );
    end
  ;

  {*
    Get function for property containsOutput
    @return true if at least one iteration has been completed, else false
  }
  function TModelDatabase.getContainsOutput(): boolean;
    begin
      result := ( 0 < completedIterations() );
    end
  ;


  {*
    Helper function that checks the database for the number of iterations completed
    @return the number of completed iterations
  }
  function TModelDatabase.completedIterations(): integer;
    var
      row: TSqlRow;
    begin
      if( nil = _sqlResult ) then
        createSqlResult()
      ;

      _sqlResult.runQuery( 'SELECT completedIterations FROM outGeneral' );
      row := _sqlResult.fetchArrayFirst();

      if( nil = row ) then
        result := 0
      else if( null = row.field('completedIterations') ) then
        result := 0
      else
        result := row.field('completedIterations')
      ;
    end
  ;

  {*
    Update to a single record field value is performed by seeking the specified index
    @param indexValue Index value of table key being sought
    @param fieldName Name of field to be updated
    @param newValue Updated value
  }
  procedure TModelDatabase.quickUpdate( indexValue: variant; fieldName: string; newValue: variant );
    begin
      workingDBHasChanged := true;
      inherited quickUpdate( indexValue, fieldName, newValue );
    end
  ;

  {*
    Inserts a single record into a table of the database
    @param dict object holding the values for the record fields
    @param keyField field name of table's primary key, optional
    @return Returns the value for keyField or an empty string if keyField is an empty string
  }
  function TModelDatabase.quickInsert( dict: TQStringVariantMap; keyField: string = '' ): integer;
    begin
      workingDBHasChanged := true;
      result := inherited quickInsert( dict, keyField );
    end
  ;


  {*
     Removes any existing daily and iteration output data from the database
  }
  procedure TModelDatabase.initializeAllOutputRecords();
    var
      q: string;
    begin
    	// Clear status times and other general properties
      //------------------------------------------------
      q := 'UPDATE `outGeneral` SET'
      	+ ' `simulationStartTime` = NULL,'
        + ' `simulationEndTime` = NULL,'
        + ' `completedIterations` = 0,'
        + ' `version` = NULL'
      ;
			execute( q );

      if( self.tableExists( 'outIteration' ) ) then
        execute( 'DELETE FROM `outIteration`' )
      ;
      if( self.tableExists( 'outDaily' ) ) then
        execute( 'DELETE FROM `outDaily`' )
      ;
    end
  ;

  {*
    Sets the simulation start time and version number in the database
    @param versionNumber schema version number
  }
  procedure TModelDatabase.recordStartTime( const versionNumber: string );
    var
      q: string;
    begin
      q := 'UPDATE `outGeneral` SET `simulationStartTime` = NOW(), `version` = "' + versionNumber + '"';
    	execute( q );
    end
  ;


  /// Updates the database with the number of completed iterations
  procedure TModelDatabase.incrementCompleteIterations();
    begin
      execute( 'UPDATE `outGeneral` SET `completedIterations` = `completedIterations` + 1' );
    end
  ;

  /// Updates the database with the simulation end time
  procedure TModelDatabase.recordEndTime();
  	var
    	q: string;
  	begin
      q := 'UPDATE `outGeneral` SET `simulationEndTime` = NOW()';
    	execute( q );
    end
  ;
//-----------------------------------------------------------------------------


end.
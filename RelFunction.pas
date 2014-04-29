{*
RelFunction.pas
----------------
Begin: 2005/02/26
Last revision: $Date: 2011-01-26 04:38:04 $ $Author: areeves $
Version number: $Revision: 1.24 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2005 - 2011 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

----------------------------------------------------
Relational functions (Rels) are used where one variable is a function of another.
For example, the probability of detecting an infectious herd is a function of
time since the herd was infected. This unit defines the class for creating and
managing relational functions.
}

  (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator

    Base classes have been documented pretty well, the derived classes not so much
  *)

unit RelFunction;

{$INCLUDE Defs.inc}

interface

  uses
    ChartFunction,
    SysUtils,

    {$IFDEF DATABASE_ENABLED}
    SqlClasses,
    {$ENDIF}

    Sdew,
    Points
  ;

  const
    /// An extremely high (2000) Y-axis value to initialize an unlimited resource
    UNLIMITEDRELFUNCTION_Y: integer = 2000 
  ;


  type
  /// Class for creating and managing relational function charts
  TRelFunction = class( TChartFunction )
    protected
     _points: RPointArray;  /// to hold the point (X,Y) coordinates of the Rel

      function getPointCount(): integer;
      function getPoint( i: integer ): RPoint;
      procedure setPoint( i: integer; pt: RPoint );
      function getPointArray(): RPointArray;

      function getYStartsAtZero(): boolean; override;

      function getDescr(): string; override;

    public
      constructor create(); overload; override;
      constructor create( const src: TRelFunction ); reintroduce; overload;
      constructor create( xUnits, yUnits: TChartUnitType; unlimited: boolean = false ); overload;
      constructor create( pnts: RPointArray; xUnits: TChartUnitType; yUnits: TChartUnitType ); overload;
      constructor create( pnts: T2DPointList; xUnits: TChartUnitType; yUnits: TChartUnitType ); overload;
      {$IFDEF DATABASE_ENABLED}
      constructor create( db: TSqlDatabase; chartID: integer; xu: TChartUnitType; yu: TChartUnitType ); overload;
      {$ENDIF}
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ): boolean; override;

      {Converts y values from proportions and probabilities (0 to 1) to percentages (0 to 100) and vice versa }
      function convertYUnitsTo( const newUnits: TChartUnitType ): boolean;

      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createRecordPointArray(): RPointArray;

      function ssXml( const indent: integer ): string; override;

      { Generates a plain-text representation of the chart. }
      function asCsv(): string; override;

      function createCopy(): TChartFunction; override;

      procedure addPoint( x, y: double );

      function y( const x: double ): double;

      procedure debug(); override;

      property pointCount: integer read getPointCount; /// Number of points in the Rel, read-only
      property points[i: integer]: RPoint read getPoint write setPoint; /// Access to a value in _points
      //property pointArray: RPointArray read getPointArray;
    end
  ;

  // Global helper functions
  //-----------------------------------------------
  function createRelFromXml( xml: string ): TRelFunction; overload;
  function createRelFromXml( element: pointer; sdew: TSdew ): TRelFunction; overload;
  function createRelFromXml( Element: pointer; Sdew:TSdew; extData: Pointer ): TObject; overload;

  // FIX ME: I'm not sure if this function is still necessary.
  function createRelFromBaseObject( rel: TRelFunction ): TChartFunction;

implementation

  uses
    Math,
    StrUtils,
    
    MyStrUtils,
    DebugWindow,
    I88n

    {$IFDEF DATABASE_ENABLED}
    , FunctionEnums
    {$ENDIF}
  ;

  const
    DBRELFUNCTION: boolean = false; /// set to true to enable debugging messages for this unit.


//-----------------------------------------------
// Global helper functions
//-----------------------------------------------
  // FIX ME: I'm not sure if this function is still necessary. rbh20101026: Still used by FrameFunctionEditor
  {*
    Creates an RelFunction object from rel
    @param Rel object from which properties are copied
    @return ChartFunction object, the base class for Rels and PDFs
    @comment This function and similar one in the pdf class allows the caller to
    create an Rel or a PDF as necessary, since the ancestor class is returned.
  }
  function createRelFromBaseObject( rel: TRelFunction ): TChartFunction;
    var
      newobj: TChartFunction;
    begin
      newObj := TRelFunction.create( rel );
      newobj.name := rel.name;
      newobj.id := rel.id;
      newobj.dbField := rel.dbField;

      result := newobj;
    end
  ;
//-----------------------------------------------

  /// Creates an empty Rel object
  constructor TRelFunction.create();
    begin
      inherited create();
      setChartType( CTRel );

      xUnits := UUnknown;
      yUnits := UUnknown;
    end
  ;

  {*
    Creates an Rel object having units for the x and y axes but not necessarily data points.
    @param xUnits enumeration delimiting units for the x-axis
    @param yUnits enumeration delimiting units for the y-axis
    @param unlimited if true two points are added having y-axis values of UNLIMITEDRELFUNCTION_Y, else the Rel has no data.
  }
  constructor TRelFunction.create( xUnits, yUnits: TChartUnitType; unlimited: boolean = false );
    begin
      inherited create();
      setChartType( CTRel );

      _xUnits := xUnits;
      _yUnits := yUnits;
    
      if( unlimited ) then
        begin
          setLength( _points, 0 );
          // FIX ME: can SharcSpread get by with a single point?
          addPoint( 0, UNLIMITEDRELFUNCTION_Y );
          addPoint( 1, UNLIMITEDRELFUNCTION_Y );
        end
      ;
      (*
      // Hmm.  Should this really be the default action?
      else
        begin
          setLength( _points, 2 );
          addPoint( 0, 0 );
          addPoint( 1, 10 );
        end
      ;
      *)
    end
  ;


  {*
    Creates a copy of src
    @param src source RelFunction object for copying
  }
  constructor TRelFunction.create( const src: TRelFunction );
    begin
      dbcout( '*** Calling TRelFunction copy constructor.', DBRELFUNCTION );
      
      inherited create( src );

      setChartType( CTRel );

      setLength( _points, length( src._points ) );
      _points := copy( src._points, 0, length( src._points ) );

      _xUnits := src._xUnits;
      _yUnits := src._yUnits;

      _name := src._name;
      _id := src._id;
      _dbField := src._dbField;
    end
  ;


  {*
    Creates an Rel object having units for the x and y axes and data points.
    @param pnts a RPointArray data structure of x and y axis coordinate data values
    @param xUnits enumeration delimiting units for the x-axis
    @param yUnits enumeration delimiting units for the y-axis
  }
  constructor TRelFunction.create( pnts: RPointArray; xUnits: TChartUnitType; yUnits: TChartUnitType );
    begin
      inherited create();
      setChartType( CTRel );

      setLength( _points, length( pnts ) );
      _points := copy( pnts, 0, length( pnts ) );

      _xUnits := xUnits;
      _yUnits := yUnits;
    end
  ;


  {*
    Creates an Rel object having units for the x and y axes and data points.
    @param pnts a T2DPointList data structure of x and y axis coordinate data values
    @param xUnits enumeration delimiting units for the x-axis
    @param yUnits enumeration delimiting units for the y-axis
  }
  constructor TRelFunction.create( pnts: T2DPointList; xUnits: TChartUnitType; yUnits: TChartUnitType );
    var
      i: integer;
    begin
      inherited create();
      setChartType( CTRel );
      _xUnits := xUnits;
      _yUnits := yUnits;

      setLength( _points, pnts.Count );

      for i := 0 to pnts.Count - 1 do
        begin
          _points[i].x := pnts[i].x;
          _points[i].y := pnts[i].y;
        end
      ;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  
  {*
    Creates an Rel object having units and fills it with x,y coordinates from the database
    @param db object representing the database
    @param chartID key value for the Rel chart data in the table inChartDetail
    @param xu enumeration delimiting units for the x-axis
    @param yu enumeration delimiting units for the y-axis
  }
  constructor TRelFunction.create( db: TSqlDatabase; chartID: integer; xu: TChartUnitType; yu: TChartUnitType );
    var
      q: string;
      res: TSqlResult;
      row: TSqlRow;
    begin
      inherited create();

      q := 'SELECT [x], [y] FROM [inChartDetail] WHERE [chartID] = ' + intToStr( chartID ) + ' ORDER BY [pointOrder]';
      res := TSqlResult.create( q, db );

      row := res.fetchArrayFirst();

      while( row <> nil ) do
        begin
          addPoint( double( row.field(0) ), double( row.field(1) ) );
          row := res.fetchArrayNext();
        end
      ;

      setXUnits( xu );
      setYUnits( yu );

      freeAndNil( res );
    end
  ;
  {$ENDIF}


  /// Destroys the object and frees it from memory
  destructor TRelFunction.destroy();
    begin
      dbcout( 'DESTRUCTOR Destroying TRelFunction ' + self.name, DBRELFUNCTION );
      setLength( _points, 0 );
      dbcout( 'Points are gone!', DBRELFUNCTION );
      inherited destroy();
    end
  ;

  /// Returns a copy of self.
  function TRelFunction.createCopy(): TChartFunction;
    begin
      result := TRelFunction.create( self );
    end
  ;


  {$IFDEF DATABASE_ENABLED}

  {*
    Copies the data from self to the database
    @param db object representing the database
    @param update if true performs an update query, else an insert query
    @return the key value for this record in table inChartDetail
  }
  function TRelFunction.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
    var
      q, q2: string;
      i: integer;
      dict: TQueryDictionary;
    begin
      dict := TQueryDictionary.create();

      if( update ) then
        begin
          q := 'DELETE FROM `inChartDetail` WHERE `chartID` = ' + intToStr( id );
          db.execute( q );
        end
      ;

      dict['fieldName']   := db.sqlQuote( smChartStr( TSMChart( self.dbField ) ) );
      dict['chartName']   := db.sqlQuote( self.name );
      dict['isPdf']       := usBoolToText( false );
      dict['chartType']   := 'NULL';
      dict['mean']        := 'NULL';
      dict['stddev']      := 'NULL';
      dict['min']         := 'NULL';
      dict['mode']        := 'NULL';
      dict['max']         := 'NULL';
      dict['xAxisUnits']  := db.sqlQuote( chartUnitTypeAsSqlString( xUnits ) );
      dict['yAxisUnits']  := db.sqlQuote( chartUnitTypeAsSqlString( yUnits ) );
      dict['notes']       := db.sqlQuote( notes );

      if( not update ) then
        q := writeQuery( 'inChart', QInsert, dict )
      else
        q := writeQuery( 'inChart', QUpdate, dict, 'WHERE [chartID] = ' + intToStr( id ) );
      ;

      dict.Clear();
      freeAndNil( dict );

      //dbcout( q, DBSMRELFUNCTION );

      db.execute( q );

      if( not update ) then id := db.lastInsertID();

      // Populate the inChartDetail with the individual points
      q := 'INSERT INTO [inChartDetail] ([chartID], [pointOrder], [x], [y])';

      for i := low( _points ) to high( _points ) do
        begin
          q2 := q + ' VALUES (' + usFloatToStr( id ) + ', '+ intToStr( i ) + ', ' + usFloatToStr( _points[i].x ) + ', ' + usFloatToStr( _points[i].y ) + ')';
          db.execute( q2 );
        end
      ;

      result := id;
    end
  ;
  {$ENDIF}

  {*
    Creates XML formated string that defines the relational function 
    @param  indent the number of leading spaces to indent the lines of XML
    @return XML formated string of data that delimit the Rel
  }
  function TRelFunction.ssXml( const indent: integer ): string;
    var
      i, j: integer;
      valStr: string;
      unitsStr: string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<relational-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<relational-function>' + endl
      ;

      for i := low( _points ) to high( _points ) do
        begin
          valStr := usFloatToStr( _points[i].x );
          result := result + _indent + '  <value><x>' + valStr + '</x>';

          for j := ( length( valStr ) + 8 + 9 ) to 28 do
            begin
              result := result + ' ';
            end
          ;
          result := result + '  ';

          if( yUnits in [ UPercent, UPercentProbability ] ) then
            result := result + '<y>' + usFloatToStr( _points[i].y / 100 ) + '</y></value>' + endl
          else
            result := result + '<y>' + usFloatToStr( _points[i].y ) + '</y></value>' + endl
          ;
        end
      ;

      unitsStr := ansiReplaceStr( chartUnitTypeAsSSXml( xUnits ), '<units>', '<x-units>' );
      unitsStr := ansiReplaceStr( unitsStr, '</units>', '</x-units>' );
      result := result + _indent + '  ' + unitsStr + endl;

      if( UPercent = yUnits ) then
        begin
          // The NAADSM core model expects proportions, not percentages.
          unitsStr := ansiReplaceStr( chartUnitTypeAsSSXml( UProportion ), '<units>', '<y-units>' );
        end
      else if( UPercentProbability = yUnits ) then
        begin
          // The NAADSM core model expects proportions, not percentages.
          unitsStr := ansiReplaceStr( chartUnitTypeAsSSXml( UProbability ), '<units>', '<y-units>' );
        end
      else
        unitsStr := ansiReplaceStr( chartUnitTypeAsSSXml( yUnits ), '<units>', '<y-units>' )
      ;
      unitsStr := ansiReplaceStr( unitsStr, '</units>', '</y-units>' );
      result := result + _indent + '  ' + unitsStr + endl;

      result := result + _indent + '</relational-function>' + endl;
    end
  ;


  {*
    Creates a list of x,y point coordinates that delimit the relational function
    @param  indent the number of leading spaces to indent the lines of XML
    @return list of coordinates, one pair per line separated by a comma or other locale character
    @comment the units for the x and y values are not defined
  }
  function TRelFunction.asCsv(): string;
    var
      i: integer;
    begin
      result := '';

      if( self.validate() ) then
        begin
          try
            result := 'x' + csvListSep() + ' y' + endl;

            for i := 0 to length( _points ) - 1 do
              begin
                result :=
                  result
                  + csvFloatToStr( _points[i].x )
                  + csvListSep() + ' '
                  + csvFloatToStr( _points[i].y )
                  + endl
                ;
              end
            ;
          except
            raise exception.Create( 'An exception occurred in TRelFunction.asCsv()' );
            result := '';
          end;
        end
      ;
    end
  ;

  {*
     Compares Chart to self.
     @param Chart generally another RelFunction object
     @return true if the two have the same axes units and coordinate values, else false
  }
  function TRelFunction.compare( Chart: TChartFunction ): boolean;
    var
      i, count: integer;
    begin
      result := inherited compare( chart );

      if( result ) then
        begin
          if( self.xUnits <> Chart.xUnits ) then
            result := false
          else if( self.yUnits <> Chart.yUnits ) then
            result := false
          else if( length( self._points ) <> length( (Chart as TRelFunction)._points ) ) then
            result := false
          else
            begin
              count := length( self._points );

              for i := 0 to count - 1 do
                begin
                  if( ( self._points[I].x <> (Chart as TRelFunction)._points[I].x ) or ( self._points[I].y <> (Chart as TRelFunction)._points[I].y ) ) then
                    begin
                      result := false;
                      break;
                    end
                  ;
                end
              ;
            end
          ;
        end
      ;
    end
  ;


  {*
    Add a new point to the chart
    @param x value on x-axis
    @param y value on y-axis
  }
  procedure TRelFunction.addPoint( x, y: double );
    begin
      setLength( _points, length(_points)+1 );
      _points[high(_points)].x := x;
      _points[high(_points)].y := y;
    end
  ;


  {*
     Returns a copy of the chart coordinates in an RPointArray structure
  }
  function TRelFunction.createRecordPointArray(): RPointArray;
    var
      recArray: RPointArray;
      i: integer;
    begin
      setLength( recArray, length( _points ) );

      for i := low( _points ) to high( _points ) do
        begin
          recArray[i].x := _points[i].x;
          recArray[i].y := _points[i].y;
        end
      ;

      result := recArray;
    end
  ;


  {*
     Get function for property points[], returning the x,y coordinates of chart point i
     @param i array element index
     @return reference to the contentsof element i
     @comment currently no bound checking that i is a valid array element
  }
  function TRelFunction.getPoint( i: integer ): RPoint;
    begin
      result := _points[i];
    end
  ;

  {*
     Sets the x,y coordinates of chart point i for property points[]
     @param i array element index
     @param pt structure holding the x and y coordinates of the point
     @comment currently no bound checking that i is a valid array element
  }
  procedure TRelFunction.setPoint( i: integer; pt: RPoint );
    begin
      _points[i] := pt;
    end
  ;

  /// Returns a reference to _points
  function TRelFunction.getPointArray(): RPointArray;
    begin
      result := _points;
    end
  ;

  /// Returns the number of chart coordinates
  function TRelFunction.getPointCount(): integer;
    begin
      result := length( _points );
    end
  ;


  /// Returns the text "Relational function"
  function TRelFunction.getDescr(): string;
    begin
      result := 'Relational function';
    end
  ;


  /// True if the first chart point y-axis coordinate value is 0, else false
  function TRelFunction.getYStartsAtZero(): boolean;
    begin
      if( 0 = pointCount ) then
        result := false
      else
        result := ( points[0].y = 0 )
      ;
    end
  ;

  {*
     Provides an estimate of the value of y for the input x value
     @param x value on the x-axis of the chart
     @return an estimate of the value of y corresponding to the inputted x
     @comment The result is an approximation when x is outside the chart range or does not exactly match an existing chart value.
  }
  function TRelFunction.y( const x: double ): double;
    var
      i: integer;
      p1, p2: RPoint;
      m: double; // Slope for the line segment of interest
    begin
      // Deal with easy cases first
      //---------------------------
      if( 0 = pointCount ) then
        begin
          raise exception.create( 'There are too few points to specify this function in TRelFunction.y()' );
          result := NaN;
        end
      else if( x <= points[0].x ) then
        result := points[0].y
      else if( x >= points[pointCount - 1].x ) then
        result := points[pointCount - 1].y

      // Things are a little more complicated here
      //------------------------------------------
      else
        begin
          // Determine slope for the line segment of interest
          //-------------------------------------------------
          for i := 0 to pointCount - 1 do
            begin
              if( x <= points[i].x ) then
                begin
                  p1 := points[ i - 1 ];
                  p2 := points[ i ];
                  break;
                end
              ;
            end
          ;

          if( x = p2.x ) then
            result := p2.y
          else
            begin
              // Perform the simple linear interpolation
              //----------------------------------------
              m := ( p2.y - p1.y ) / ( p2.x - p1.x );
              result := p1.y + ( m * ( x - p1.x ) );
            end
          ;
        end
      ;
    end
  ;


  {*
    Validates the chart data
    @param err validation fails then an explanation message is set in the address location of err
    @return true if there is no missing data, y values are not negative, and x values consecutively increase; else false. 
  }
  function TRelFunction.validate( err: pstring = nil ): boolean;
    var
      i: integer;
    begin
      result := true;

      if( 1 > length( _points ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'There are too few points to specify this function.' );
          result := false;
          exit;
        end
      ;

      for i := 0 to length( _points ) - 1 do
        begin

          if( isNaN( _points[i].y ) or isNaN( _points[i].x ) ) then
            begin
              result := false;
              if( nil <> err ) then err^ := tr( 'All values must be specified.' );
              break;
            end
          ;

          // FIX ME: This may be the case for NAADSM, but I don't think that it is necessarily always the case.
          // All y (probability) values must be greater than or equal to 0
          if( 0.0 > _points[i].y ) then
            begin
              result := false;
              if( nil <> err ) then err^ := tr( 'All values must be at least 0.' );
              break;
            end
          ;

          // Each X must be greater than the preceding one
          if( i > 0 ) then
            begin
              if( _points[i - 1].x >= _points[i].x ) then
                begin
                  result := false;
                  if( nil <> err ) then err^ := tr( 'Each X value must be greater than the preceding one.' );
                  break;
                end
              ;
            end
          ;
        end
      ;
    end
  ;


  {*
    Converts y values from probabilities or proportions (0 to 1) to percentages (0 to 100), and vice versa
    @param newUnits the units to which the y values should be converted
    @comment additionally the y-axis unit is set to UPercentProbability
    @return True on successful conversion
  }
  function TRelFunction.convertYUnitsTo( const newUnits: TChartUnitType ): boolean;
    var
      i: integer;
      multiplier: double;
    begin
      result := true;

      // Proportions and percentages are interchangeable.  Try the conversions here.
      if( ( UPercent = self.yUnits ) and ( UProportion = newUnits ) ) then // Convert percent to proportion: divide by 100
        multiplier := 1.0 / 100.0
      else if( ( UProportion = self.yUnits ) and ( UPercent = newUnits ) ) then // Convert proportion to percent: multiply by 100
        multiplier := 100.0
      else if( ( UPercentProbability = self.yUnits ) and ( UProbability = newUnits ) ) then // Convert percent to proportion: divide by 100
        multiplier := 1.0 / 100.0
      else if( ( UProbability = self.yUnits ) and (  UPercentProbability = newUnits ) ) then // Convert proportion to percent: multiply by 100
        multiplier := 100.0
      else
        begin
          result := false;
          multiplier := -1.0;
        end
      ;

      if( result ) then
        begin
          // Convert percent to proportion: divide by 100
          for i := low( _points ) to high( _points ) do
            _points[i].y := _points[i].y * multiplier
          ;

          self.yUnits := newUnits;
        end
      ;
    end
  ;


  /// Sends information about the Rel to the debug output window
  procedure TRelFunction.debug();
    var
      i: integer;
    begin
      dbcout( 'Relational function: ' + name + '(id ' + intToStr( id ) + ')', true );
      dbcout( 'Number of points: ' + intToStr( length( _points ) ), true );

      for i := low( _points ) to high( _points ) do
        dbcout( '(' + usFloatToStr(_points[i].x) + ', ' + usFloatToStr(_points[i].y) + ')', true )
      ;
      
      dbcout( 'x units of ' + xUnitsString() + ', y units of ' + yUnitsString(), true );
      dbcout( endl, true );
    end
  ;


  {*
    Creates an Rel object from an XML file that was read into xml
    @param xml string containing the contents of a well-formed XML formatted file
    @return Rel chart object
  }
  function createRelFromXml( xml: string ): TRelFunction;
    var
      sdew: TSdew;
      element: Pointer;
    begin
      xml := stripXmlHeaders( xml );
      sdew := TSdew.createFromString( xml );
      element := sdew.GetRootTree();

      try
        result := createRelFromXml( element, sdew, nil ) as TRelFunction;
      except
        result := nil;
      end;
      
      sdew.Free();
    end
  ;


  {*
    Creates an Rel object from an XML file that has already been pre-processed by sdew
    @param element pointer to the DOM XML data structure returned from sdew.GetRootTree()
    @param instance of the Dephi interface into the scew/expat libraries for XML parsing
    @return Rel chart object
  }
  function createRelFromXml( element: pointer; sdew: TSdew ): TRelFunction;
    begin
      result := createRelFromXml( element, sdew, nil ) as TRelFunction;
    end
  ;



  {*
    Helper function used by createRelFromXml() to populate the Rel coordinates from a deprecated format
    @param rel  the relational function object being created
    @param sdew SDEW XML parsing object already holding the xml file
    @comment an rel function element in the former format looks like:
    <movement-control>
      <value>0.5415</value>        <value>1</value>
      <value>1056.3</value>        <value>1</value>
      <units><xdf:unit>day</xdf:unit></units>
      <units><xdf:unitless /></units>
    </movement-control>
  }
  procedure fillValuesFromOldXml( rel: TRelFunction; element: pointer; sdew: TSdew );
    var
      x, y, temp: double;
      i: integer;
      count: integer;
      gotX: boolean;
    begin
      Count := Sdew.GetElementCount( Element );
      gotX := false;
      x := 0.0;
    
      for i := 0 to (Count-1) do
        begin
          if ( 'value' = Sdew.GetElementName(Sdew.GetElementByIndex( Element, i) ) ) then
            begin
              temp := usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByIndex( Element, i)) );
              if ( gotX ) then
                begin
                  y := temp;
                  gotX := false;
                  rel.addPoint( x, y );
                end
              else
                begin
                  x := temp;
                  gotX := true;
                end
              ;
            end
          ;
        end
      ;
    end
  ;


  {*
    Helper function used by createRelFromXml() to populate the Rel coordinates from the current format
    @param rel  the relational function object being created
    @param sdew SDEW XML parsing object already holding the xml file
    @comment an rel function element in the current format looks like:
    <movement-control>
      <relational-function name="Unrestricted movement">
        <value><x>0.5415</x>        <y>1</y></value>
        <value><x>1056.3</x>        <y>1</y></value>
        <x-units><xdf:unit>day</xdf:unit></x-units>
        <y-units><xdf:unit>proportion</xdf:unit></y-units>
      </relational-function>
    </movement-control>
  }
  procedure fillValuesFromNewXml( rel: TRelFunction; element: pointer; sdew: TSdew );
    var
      x, y: double;
      i: integer;
      e, xe, ye: pointer;
    begin
      for i := 0 to sdew.getElementCount( element ) - 1 do
        begin
          e := sdew.GetElementByIndex( element, i );
          if( 'value' = sdew.GetElementName( e ) ) then
            begin
              xe := sdew.GetElementByName( e, 'x' );
              ye := sdew.GetElementByName( e, 'y' );

              if( nil <> xe ) then
                x := usStrToFloat( sdew.GetElementContents( xe ), NaN )
              else
                x := NaN
              ;
              if( nil <> ye ) then
                y := usStrToFloat( sdew.GetElementContents( ye ), NaN )
              else
                y := NaN
              ;

              rel.addPoint( x, y );
            end
          ;
        end
      ;
    end
  ;


  {*
    Creates an Rel object from an XML file that has already been pre-processed by sdew
    @param element pointer to the DOM XML data structure returned from sdew.GetRootTree()
    @param instance of the Dephi interface into the scew/expat libraries for XML parsing
    @param extData Currently not used
    @return Rel chart object
    @comment This function is the workhorse used by the other XML parsing functions to
    actually parse the XML file contents and create and populate the Rel function object
  }
  function createRelFromXml( Element: pointer; Sdew:TSdew; extData: Pointer ): TObject;
    var
      xUnits: TChartUnitType;
      yUnits: TChartUnitType;
      unitElement: Pointer;
      nextUnitElement: Pointer;

      oldStyle, newStyle: boolean;

      ret_val: TRelFunction;

      tmp: pointer;
    begin
      // Eventually, one of these must be true...
      oldStyle := false;
      newStyle := false;

      ret_val := TRelFunction.create();

      if( 'relational-function' = sdew.GetElementName( element ) ) then
        begin
          // This is new-style XML.
          newStyle := true;
        end
      else
        begin
          tmp := sdew.GetElementByName( element, 'relational-function' );
          if( nil <> tmp ) then
            begin
              element := tmp;
              newStyle := true;
            end
          else
            oldStyle := true
          ;
        end
      ;

      // Determining the function name
      //------------------------------
      // Old-style XML had an optional "name" tag.
      if( nil <> Sdew.GetElementByName( Element, 'name' ) ) then
        begin
          ret_val.name := Sdew.getElementContents( Sdew.GetElementByName( Element, 'name' ) );
          oldStyle := true;
        end

      // New-style XML has an optional "name" attribute.
      else if( nil <> sdew.getAttributeByName( Element, 'name' ) ) then
        begin
          ret_val.name := sdew.GetElementAttribute( element, 'name' );
          newStyle := true;
        end
      ;

      if( oldStyle and newStyle ) then
        raise exception.create( 'Multiple XML schema versions seem to be mixed in createRelFromXml()' )
      ;


      // Determining the function units
      //-------------------------------
      xUnits := UUnknown;
      yUnits := UUnknown;

      // FIX ME: The "per day" units aren't handled correctly here.

      // Old-style XML had two tags, both called "units", the first of which was for the x units and the second for y units.
      // If a "units" tag exists, this is an old-style XML block.
      unitElement := Sdew.GetElementByName( Element, 'units' );
      if ( nil <> unitElement ) then
        begin
          oldStyle := true;

          if ( nil <> Sdew.GetElementByName( unitElement, 'xdf:unit') ) then
            xUnits := chartUnitTypeFromXml( Sdew.GetElementContents( Sdew.GetElementByName( unitElement, 'xdf:unit' ) ) )
          else if ( nil <> Sdew.GetElementByName( unitElement, 'xdf:unitless') ) then
            xUnits := UUnitless
          else
            xUnits := UUnknown
          ;

          nextUnitElement := Sdew.GetNextElement( Element, unitElement );
          if ( nil <> nextUnitElement ) then
            begin
              if ( nil <> Sdew.GetElementByName( nextUnitElement, 'xdf:unit') ) then
                yUnits := chartUnitTypeFromXml( Sdew.GetElementContents( Sdew.GetElementByName( nextUnitElement, 'xdf:unit' ) ) )
              else if ( nil <> Sdew.GetElementByName( nextUnitElement, 'xdf:unitless') ) then
                yUnits := UUnitless
              else
                yUnits := UUnknown
              ;
            end
          else
            yUnits := UUnknown
          ;
        end
      // If there is no "units" tag, look for the new-style "x-units" and "y-units" tags.
      else
        begin
          unitElement := Sdew.GetElementByName( Element, 'x-units' );
          if ( nil <> unitElement ) then
            begin
              newStyle := true;

              if ( nil <> Sdew.GetElementByName( unitElement, 'xdf:unit') ) then
                xUnits := chartUnitTypeFromXml( Sdew.GetElementContents( Sdew.GetElementByName( unitElement, 'xdf:unit' ) ) )
              else if ( nil <> Sdew.GetElementByName( unitElement, 'xdf:unitless') ) then
                xUnits := UUnitless
              else
                xUnits := UUnknown
              ;
            end
          ;

          unitElement := Sdew.GetElementByName( Element, 'y-units' );
          if ( nil <> unitElement ) then
            begin
              newStyle := true;

              if ( nil <> Sdew.GetElementByName( unitElement, 'xdf:unit') ) then
                yUnits := chartUnitTypeFromXml( Sdew.GetElementContents( Sdew.GetElementByName( unitElement, 'xdf:unit' ) ) )
              else if ( nil <> Sdew.GetElementByName( unitElement, 'xdf:unitless') ) then
                yUnits := UUnitless
              else
                yUnits := UUnknown
              ;
            end
          ;
        end
      ;

      // Remember that the NAADSM core library likes proportions, but the interface likes percentages.
      if( UProportion = yUnits ) then
        ret_val.convertYUnitsTo( UPercent )
      else if( UProbability = yUnits ) then
        ret_val.convertYUnitsTo( UPercentProbability )
      ;

      ret_val.xUnits := xUnits;
      ret_val.yUnits := yUnits;

      if( oldStyle and newStyle ) then
        raise exception.create( 'Multiple XML schema versions seem to be mixed in createRelFromXml()' )
      ;


      // Determining the function values
      //--------------------------------
      if( oldStyle ) then
        fillValuesFromOldXml( ret_val, element, sdew )
      else if( newStyle ) then
        fillValuesFromNewXml( ret_val, element, sdew )
      else
        raise exception.create( 'XML schema version cannot be determined in createRelFromXml()' )
      ;


      // Did we actually create anything?  If not, return a nil pointer.
      //----------------------------------------------------------------
      if( 0 = ret_val.pointCount ) then
        begin
          ret_val.Free();
          ret_val := nil;
        end
      ;

      result := ret_val;
    end
  ;


end.
 
{*
ChartFunction.pas
------------------
Begin: 2005/03/15
Last revision: $Date: 2011-01-17 06:48:44 $ $Author: areeves $
Version number: $Revision: 1.18 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2005 - 2011 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This unit deals with the classes associated with input parameters that can be charted as functions by
a series of one or more x and y coordinates that describe a relationship between two variables or a distribution.
}

 (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator
  *)


unit ChartFunction;

{$INCLUDE Defs.inc}

interface

  uses
    StrUtils,
    Sdew
    {$IFDEF DATABASE_ENABLED}
    , SqlClasses
    {$ENDIF}
  ;

  
  type
  /// Enumerations delimiting the allowable set of units for the axes of a chart
  TChartUnitType = (
  	UUnknown,
    UDays,
    UKilometers,
    UPerDay,
    UDegrees,
    UHerdsPerDay,
    UUSDollars,
    UPercent,
    UPercentProbability,
    UUnitless,
    UIndividuals,
    UProportion,
    UProbability,
    UEventCount,
    UCustom
  );

  type
  /// Enumerations delimiting chart type - unspecified, probability density function, or relational
  TChartType = ( CTUnspecified, CTPdf, CTRel );

	function chartTypeAsString( val: TChartType ): string;

  function chartUnitTypeAsString( val: TChartUnitType ): string;
  function chartUnitTypeAsSqlString( val: TChartUnitType ): string;

  function chartUnitTypeFromString( val: string ): TChartUnitType;
  // FIX ME: consider deprecating these and replacing with chartUnitTypeFromString()
  function chartUnitTypeFromXml( val: string ): TChartUnitType;
  function chartUnitTypeFromDatabase( val: string ): TChartUnitType;

  // Currently NAADSM-specific, but it doesn't hurt
  // anything to have it in this unit.
  function chartUnitTypeAsSSXml( val: TChartUnitType ): string;

  
  type
  /// Base class for PDFs and relational functions
  TChartFunction = class
  	protected
    	_xUnits: TChartUnitType; /// x-axis units
      _yUnits: TChartUnitType; /// y-axis units
      _name: string;           /// user-defined chart name
      _id: integer;            /// chart ID stored by the database
      _dbField: word;          /// database field
      _chartType: TChartType;  /// chart type - pdf or relational
      _notes: string;          /// comments about the origin and use of the function

      _indent: string;         /// the number of spaces to indent a line of XML

      procedure setIndent( const indentLevel: integer );

      procedure initialize();

      function getChartType(): TChartType;
      procedure setChartType( val: TChartType );

      function getXUnits() : TChartUnitType;
      procedure setXUnits( val: TChartUnitType );
      function getYUnits() : TChartUnitType;
      procedure setYUnits( val: TChartUnitType );

      function getName(): string;
      procedure setName( Val: string );

      function xUnitsString(): string;
      function yUnitsString(): string;
      
      function getValid(): boolean;

      function getID(): integer;
      procedure setID( val: integer );

      function getNotes(): string;
      procedure setNotes( val: string );

      /// Implemented by decendant classes
      function getDescr(): string; virtual; abstract;

      /// Implemented by decendant classes
      function getYStartsAtZero(): boolean; virtual; abstract;

      procedure setDbField( val: integer );
			function getDbField(): integer;

    public
    	constructor create(); overload; virtual;
      constructor create( const src: TChartFunction ); overload; virtual;

      destructor destroy(); override;

      /// Implemented by decendant classes
      function createCopy(): TChartFunction; virtual; abstract;

      /// Implemented by decendant classes
      function validate( err: pstring = nil ): boolean; virtual; abstract;
      function compare( Chart: TChartFunction ): boolean; virtual;

      /// Implemented by decendant classes
      procedure debug(); virtual; abstract;

      // Return value of -1 indicates failure
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; virtual;
      {$ENDIF}

      function ssXml( const indent: integer ): string; virtual;
      function importXml( element: Pointer; sdew: TSdew ): boolean; virtual;

      /// Implemented by decendant classes - to generate an APPROXIMATE plain-text representation of the chart.
      function asCsv(): string; virtual; abstract;

      property chartType: TChartType read getChartType write setChartType;  /// chart type as TChartType
      property xUnits: TChartUnitType read getXUnits write setXUnits;       /// x-axis units as TChartUnitType
      property yUnits: TChartUnitType read getYUnits write setYUnits;       /// y-axis units as TChartUnitType
      property name: string read getName write setName;                     /// chart name
			property dbField: integer read getDbField write setDbField;           /// database field associated with chart
      property id: integer read getID write setID;                          /// chart ID
			property isValid: boolean read getValid;                              /// chart is composed properly
      property yStartsAtZero: boolean read getYStartsAtZero;                /// whether the first y-axis value is 0

      /// A one-line description of the function, e.g., "Triangular( 3, 7, 15 )" or "Relational function"
      property descr: string read getDescr;

      /// Intended for detailed notes regarding the origin and use of this function
      property notes: string read getNotes write setNotes;
    end
  ;


  
  type
  {*
    Special function pointer type for when functions are updated. This type is used 
    to reference a method of an instance object, so known as a method pointer.
    @param sender instance of form containing updated chart function
    @param fn instance of chart function that was updated
    @comment Delphi topic group procedure types, even though this is a method pointer type
  }
  TUpdateParamsEvent = procedure( sender: TObject; fn: TChartFunction ) of object;


implementation

	uses
  	TypInfo,
    SysUtils,
    I88n,
    
    MyStrUtils,
    DebugWindow,

    ProbDensityFunctions
  ;


  const
  	DBCHARTFUNCTION: boolean = false; /// set to true to enable debugging messages for this unit.


  {*
    Returns the name of the chart type enumeration for val
    @param val the value of a variable of TChartType
    @return the name of the chart type enumeration
  }
  function chartTypeAsString( val: TChartType ): string;
  	begin
			result := getEnumName( TypeInfo( TChartType ), ord( val ) );
    end
  ;


  {*
     Provides a descriptive, untranslated unit type name for val
     @param val enumerated value of a chart unit type
     @return chart unit name or 'Axis units unspecified' if val is unrecognized
     @comment currently no difference from chartUnitTypeAsString
  }
  function chartUnitTypeAsSqlString( val: TChartUnitType ): string;
    begin
      // DO NOT translate the values shown here!!
      case val of
 				UDays: result := 'Days';
        UKilometers: result := 'Kilometers';
        UUnitless: result := '(Unitless)';
        UPerDay: result := 'Per day';
        UDegrees: result := 'Degrees';
        UHerdsPerDay: result := 'Units per day';
        UUSDollars: result := 'US dollars';
        UPercent: result := 'Percent';
        UPercentProbability: result := 'Percent probability';
        UIndividuals: result := 'Individuals';
        UProportion: result := 'Proportion';
        UProbability: result := 'Probability';
        UEventCount: result := 'Number of events';
        UCustom: result := 'Custom units';
        else
        	begin
         		//raise exception.Create( 'Unrecognized TChartUnitType in chartUnitTypeAsString' );
            result := 'Axis units unspecified';
          end
        ;
      end;
    end
  ;


  {*
     Provides a descriptive, untranslated unit type name for val
     @param val enumerated value of a chart unit type
     @return chart unit name or 'Axis units unspecified' if val is unrecognized
  }
  function chartUnitTypeAsString( val: TChartUnitType ): string;
  	begin
      case val of
 				UDays: result := tr( 'Days' );
        UKilometers: result := tr( 'Kilometers' );
        UUnitless: result := tr( '(Unitless)' );
        UPerDay: result := tr( 'Per day' );
        UDegrees: result := tr( 'Degrees' );
        UHerdsPerDay: result := tr( 'Units per day' );
        UUSDollars: result := tr( 'US dollars' );
        UPercent: result := tr( 'Percent' );
        UPercentProbability: result := tr( 'Percent probability' );
        UIndividuals: result := tr( 'Individuals' );
        UProportion: result := tr( 'Proportion' );
        UProbability: result := tr( 'Probability' );
        UEventCount: result := tr( 'Number of events' );
        UCustom: result := tr( 'Custom units' );
        else
        	begin
         		//raise exception.Create( 'Unrecognized TChartUnitType in chartUnitTypeAsString' );
            result := tr( 'Axis units unspecified' );
          end
        ;
      end;
    end
  ;

  {*
     Returns the chart unit type for the value of an XML unit element
     @param val one of the accepted units associated with a chart unit type
     @return the chart unit type for val
     @comment 'km', parsed from <units><xdf:unit>km</xdf:unit></units>, returns UKilometers
  }
  function chartUnitTypeFromXml( val: string ): TChartUnitType;
    begin
      result := chartUnitTypeFromString( val );
    end
  ;


  {*
     Returns the chart unit type for a unit name queried from the database
     @param val one of the accepted units associated with a chart unit type
     @return the chart unit type for val
     @comment the val 'km', queried from inChart.XAxisUnits returns UKilometers
  }
  function chartUnitTypeFromDatabase( val: string ): TChartUnitType;
    begin
      result := chartUnitTypeFromString( val );
    end
  ;


  {*
    Returns the chart unit type for val
    @param val one of the accepted units associated with a chart unit type
    @return the chart unit type for val
    @comment the val 'km' or 'kilometers' returns UKilometers
  }
  function chartUnitTypeFromString( val: string ): TChartUnitType;
  	begin
    	val := fixup( val );

      if( ( val = 'days' ) or ( val = 'day' ) ) then
      	result := UDays
      else if( ( val = 'kilometers' ) or ( val = 'km' ) ) then
      	result := UKilometers
      else if( ( val = '(unitless)' ) or ( val = 'unitless' ) ) then
      	result := UUnitless
      else if( val = 'per day' ) then
      	result := UPerDay
      else if( val = 'degrees' ) then
      	result := UDegrees
      else if( ( val = 'units per day' ) or ( val = 'herds per day' ) ) then
      	result := UHerdsPerDay
      else if( val = 'us dollars' ) then
      	result := UUSDollars
      else if( val = 'percent' ) then
      	result := UPercent
      else if( val = 'percent probability' ) then
        result := UPercentProbability
      else if( val = 'units per day' ) then
      	result := UHerdsPerDay
      else if( val = 'individuals' ) then
        result := UIndividuals
      else if( val = 'proportion' ) then
        result := UProportion
      else if( val = 'probability' ) then
        result := UProbability
      else if( val = 'number of events' ) then
        result := UEventCount
      else if( val = 'custom units' ) then
        result := UCustom
      else
      	result := UUnknown
      ;
    end
  ;


  {*
     Returns the XML formatted line vor val
     @param val a TChartUnitType enumeration
     @return the name for the unit type embedded in XML elements recognizable to the NAADSM parameters XML schema
     @comment UKilometers returns '<units><xdf:unit>km</xdf:unit></units>'
     @comment an exception is raised and 'unknown' returned if val is not recognized
  }
  function chartUnitTypeAsSSXml( val: TChartUnitType ): string;
  	begin
    	case val of
      	UDays: result := '<units><xdf:unit>day</xdf:unit></units>';
        UKilometers: result := '<units><xdf:unit>km</xdf:unit></units>';
        UUnitless: result := '<units><xdf:unitless /></units>';
        UPerDay: result := '<units><xdf:unit power="-1">day</xdf:unit></units>';
        UDegrees: result := '<units><xdf:unit>degree</xdf:unit></units>';
        UHerdsPerDay: result := '<units><xdf:unit>herd</xdf:unit><xdf:unit power="-1">day</xdf:unit></units>';
        UUSDollars: result := '<units><xdf:unit>USD</xdf:unit></units>';
        UPercent: result := '<units><xdf:unit>percent</xdf:unit></units>';
        UPercentProbability: result := '<units><xdf:unit>percent probability</xdf:unit></units>';
        UIndividuals: result := '<units><xdf:unit>individuals</xdf:unit></units>';
        UProportion: result := '<units><xdf:unit>proportion</xdf:unit></units>';
        UProbability: result := '<units><xdf:unit>probability</xdf:unit></units>';
        UEventCount: result := '<units><xdf:unit>event count</xdf:unit></units>';
        UCustom: result := '<units><xdf:unit>custom units</xdf:unit></units>';
      else
      	begin
        	raise exception.Create( 'Unrecognized TChartUnitType (#' + intToStr( integer( val ) ) + ') in chartUnitTypeAsSSXml.' );
      		result := 'unknown';
        end
      end;
    end
  ;

  {*
    Creates a safe (initialized) but empty chart function object
  }
  constructor TChartFunction.create();
    begin
    	inherited create();
      initialize();
    end
  ;


  {*
    Creates a copy of src by copying it's parameters
    @param src source ChartFunction object
  }
  constructor TChartFunction.create( const src: TChartFunction );
    begin
      dbcout( '*** Calling TChartFunction copy constructor', DBCHARTFUNCTION );

      inherited create();

      setChartType( src._chartType );
      setXUnits( src._xUnits );
      setYUnits( src._yUnits );
      setName( src._name );

      setDbField( src._dbField );
      setID( src._id );
      
      setNotes( src._notes );
    end
  ;


  {*
     Sets the value of private members to appropriate initial conditions
  }
  procedure TChartFunction.initialize();
    begin
      setChartType( CTUnspecified );
      setXUnits( UUnknown );
      setYUnits( UUnknown );
      setName( '' );
      setNotes( '' );
    end
  ;


  {*
    Frees the object from memory
  }
  destructor TChartFunction.destroy();
  	begin
    	dbcout( 'Chart function destroyed: ' + name, DBCHARTFUNCTION );
    	inherited destroy();
    end
  ;


  {*
    Compares some key attributes of self and Chart
    @param Chart chart function object to compare with self
    @return true if key attributes have the same values else false
    @comment chart type, chart name, and axes units are compared
  }
  function TChartFunction.compare( Chart:TChartFunction ):boolean;
    var
      myNameForComparison, otherNameForComparison: string;

      // FIX ME: should this go into MyStrUtils?
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
      result := false;

      if( Assigned( Chart ) ) then
        begin
          myNameForComparison := stripParens( self.name );
          otherNameForComparison := stripParens( Chart.name );

          result := (
            ( self.chartType = Chart.chartType )
            //and ( self._dbField = Chart._dbField ) // The dbField won't necessarily be the same: some fields share charts, and this behavior is allowed.
            and ( self.xUnits = Chart.xUnits )
            and ( self.yUnits = Chart.yUnits )
            and ( myNameForComparison = otherNameForComparison )
          );
        end
      ;
    end
  ;

  {*
    Returns the result of the abstract method validate, implemented by the decendant class
    @return true if the chart is valid else false
  }
  function TChartFunction.getValid(): boolean;
  	begin
   		result := validate();
    end
  ;


  {*
     Returns the name/label of the x-axis chart unit type
     @return name of the x-axis chart unit type
  }
  function TChartFunction.xUnitsString(): string;
    begin
    	result := chartUnitTypeAsString( xUnits );
    end
  ;


  {*
     Returns the name/label of the y-axis chart unit type
     @return name of the y-axis chart unit type
  }
  function TChartFunction.yUnitsString(): string;
    begin
    	result := chartUnitTypeAsString( yUnits );
    end
  ;

  /// get function for property chartType returning the value of _chartType
  function TChartFunction.getChartType(): TChartType; begin result := _chartType; end;
  /// set function for property chartType setting the value of _chartType from val
  procedure TChartFunction.setChartType( val: TChartType ); begin _chartType := val; end;

  /// get function for property xUnits returning the value of _xUnits
  function TChartFunction.getXUnits() : TChartUnitType; begin Result := _xUnits; end;
  /// set function for property xUnits setting the value of _xUnits from val
  procedure TChartFunction.setXUnits( val: TChartUnitType ); begin _xUnits := val; end;

  /// get function for property yUnits returning the value of _yUnits
  function TChartFunction.getYUnits() : TChartUnitType; begin Result := _yUnits; end;
  /// set function for property yUnits setting the value of _yUnits from val
  procedure TChartFunction.setYUnits( val: TChartUnitType ); begin _yUnits := val; end;

  /// get function for property name returning the value of _name
  function TChartFunction.getName(): string; begin result := _name; end;
  /// set function for property name setting the value of _name from val
  procedure TChartFunction.setName( val: string ); begin _name := val; end;

  /// get function for property id returning the value of _id
  function TChartFunction.getID(): integer; begin result := _id; end;
  /// set function for property id setting the value of _id from val
  procedure TChartFunction.setID( val: integer ); begin _id := val; end;

  /// get function for property dbField returning the value of _dbField
  function TChartFunction.getDbField(): integer; begin Result := _dbField; end;
  /// set function for property dbField setting the value of _dbField from val
	procedure TChartFunction.setDbField( val: integer ); begin _dbField := val; end;

  /// get function for property notes returning the value of _notes
  function TChartFunction.getNotes(): string; begin result := _notes; end;
  /// set function for property notes setting the value of _notes from val
  procedure TChartFunction.setNotes( val: string ); begin _notes := val; end;

  {*
    set function for _indent setting the value from indentLevel
    @comment method currently not used
  }
  procedure TChartFunction.setIndent( const indentLevel: integer );
    var
      i: integer;
    begin
      _indent := '';

      for i := 0 to indentLevel - 1 do
        begin
          _indent := _indent + '  ';
        end
      ;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  /// Safety measure that does nothing but raise an exception for calling this function from the base class
  function TChartFunction.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
    begin
      raise exception.create( '"Abstract" error: call made to TChartFunction.populateDatabase()' );
      result := 0;
    end
  ;
  {$ENDIF}

  /// Safety measure that does nothing but raise an exception for calling this function from the base class
  function TChartFunction.ssXml( const indent: integer ): string;
    begin
      raise exception.create( '"Abstract" error: call made to TChartFunction.ssXml()' );
      result := '';
    end
  ;

  /// Safety measure that does nothing but raise an exception for calling this function from the base class
  function TChartFunction.importXml( element: Pointer; sdew: TSdew ): boolean;
    begin
      raise exception.create( '"Abstract" error: call made to TChartFunction.importXml()' );
      result := false;
    end
  ;

  
end.

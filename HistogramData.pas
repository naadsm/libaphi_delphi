{*
HistogramData.inc
-----------------
Begin: 2010/03/28
Last revision: $Date: 2011-01-24 02:00:38 $ $Author: areeves $
Version: $Revision: 1.6 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2010 - 2011 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This unit deals with generating a histogram data structure and several
forms of output from a vector containing data. The bin number of the histogram
can be calculated several different ways.
}

  (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator
  *)


unit HistogramData;

interface

  uses
    MyDelphiArrayUtils,

    QVectors,
    ProbDensityFunctions
  ;

  type
  /// Enumeration for methods to calculate the number of histogram bins
  THistBreakType = (
    Unspecified,
    UserSpecified,
    Sturges,          /// bases the bin number on the range of the data 
    Scott,            /// bin width calculation uses the sample standard deviation
    QuarterScott,     /// uses 1/4 as many bins as the Scott algorithm
    SquareRoot,       /// bin number based on the square root of the number of data points
    FreedmanDiaconis  /// bin width calculation based on the interquartile range
  );

  {*
    This class generates an R-style histogram data structure from a vector containing data.
    Several algorithms (listed in THistBreakType) aare supported for determining break points.
    The following websites/references describe these algorithms:
    
    http://en.wikipedia.org/wiki/Histogram
    http://stat.ethz.ch/R-manual/R-patched/library/grDevices/html/nclass.html
    http://www.vosesoftware.com/ModelRiskHelp/index.htm#Presenting_results/Histogram_plots/Sturges_rule.htm
    http://www.vosesoftware.com/ModelRiskHelp/index.htm#Presenting_results/Histogram_plots/Determining_the_width_of_histogram_bars.htm
    
    Versions of these references, current as of the this writing, are saved in the HistogramDemo/Docs folder 
    associated with this file. 
    
    This class does not generate a visual representation of the histogram.  To display the histogram,
    use TFrameArrayHistogram (see FrameArrayHistogram.pas/dfm in libaphi_delphi_gui).
  *}
  type
  /// This class generates an R-style histogram data structure from a vector containing data.
  THistogramData = class( TObject )
    protected
      _data: TQDoubleVector;       /// array of data values
      _breakType: THistBreakType;  /// method used to calculate the number of histogram bins
      _breaks: TQDoubleVector;     /// x-axis bin interval values
      _densities: TQDoubleVector;  /// the proportion of cases in each bin (y-axis category values)
      _counts: TQDoubleVector;     /// the number of cases in each bin (y-axis category values)
      _nBins: integer;             /// number of bins (categories)
      _histGenerated: boolean;     /// indicates whether the histogram was generated

      procedure initialize();

      function getIsEmpty(): boolean;
      function getNBins(): integer;

      function getBreaks(): TQDoubleVector;
      function getCounts(): TQDoubleVector;
      function getDensities(): TQDoubleVector;

      procedure generateHistogram();
      procedure calculateBreaks();
      function dataChanged(): boolean;
      function breaksChanged( const breaks: THistBreakType ): boolean;

      procedure setBreakType1( const breaks: THistBreakType );

    public
      constructor create( const breaks: THistBreakType = Scott; const nBins: integer = 50 ); overload;
      constructor create( const data: TQDoubleVector; const breaks: THistBreakType = Scott; const nBins: integer = 50 ); overload;
      constructor create( const data: TARDoubleArray; const breaks: THistBreakType = Scott; const nBins: integer = 50 ); overload;
      constructor create( const data: TARIntArray; const breaks: THistBreakType = Scott; const nBins: integer = 50 ); overload;
      constructor create( const data: TQIntegerVector; const breaks: THistBreakType = Scott; const nBins: integer = 50 ); overload;

      destructor destroy(); override;

      function createPdf(): TPdfHistogram; overload;
      function createPdf( const breakType: THistBreakType ): TPdfHistogram; overload;
      //function createPdf( breaks: TQDoubleVector ): TPdfHistogram; overload;

      function asNaadsmCsv(): string; overload;
      function asNaadsmCsv( const breakType: THistBreakType ): string; overload;
      //function asNaadsmCsv( breaks: TQDoubleVector ): string; overload;

      procedure clear();

      procedure setBreakType( const breakType: THistBreakType; const nBins: integer );

      property data: TQDoubleVector read _data;             /// read-only access to _data
      property breaks: TQDoubleVector read getBreaks;       /// read-only access to _breaks
      property densities: TQDoubleVector read getDensities; /// read-only access to _densities
      property counts: TQDoubleVector read getCounts;       /// read-only access to _counts
      property isEmpty: boolean read getIsEmpty;            /// true if _data.count is 0, else false
      property breakType: THistBreakType read _breakType write setBreakType1;  /// bin calculation method
      property nBins: integer read getNBins;                /// read-only access to _nBins
    end
  ;

implementation

  uses
    Math,
    SysUtils,

    DebugWindow,
    MyStrUtils,
    RoundToXReplacement_3c,

    ChartFunction
  ;

  {*
     Creates an empty histogram data structure
     @param breaks bin calculation method, if unspecified defaults to Scott
     @param nBins number of bins, if unspecified defaults to 50
  }
  constructor THistogramData.create( const breaks: THistBreakType = Scott; const nBins: integer = 50 );
    begin
      inherited create();
      initialize();
      _nBins := max( nBins, 1 );
      setBreakType1( breaks );
    end
  ;

  {*
     Creates an histogram data structure for floating point data
     @param data input data contained in a TQDoubleVector data structure
     @param breaks bin calculation method, if unspecified defaults to Scott
     @param nBins number of bins, if unspecified defaults to 50
  }
  constructor THistogramData.create( const data: TQDoubleVector; const breaks: THistBreakType = Scott; const nBins: integer = 50 );
    begin
      inherited create();
      initialize();
      
      _data.assign( data );
      _data.sort();

      _nBins := max( nBins, 1 );
      setBreakType1( breaks );
    end
  ;


  {*
     Creates an histogram data structure for floating point data
     @param data input data contained in a TARDoubleArray data structure
     @param breaks bin calculation method, if unspecified defaults to Scott
     @param nBins number of bins, if unspecified defaults to 50
  }
  constructor THistogramData.create( const data: TARDoubleArray; const breaks: THistBreakType = Scott; const nBins: integer = 50 );
    var
      i: integer;
    begin
      inherited create();
      initialize();

      for i := 0 to length( data ) - 1 do
        _data.append( data[i] )
      ;
      _data.sort();

      _nBins := max( nBins, 1 );
      setBreakType1( breaks );
    end
  ;


  {*
     Creates an histogram data structure for integer data
     @param data input data contained in a TARIntArray data structure
     @param breaks bin calculation method, if unspecified defaults to Scott
     @param nBins number of bins, if unspecified defaults to 50
  }
  constructor THistogramData.create( const data: TARIntArray; const breaks: THistBreakType = Scott; const nBins: integer = 50 );
    var
      i: integer;
    begin
      inherited create();
      initialize();

      for i := 0 to length( data ) - 1 do
        _data.append( data[i] )
      ;
      _data.sort();

      _nBins := max( nBins, 1 );
      setBreakType1( breaks );
    end
  ;


  {*
     Creates an histogram data structure for integer data
     @param data input data contained in a TQIntegerVector data structure
     @param breaks bin calculation method, if unspecified defaults to Scott
     @param nBins number of bins, if unspecified defaults to 50
  }
  constructor THistogramData.create( const data: TQIntegerVector; const breaks: THistBreakType = Scott; const nBins: integer = 50 );
    var
      i: integer;
    begin
      inherited create();
      initialize();

      for i := 0 to data.count - 1 do
        _data.append( data[i] )
      ;
      _data.sort();

      _nBins := max( nBins, 1 );
      setBreakType1( breaks );
    end
  ;


  /// Helper method called by the constructors to create and initalize the private members
  procedure THistogramData.initialize();
    begin
      _data := TQDoubleVector.create();
      _breaks := TQDoubleVector.create();
      _densities := TQDoubleVector.create();
      _counts := TQDoubleVector.create();

      _histGenerated := false;

      _breakType := Unspecified;
      _nBins := 50;
    end
  ;


  /// Frees object and private data structures from memory
  destructor THistogramData.destroy();
    begin
      _data.Free();
      _breaks.Free();
      _densities.Free();
      _counts.Free();

      inherited destroy();
    end
  ;


  {*
    Helper method that regenerates the historgram if breaks specifies a different method than _breakType
    @param breaks bin calculation method
  }
  procedure THistogramData.setBreakType1( const breaks: THistBreakType );
    var
      regenerate: boolean;
    begin
      regenerate := ( breaksChanged( breaks ) or dataChanged() or ( not _histGenerated ) );
      _breakType := breaks;

      if( _data.isEmpty ) then
        clear()
      else if( regenerate ) then
        generateHistogram()
      ;
    end
  ;


  {*
    Set method for property breakType, also setting _nBins to nBins and regenerating the histogram
    @param breakType one of the enumerated bin calculation methods
    @param nBins suggested number of bins
  }
  procedure THistogramData.setBreakType( const breakType: THistBreakType; const nBins: integer );
    begin
      // This will force breaks to be reset based on the new number of bins.
      if( ( UserSpecified = breakType ) and ( nBins <> _nBins ) ) then
        _breakType := Unspecified
      ;

      _nBins := nBins;
      setBreakType1( breakType );
    end
  ;


  {*
     Calculates the number of bins and bin width based on the selected method
     and then generates the necessary values for _breaks.
  }
  procedure THistogramData.calculateBreaks();
    var
      iqr: double;
      binWidth: double;
      dataRange: double;
      i: integer;
      dataMin: double;
      dataMax: double;
    begin
      _breaks.clear();

      _data.sort();
      dataMin := _data.at( 0 );
      dataMax := _data.at( _data.count - 1 );

      if( dataMin = dataMax ) then
        begin
          _nBins := 1;
          _breaks.append( dataMax - 0.5 );
          _breaks.append( dataMax + 0.5 );
        end
      else
        begin
          dataRange := dataMax - dataMin;
      
          case _breakType of
            Sturges:
              begin
                _nBins := ceil( log2( _data.count ) + 1 );
              end
            ;
            Scott:
              begin
                binWidth := ( ( 3.5 * _data.stddev() )/( power( _data.count, 0.33333 ) ) );
                binWidth := max( binWidth, 1.0 );
                _nBins := ceil( dataRange / binWidth );
              end
            ;
            QuarterScott:
              begin
                binWidth := ( ( 3.5 * _data.stddev() )/( power( _data.count, 0.33333 ) ) );
                binWidth := max( binWidth, 1.0 );
                _nBins := ceil( ( dataRange / binWidth ) / 4 );
              end
            ;
            SquareRoot:
              begin
                _nBins := ceil( sqrt( _data.count ) );
              end
            ;
            FreedmanDiaconis:
              begin
                iqr := _data.quantile( 0.75 ) - _data.quantile( 0.25 );
                binWidth := ( ( 2 * iqr )/( power( _data.count, 0.33333 ) ) );
                binWidth := max( binWidth, 1.0 );
                _nBins := ceil( dataRange / binWidth );
              end
            ;
            UserSpecified:
              begin
                _nBins := _nBins; // Cool, huh?
              end
            ;
            else
              begin
                raise exception.create( 'Unsupported break type in THistogramData.calculateBreaks()' );
                _nBins := 1;
              end
            ;
          end;

          _nBins := max( _nBins, 1 );

          binWidth := dataRange / _nBins;

          for i := 0 to _nBins - 1 do
            _breaks.append( ( i * binWidth ) + dataMin )
          ;
          _breaks.append( dataMax );
        end
      ;
    end
  ;


  {*
     Calculates the frequencies and densities for each bin
  }
  procedure THistogramData.generateHistogram();
    var
      i: integer;
      val: double;
      bin: integer;
      count: integer;
      brk: double;
    begin
      _breaks.clear();
      _counts.clear();
      _densities.clear();

      if( data.isEmpty ) then
        raise exception.Create( 'There is no data in THistogramData.generateHistogram()' )
      ;

      _data.sort();

      calculateBreaks();

      bin := 1;
      count := 0;
      brk := _breaks.at( bin );
      i := 0;

      // If there's only 1 bin, we don't have to sort anything into individual bins.
      if( 2 = _breaks.count ) then
        _counts.append( _data.count )
      else
        begin
          while( i < _data.count ) do
            begin
              val := _data.at( i );
              //dbcout2( 'brk: ' + usFloatToStr( brk ) + ', val: ' + usFloatToStr( val ) + ', count: ' + intToStr( count ) );
              if( val <= brk ) then
                begin
                  inc( count );
                  inc( i );
                end
              else
                begin
                  _counts.append( count );
                  count := 0;
                  if( bin < _breaks.count - 2 ) then
                    begin
                      inc( bin );
                      brk := _breaks.at( bin );
                    end
                  else
                    break
                  ;
                end
              ;
            end
          ;
          _counts.append( _data.count - i ); // Don't forget the last bin.
        end
      ;

      // Calculate densities
      for i := 0 to _counts.count - 1 do
        _densities.append( _counts.at(i) / _data.count )
      ;

      _histGenerated := true;
    end
  ;


  {*
     Indicates if the input data has changed
     @param true if input has changed else false
  }
  function THistogramData.dataChanged(): boolean;
    begin
      result := not( _data.sorted );
    end
  ;


  {*
    Indicates if breaks is the same bin calculation method as _breakType
    @return true if the values are the same, else false
  }
  function THistogramData.breaksChanged( const breaks: THistBreakType ): boolean;
    begin
      result := ( breaks <> _breakType );
    end
  ;


  {*
     Creates an probaility density function object based on the histogram values calculated with breakType
     @param breakType one of the enumerated bin calculation methods
     @return a unitless PDF object holding the x and y values of the histogram
  }
  function THistogramData.createPdf( const breakType: THistBreakType ): TPdfHistogram;
    begin
      setBreakType1( breakType );
      result := createPdf();
    end
  ;


  {*
     Creates an probability density function object based on the histogram values
     @return a unitless PDF object holding the x and y values of the histogram
  }
  function THistogramData.createPdf(): TPdfHistogram;
    var
      rangesArr, countsArr: TARDoubleArray;
    begin
      if( nil = _data ) then
        result := nil
      else if( _data.isEmpty ) then
        result := nil
      else
        begin
          rangesArr := breaks.createDoubleArray();
          countsArr := counts.createDoubleArray();

          breaks.debug();
          counts.debug();

          result := TPdfHistogram.create( rangesArr, countsArr, UUnitless );

          setLength( rangesArr, 0 );
          setLength( countsArr, 0 );
        end
      ;
    end
  ;

  /// Clears all the private member data structures
  procedure THistogramData.clear();
    begin
      _data.clear();
      _breaks.clear();
      _densities.clear();
      _counts.clear();
      _histGenerated := false;
    end
  ;

  {*
    Uses breakType to generate the histogram formatted as a list of comma-separated x,y values
    @param breakType one of the enumerated bin calculation methods
    @return histogram values in NAADSM CSV format
  }
  function THistogramData.asNaadsmCsv( const breakType: THistBreakType ): string;
    begin
      setBreakType1( breakType );
      result := asNaadsmCsv();
    end
  ;


  {*
    Generates the histogram formatted as a list of comma-separated x,y values
    @return histogram values in NAADSM CSV format
  }
  function THistogramData.asNaadsmCsv(): string;
    var
      i: integer;
    begin
      //breaks.debug();
      //counts.debug();
      //densities.debug();

      result := 'x, y' + endl;
      for i := 0 to _counts.count - 1 do
        result := result + uiFloatToStr( breaks.at(i) ) + ', ' + intToStr( roundDbl( _counts.at(i) ) ) + endl
      ;
      result := result + uiFloatToStr( breaks.at( breaks.count - 1 ) ) + ', 0' + endl
    end
  ;

  {*
    Indicates whether the input data structure is empty
    @return true if _data contains no data else false
  }
  function THistogramData.getIsEmpty(): boolean;
    begin
      result := ( 0 = _data.count );
    end
  ;


  {*
    Get method for property nBins
    @return number of bins in the histogram
  }
  function THistogramData.getNBins(): integer;
    begin
      result := _nBins;
    end
  ;


  {*
    Get method for property breaks
    @return bin locations along the x-axis of the histogram
  }
  function THistogramData.getBreaks(): TQDoubleVector;
    begin
      if( _breaks.isEmpty ) then
        generateHistogram()
      ;
      result := _breaks;
    end
  ;


  {*
    Get method for property counts
    @return bin frequencies, which can comprise the y-axis of the histogram
  }
  function THistogramData.getCounts(): TQDoubleVector;
    begin
      if( _counts.isEmpty ) then
         generateHistogram()
      ;
      result := _counts;
    end
  ;


  {*
    Get method for property densities
    @return bin proportions, which can comprise the y-axis of the histogram
  }
  function THistogramData.getDensities(): TQDoubleVector;
    begin
      if( _densities.isEmpty ) then
        generateHistogram()
      ;
      result := _densities;
    end
  ;


end.

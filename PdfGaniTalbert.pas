{*
PdfGaniTalbert.pas
------------------
Begin: 2009/01/29
Last revision: $Date: 2010-10-28 21:17:00 $ $Author: rhupalo $
Version: $Revision: 1.17 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2009 - 2010 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
--------------------------------------------------
This distribution, described by Gani and revisited by Talbert based on work by Stirling, can
be used to address questions like:
  Suppose that we have a set of n items.  We will select k items from this set, replacing each
  selected item, so that it is returned to the set before the next selection is made.
  Each subsequent selection is made at random, without considering whether an
  item had been previously selected.  In the end, how many unique items (call this number j)
  will we have selected?  In other words, how many items will have been selected at least once?

}
  (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator

    Base classes have been documented pretty well, the derived classes not so much
  *)

unit PdfGaniTalbert;

interface

  uses
    QVectors,
    TwoDIntArray,

    ChartFunction,
    ProbDensityFunctions
  ;


  type
  /// Enumrations for the algorithms that cab be used to compute the distribution
  TGaniTalbertAlgorithm = (
    TDBest,
    TDStirling32,
    TDExpMatrix, // Should not be used, except for demonstration: see comments below.
    TDNormalApprox,
    TDStirlingGmp
  );

  {*
    Suppose that we have a set of n items.  We will select k items from this set, replacing each
    selected item, so that it is returned to the set before the next selection is made.
    Each subsequent selection is made at random, without considering whether an
    item had been previously selected.  In the end, how many unique items (call this number j)
    will we have selected?  In other words, how many items will have been selected at least once?

    This distribution, described by Gani and revisited by Talbert based on work by Stirling, can
    be used to address questions like the one above.

    Several algorithms for this calculation are demonstrated:
      - Use of an exponent matrix, as shown in FIX ME: reference.  This algorithm gives an
        "exact" result, at least within computational limits.  This method is tractable, but slow,
        for values of k up to about 30.
      - Use of Stirling numbers, as shown in FIX ME: reference.  This algorithm gives the same
        result as the use of an exponent matrix, but the computation is far less intensive.  Like
        the exponent matrix algorithm, this approach is tractable for values up to about 30.
      - Use of a Normal approximation.  For k greater than 30, the exact methods cannot be used with 32-bit math:
        The precision of double-precision floating point math is inadequate.  The distribution can be
        approximated, sort of, with a Normal distribution having a mean and variance as calculated.  The
        Normal approximation is not the best available approximation (see FIX ME: reference), but it is
        computationally efficient and easy to understand.
      - Use of multiple precision arithmetic to calculate Stirling numbers.  This approach should
        work for values far greater than 30 and should produce an exact (or nearly exact) answer.

    By default, for k < 30, the Stirling number method is used to calculate a result.  For other values of
    k, the Normal approximation is used.
    FIX ME: Figure out when it's practical to use multiple precision arithmetic.

    For additional information, please see the following sources:
    FIX ME: Include key references here.
    Fix Me: The DelphiDoc algorithm method descriptions need review, especially those related to the exponent matrix algorithm.
            Without the paper it is hard to figure out and accurately describe what some of the methods do.
  }
  type TPdfGaniTalbert = class( TPdfDiscrete )
    protected
      _k: integer;  /// the number of items to be selected from the set
      _n: integer;  /// total number of items in the set

      _mean: double;     /// the mean number of unique items selected
      _variance: double; /// variance
      _stdev: double;    /// standard deviation around the mean

      function getMin(): double; override;
      function getMax(): double; override;

      procedure calcMeanAndVariance();

      // Functions for the Stirling number method
      //------------------------------------------
      function probStirling( const j: integer ): double;
      function probStirlingGmp( const j: integer ): double;
      function stirlingNumber( const i, j: integer ): double;

      // Functions for the exponent matrix method
      // This code is included for demonstration, and should
      // not be used. The result is identical to the Stirling
      //  method, but the calculation is much, much slower.
      //-----------------------------------------------------
      function probExpMatrix( const j: integer ): double;
      function calcSum( expMatrix: TTwoDIntArray; const j: integer ): double;
      function createExpMatrix( const j: integer ): TTwoDIntArray;
      procedure fillMatrixRecurse(
        arr: TTwoDIntArray;   // The matrix to which completed arrays will be appended as columns
        v: TQIntegerVector;   // The array to fill
        arrPos: integer;      // Array position
        val: integer;         // The initial value
        const target: integer // The desired sum
      );

      // Functions for the Normal approximation
      //---------------------------------------
      function probNormalApprox( const j: integer ): double;

      // Reimplemented functions
      //------------------------
      procedure initialize(); override;
      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      // Properties
      //-----------
      function getVariance(): double;
      function getStDev(): double;

    public
      // Creation/initialization/destruction
      //------------------------------------
      constructor create(); overload; override;
      constructor create( const n, k: integer; u: TChartUnitType = UUnknown ); overload;
      constructor create( const src: TPdfGaniTalbert ); overload;
      destructor destroy(); override;

      procedure setParams( const n, k: integer; u: TChartUnitType = UUnknown );

      // Housekeeping functions
      //-----------------------
      function validate( err: pstring = nil ): boolean; override;
      function compare( chart: TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      // FIX ME: If this is protected (as it should be), I sometimes get 'abstract' errors.
      function getHistValues(): TMassHistogramValues; override;

      procedure debug(); override;

      // Workhorse functions
      //--------------------
      { What is the probability of having exactly j unique selections, as calculated by the indicated algorithm? }
      function probDensity( const j: integer; const calcType: TGaniTalbertAlgorithm = TDBest ): double; reintroduce;
      function invCumulative( prob: double ): integer; override;
      function rand(): double; override;

      /// read-only access to _k
      property k: integer read _k;
      /// read-only access to _n
      property n: integer read _n;
      /// read-only ccess to _variance
      property variance: double read getVariance;
      /// read-only access to _stDev
      property stDev: double read getStDev;
    end
  ;


  // Library management
  //-------------------
  function libAphiGTLoaded(): boolean;

  /// Empty string if the required library was loaded and can be used.  Otherwise, contains an error message.
  function libAphiGTLoadErrors(): string;


implementation

  uses
    Windows, // Defines THandle
    Forms, // For Application object
    Math,
    SysUtils,

    FunctionPointers,
    ARMathAdvanced,
    RoundToXReplacement_3c,
    DebugWindow,
    MyStrUtils,
    I88n,

    AphiRng
  ;

  const
    REQUIRED_APHI_DLL_VERSION: string = '0.5';  /// Required version of the APHI library
    DBSHOWMSG = false; /// Set to true to enable debugging messages for this unit

  var
    // For the DLL
    /// indicates whether the APHI library is loaded
    _libAphiGTLoaded: boolean;
    /// empty if the library loaded sucessfully, otherwise it contains the load error messages
    _libAphiGTLoadErrors: string;
    /// function pointer to pass n, k, j parameters to a multiple precision Stirling algorithm in the dll; a probability is returned
    aphi_gani_talbert_probability: TCFnDouble_3_Int_Int_Int;


//-----------------------------------------------------------------------------
// Library management
//-----------------------------------------------------------------------------
  {*
    Reads _libAphiGTLoaded
    @return true if the required libraries were loaded and can be used, else false
  }
  function libAphiGTLoaded(): boolean;
    begin
      result := _libAphiGTLoaded;
    end
  ;

  {*
    Reads _libAphiGTLoadErrors
    @return empty if APHI.dll loaded, otherwise contains an erro message
  }
  function libAphiGTLoadErrors(): string;
    begin
      result := _libAphiGTLoadErrors;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Gani-Talbert distribution
//-----------------------------------------------------------------------------

  {*
    Creates the PDF object and initializes with safe but unuseful values
  }
  constructor TPdfGaniTalbert.create();
    begin
      inherited create();
      initialize();
      _n := -1;
      _k := -1;
    end
  ;


  {*
    Creates the PDF object and lets you specify some parameters
    @param n the total number of items in the set
    @param k the number of items to be selected from the set
    @param u enumeration for the type of units for the x axis of the chart
  }
  constructor TPdfGaniTalbert.create( const n, k: integer; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      initialize();
      _n := n;
      _k := k;

      calcMeanAndVariance();
      
      xUnits := u;
    end
  ;


  {*
    Creates the PDF object and copies parameters from src
    @param src source PdfGaniTalbert object
  }
  constructor TPdfGaniTalbert.create( const src: TPdfGaniTalbert );
    begin
      inherited create();
      initialize();
      _n := src.n;
      _k := src.k;

      _mean := src._mean;
      _variance := src._variance;
      _stdev := src._stdev;

      xUnits := src.xUnits;
    end
  ;


  /// Used by the constructors to initialize their private members
  procedure TPdfGaniTalbert.initialize();
    begin
      setPdfType( PdfSpGaniTalbert );

      _infLeftTail := false;
      _infRightTail := false;

      freeHistogramValues();

      _mean := -1.0;
      _variance := -1.0;
      _stdev := -1.0;
    end
  ;


  /// Destroys the object and frees memory
  destructor TPdfGaniTalbert.destroy();
    begin
      inherited destroy();
    end
  ;


  {*
    Resets the parameters and recalculates the chart and distribution statistics
    @param n the total number of items in the set
    @param k the number of items to be selected from the set
    @param u enumeration for the type of units for the x axis of the chart
  }
  procedure TPdfGaniTalbert.setParams( const n, k: integer; u: TChartUnitType = UUnknown );
    begin
      _n := n;
      _k := k;

      xUnits := u;
      freeHistogramValues();

      calcMeanAndVariance();
    end
  ;


  {*
    Creates and returns a copy of itself
    @return ChartFunction object, the base class for Rels and PDFs
    @comment This function and similar ones in the TPDF and TRelFuntion classes
    allows the caller to create an Rel or a PDF as necessary, since the ancestor class is returned.
  }
  function TPdfGaniTalbert.createCopy(): TChartFunction;
    begin
      result := TPdfGaniTalbert.create( self );
    end
  ;


  /// Send some parameter values to the debug window
  procedure TPdfGaniTalbert.debug();
    var
      str: string;
    begin
      str := name + endl;
      str := str
        + 'Gani-Talbert PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 'n=' + usFloatToStr( n )
        + ', k=' + usFloatToStr( k )
      ;

      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  /// Returns 'Gani-Talbert' and the values for n and k
  function TPdfGaniTalbert.getDescr(): string;
    begin
      result := tr( 'Gani-Talbert' ) + format( ' ( %d, %d )', [n, k] );
    end
  ;


  /// Get function for property mean, returning the value of _mean
  function TPdfGaniTalbert.getMean(): double;
    begin
      result := _mean;
    end
  ;

  /// Get function for property variance, returning the value of _variance
  function TPdfGaniTalbert.getVariance(): double;
    begin
      result := _variance;
    end
  ;

  /// Get function for property standard deviation, returning the value of _stdev
  function TPdfGaniTalbert.getStDev(): double;
    begin
      result := _stdev;
    end
  ;

  /// Returns true, indicating that this distribution has a mean
  function TPdfGaniTalbert.getHasMean(): boolean;
    begin
      result := true;
    end
  ;

  {*
    Returns the minimum number of unique items
    @comment Being a discreet distribution, the minimum will typically be 1,
    unless k is 0 then the value returned is 0.
  }
  function TPdfGaniTalbert.getMin(): double;
    begin
      if( 0 = k ) then
        result := 0
      else
        result := 1
      ;
    end
  ;

  /// Returns the value of n or k, whichever is less
  function TPdfGaniTalbert.getMax(): double;
    begin
      result := Math.min( n, k );
    end
  ;


  {*
    Input parameter error checking of n, k, and x-axis units
    @param err error message if parameters do no validate, else it will be empty
    @return true if the parameter values are valid, else false
  }
  function TPdfGaniTalbert.validate( err: pstring = nil ): boolean;
    begin
      if( n < 1 ) then
        begin
          if( nil <> err ) then err^ := tr( 'n must be a greater than 0.' );
          result := false;
        end
      else if( k < 0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'k cannot be negative.' );
          result := false;
        end
       else if( ( UProportion = xUnits ) or ( UProbability = xUnits ) ) then
         begin
           if( nil <> err ) then err^ := tr( 'This probability density function cannot be used to represent a proportion or a probabilitiy.' );
           result := false;
        end
      else
        result := true
      ;
    end
  ;


  {*
    Compares self to chart
    @param chart chart function object to compare with self
    @return true if PDF type, and the values of n and k are the same, else false
  }
  function TPdfGaniTalbert.compare( chart: TChartFunction ):boolean;
    begin
      result := false;

      if( inherited compare( chart ) ) then
        begin
          if
            ( self._n = ( chart as TPdfGaniTalbert ).n )
          and
            ( self._k = ( chart as TPdfGaniTalbert ).k )
          then
            result := true
          ;
        end
      ;
    end
  ;


  {*
     Returns an XML element describing self
     @param indent the number of leading spaces to indent each line
     @return XML formatted text of chart properties
  }
  function TPdfGaniTalbert.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      result := result + _indent + '<gani-talbert>' + endl;
      result := result + _indent + '  <n>' + intToStr( n ) + '</n>' + endl;
      result := result + _indent + '  <k>' + intToStr( k ) + '</k>' + endl;
      result := result + _indent + '</gani-talbert>' + endl;
      result := result + _indent + chartUnitTypeAsSSXml( xUnits ) + endl;
    end
  ;


  {$IFDEF DATABASE_ENABLED}

  {*
    Copies the data from self to the database
    @param db object representing the database
    @param update if true performs an update query, else an insert query
    @return the key value for this record in table inChartDetail
  }
  function TPdfGaniTalbert.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
    var
      q: string;
      dict: TQueryDictionary;
    begin
      dict := createDBFieldList( db );

      if( update ) then
        begin
          q := 'DELETE FROM `inChartDetail` WHERE `chartID` = ' + intToStr( id );
          db.execute( q );
        end
      ;

      dict['chartType']   := '"Gani-Talbert"';
      dict['n']           := intToStr( n );
      dict['p']           := intToStr( k );

      if( not update ) then
        q := writeQuery( 'inChart', QInsert, dict )
      else
        q := writeQuery( 'inChart', QUpdate, dict, 'WHERE [chartID] = ' + intToStr( id ) );
      ;

      dict.Clear();
      freeAndNil( dict );

      db.execute( q );

      if( not update ) then id := db.lastInsertID();

      result := id;
    end
  ;
  {$ENDIF}

  {*
    Returns the point coordinates of self's discreet distribution
    @return histogram x,y coordinates for self from 1 to the maximum possible number of unique items
  }
  function TPdfGaniTalbert.getHistValues(): TMassHistogramValues;
    var
      j: integer;
      p: double;
      maxPossible: integer;
    begin
      if( nil = _histValues ) then
        begin
          _histValues := TMassHistogramValues.create();

          maxPossible := Math.min( n, k );

          // Start with 0 and work up to the maximum possible number of unique items,
          // calculating the probability of each outcome.
          for j := 1 to maxPossible do
            begin
              p := probDensity( j );
              _histValues.append( TMassHistogramValue.create( j, p ) );
              //dbcout2( usFloatToStr( p ) );
            end
          ;
        end
      ;

      //self.debug();
      //dbcout2( 'TPdfGaniTalbert.getHistValues...' );
      //dbcout2( 'Size: ' + intToStr( _histValues.Count ) );
      //_histValues.debug();

      buildCumulHistValues();

      result := _histValues;
    end
  ;


  {*
    Returns the probability of selecting j unique items, given n and k
    @param j the number of items selected at least once
    @param calcType the enumeration of the calculation algorithm
    @return the probability of selecting j unique items
    @comment TDBest will use the normal approximation algorithm when k > 30 otherwise the Stirling algorithm is used.
  }
  function TPdfGaniTalbert.probDensity( const j: integer; const calcType: TGaniTalbertAlgorithm = TDBest ): double;
    begin
      if( 0 = k ) then
        result := 0
      else
        begin
          case calcType of
            // TDBest is handled by the default case.
            TDStirling32:
              result := probStirling( j );
            TDExpMatrix:
              result := probExpMatrix( j );
            TDNormalApprox:
              result := probNormalApprox( j );
            TDStirlingGMP:
              result := probStirlingGmp( j );
            else
              begin
                if( 30 < k ) then
                  begin
                    //dbcout2( 'Big k: using Normal approx.' );
                    // FIX ME: Figure out when it's appropriate to use GMP version
                    //result := probStirlingGmp( j );
                    result := probNormalApprox( j );
                  end
                else
                  begin
                    //dbcout2( 'Little k: using Stirling' );
                    result :=  probStirling( j );
                  end
                ;
              end
            ;
          end;
        end
      ;
    end
  ;


  {*
    Calculates the mean, standard deviation, and variance, setting the values for _mean, _stDev, and _variance
    @comment If the variance is negative (when n is large) the value is set to 0.
  }
  procedure TPdfGaniTalbert.calcMeanAndVariance();
    (*
    var
      expr1: double;
      num, denom: double;
      _tVariance: double;
    *)
    begin
      if( 0 = k ) then
        begin
          _mean := 0.0;
          _variance := 0.0;
          _stdev := 0.0;
        end
      else
        begin
          _mean := n * ( 1 - power( ( ( n - 1 ) / n ), k ) );

          dbcout( 'GT mean = ' + usFloatToStr( _mean ), DBSHOWMSG );

          // Try calculating the exact variance, according to Gani.
          // If it is negative (it happens when n is large), then use
          // Talbert's approximation.  If that's negative too, assume that the
          // variance is 0.
          _variance :=
              ( n * (n-1) * power( (1-(2/n)), k ) )
            + ( n * power( (1-(1/n)), k ) )
            - ( power( n, 2 ) * power( (1-(1/n)), (2*k) ) )
          ;
          dbcout( 'Exact variance calculation = ' + usFloatToStr( _variance ), DBSHOWMSG );

          (*
              expr1 := power( ((n - 1 ) / n ), k );
              num := sqrt( n * k ) * power( ( 1 - expr1 ), 2 ) * expr1;
              denom := 2 - expr1;
              _tVariance := num / denom;
              dbcout( 'GT variance (approx ) = ' + usFloatToStr( _tVariance ), DBSHOWMSG );

          if( 0.0 > _variance ) then
            begin
              dbcout( 'Gani''s variance is negative: trying Talbert''s approximation.', DBSHOWMSG );
              _variance := _tVariance;
            end
          ;
          *)

          if( 0.0 > _variance ) then
            begin
              dbcout( 'Variance is negative: use 0.0', DBSHOWMSG  );
              _variance := 0.0;
            end
          ;

          dbcout( 'GT variance = ' + usFloatToStr( _variance ), DBSHOWMSG );

          _stdev := sqrt( _variance );

          dbcout( 'GT stddev = ' + usFloatToStr( _stdev ), DBSHOWMSG );
        end
      ;
    end
  ;

  {*
    Employs a quantile function to return a non-uniform random number
    @return a random number ranging from ? to ?
    @comment Fix Me: This description needs review! 
  }
  function TPdfGaniTalbert.rand(): double;
    var
      val: integer;
      r: double;
      upperLimit: integer;
    begin
      (*
      // FIX ME: Figure this out!
      // This handy little approximation will potentially save a ton of time.
      if( 10.0 < ( k / n ) ) then
        val := n
      else if( 10000 < k ) then // Use the Normal approximation
      *)
      if( 30 < k ) then // Use the Normal approximation
        begin
          //dbcout2( 'Big k: generating rand with normal approx.' );

          if( 0.0 < _stdev ) then
            begin
              val := -1;

              upperLimit := Math.min( n, k );

              while( ( 0 > val ) or ( upperLimit < val ) ) do
                begin
                  r := -1.0;
                  while( 0.0 > r ) do
                    r := TPdfGaussian.rand( _mean, _stdev )
                  ;
                  val := trunc( RoundDblTo( r, 0 ) );
                end
              ;
            end
          else
            val := trunc( _mean )
          ;
        end
      else // Use an "exact" calculation
        begin
          //dbcout2( 'Little k: using exact calculation' );
            val := invCumulative( rngRand() );
        end
      ;

      //dbcout2( val );

      result := val;
    end
  ;


  {*
     Helper function to probStirling(); returning a Stirling subset number,
     the number of ways to partition a set of i objects into j groups
     @param i a number between 0 and j
     @param j the number of items selected at least once in a set of n items
     @return the number of ways to partition a set of i objects into j groups
     @comment With values of k larger than about 30, this function produces a
     result with a large enough error that subsequent calculations will fail.

     @comment Fix Me: This description needs REVIEW!
  }
  function TPdfGaniTalbert.stirlingNumber( const i, j: integer ): double;
    begin
      // With values of k larger than about 30, this function produces
      // a result with a large enough error that subsequent calculations will
      // fail. This is true whether the logarithmic form (see below) is used or not.
      // Use of the 'straight' form seems to allow a slightly larger range of
      // k for which the function works, but the difference is not enough to
      // be of practical value.

      // Calculating without logs, native Delphi power function
      result := power( -1, i ) * choose( j, i ) * power( ( j - i ), k );

      //----------------------------------------------------------------------
      // KEEP THI CODE!
      //----------------------------------------------------------------------
      // Calculating without logs, GSL power function
      (*
      result := gslPower( -1, i ) * choose( j, i ) * gslPower( ( j - i ), k );
      *)

      // Calculating with logs
      (*
      if( 0 = ( j - i ) ) then
        result := 0
      else
        begin
          result := lnChoose( j, i );
          result := result + ( k * ln( j - i ) );
          result := exp( result );
          if( 1 = i mod 2 ) then
            result := -1 * result
          ;
        end
      ;
      *)
      //----------------------------------------------------------------------
    end
  ;


  {*
    Helper function for probDensity() when TGaniTalbertAlgorithm is TDStirlingGMP.
    Uses multiple precision arithmetic to calculate Stirling numbers.
    This approach should work for values far greater than 30.
    @param j the number of items selected at least once in a set of n items
    @return the probability of selecting j unique items, given n and k
  }
  function TPdfGaniTalbert.probStirlingGmp( const j: integer ): double;
    begin
      if( _libAphiGTLoaded ) then
        result := aphi_gani_talbert_probability( n, k, j )
      else
        begin
          raise exception.Create( 'Library not loaded in probStirlingGmp()' );
          result := 0.0;
        end
      ;
    end
  ;


  {*
    Helper function for probDensity() when TGaniTalbertAlgorithm is TDStirling32.
    This algorithm gives the same result as using the an exponent matrix,
    but tcomputation is far less intensive. Like the exponent matrix algorithm,
    this approach is tractable for values up to about 30.
    @param j the number of items selected at least once in a set of n items
    @return the probability of selecting j unique items, given n and k
  }
  function TPdfGaniTalbert.probStirling( const j: integer ): double;
    var
      i: integer;
      sum: double;
      expr: double;
    begin
      if( 0 = j ) then
        result := 0.0
      else
        begin
          sum := 0.0;

          for i := 0 to j do
            sum := sum + stirlingNumber( i, j )
          ;

          // With values of k larger than about 30, the function call above
          // has enough error that a negative sum will be generated.
          // This would be a bad thing: note that, in the code below,
          // ln( sum ) is calculated, which is not possible if sum is negative.

          // It seems to make no practical difference whether the 'straight'
          // or logarithmic forms of the code below are used: it's probably
          // all doing the same thing behind the scenes anyway.

          // Calculating with logs
          //----------------------

          expr := ln( sum ) - lnFactorial( j );
          expr := lnFactorial( n ) - lnFactorial( n - j ) + (( expr ));
          expr := ((expr)) - ( k * gslLog(n) );
          result := exp( expr );

          //----------------------------------------------------------------------
          // KEEP THIS CODE!
          //----------------------------------------------------------------------
          // Calculating without logs
          //-------------------------
          {*
          expr := sum / factorial( j );
          expr := ( factorial( n ) / factorial( n - j ) ) * ((expr));
          expr := ((expr)) / power( n, k );
          result := expr;
          *}
          //----------------------------------------------------------------------
        end
      ;
    end
  ;

  {*
    Helper method for createExpMatrix, it fills the matrix recursively and is
    part of the exponent matrix algorithm.
    @param arr the matrix to which completed arrays will be appended as columns
    @param v the array to fill
    @param arrPos array position
    @param val the initial value
    @param target the desired sum
  }
  procedure TPdfGaniTalbert.fillMatrixRecurse(
        arr: TTwoDIntArray;  // The matrix to which completed arrays will be appended as columns
        v: TQIntegerVector;     // The array to fill
        arrPos: integer;     // Array position
        val: integer;     // The initial value
        const target: integer // The desired sum
      );
    var
      i: integer;
      sum: integer;
      diff: integer;
    begin
      //dbcout2( 'fillMatrixRecurse() with val ' + intToStr( val ) );
      
      Application.ProcessMessages();
      
      if( arrPos > v.size - 1 ) then
        begin
          //dbcout2( 'No more spaces in the vector.' );
          exit;
        end
      ;

      v[arrPos] := val;

      sum := 0;
      for i := arrPos to v.size - 1 do
        sum := sum + v[i]
      ;

      if( target = sum ) then
        begin
          arr.appendCol( v );
          exit;
        end
      ;

      diff := target - sum;
      for i := diff downto 0 do
        begin
          fillMatrixRecurse(
            arr,
            v,
            ( arrPos + 1 ),
            i,
            diff
          );
        end
      ;

    end
  ;

  {*
    Creates and fills the exponent matrix. Part of the exponent matrix algorithm.
    @param j the number of items selected at least once in a set of n items
    @return an array of integers
    @comment Fix Me: a better description of the values in the matrix would be helpful
  }
  function TPdfGaniTalbert.createExpMatrix( const j: integer ): TTwoDIntArray;
    var
      v: TQIntegerVector;
      i: integer;
      target: integer;
      newMatrix: TTwoDIntArray;
    begin
      newMatrix := TTwoDIntArray.create( 0, j, 0, true );

      v := TQIntegerVector.create( j, 0 );

      target := k - j;

      for i := target downto 0 do
        begin
          fillMatrixRecurse(
            newMatrix, // The matrix to which completed arrays will be appended as columns
            v,         // The array to fill
            0,         // Array position
            i,         // The initial value
            target     // The desired sum
          );
        end
      ;

      v.Free();

      result := newMatrix;
    end
  ;


  {*
    Creates the sum of the exponent matrix. Part of the exponent matrix algorithm.

    @param j the number of items selected at least once in a set of n items
    @return sum of exponent matrix.  Fix Me: how could the result be better described??

    @comment Fix Me: This description really needs review and more detail!
  }
  function TPdfGaniTalbert.calcSum( expMatrix: TTwoDIntArray; const j: integer ): double;
    var
      prod: double;
      c, r: integer;
    begin
      result := 0.0;

      for c := 0 to expMatrix.nCols - 1 do
        begin
          prod := 1.0;
          for r := 0 to j - 1 do
            prod := prod * power( ( j - r ), expMatrix.getAt( c, r ) )
          ;
          result := result + prod;
        end
      ;
    end
  ;


  {*
    Helper function for probDensity() when TGaniTalbertAlgorithm is TDExpMatrix.
    @param j the number of items selected at least once in a set of n items
    @return the probability of selecting j unique items, given n and k
    @comment This method should not be used, except for demonstration, because it is slow.
  }
  function TPdfGaniTalbert.probExpMatrix( const j: integer ): double;
    var
      f: double;
      sum: double;
      num, denom: double;
      expMatrix: TTwoDIntArray;
    begin
      if( 0 = j ) then
        begin
          result := 0.0;
          exit;
        end
      ;

      expMatrix := createExpMatrix( j );

      sum := calcSum( expMatrix, j );

      // The formula:
      // result = [ ( n! / ( n-j)! ) * sum ] / [ n^k ]

      if( 0.0 = sum ) then
        begin
          result := 0.0;
          exit;
        end
      ;

      // Without logs
      //-------------
      (*
      f := factorial( n ) / factorial( n - j );
      num := f * sum;
      denom := power( n, k );
      result := num/denom;
      dbcout2( result );
      *)

      // With logs
      //----------
      f := lnFactorial( n ) - lnFactorial( n - j );
      num := f + ln( sum );
      denom := k * ln(n);
      result :=  exp( num - denom );

      expMatrix.free();
    end
  ;


   {*
    Helper function for probDensity() when TGaniTalbertAlgorithm is TDNormalApprox.
    @param j the number of items selected at least once in a set of n items
    @return the probability of selecting j unique items, given n and k
    @comment This method should not be used, except for demonstration, because it is slow.
  }
  function TPdfGaniTalbert.probNormalApprox( const j: integer ): double;
    begin
      if( 0.0 < _stdev ) then
        result := TPdfGaussian.probDensity( j, _mean, _stdev )
      else
        begin
          if( j = k ) then
            result := 1.0
          else
            result := 0.0
          ;
        end
      ;

      if( 1.0 < result ) then
        begin
          dbcout( 'TPdfGaniTalbert.probNormalApprox() is greater than 1', true );
          result := 1.0;
        end
      ;
    end
  ;

  {*
     For a given probability prob, the value that j will be at, or below.
     @param prob a probability between 0 and 1.0

     @comment Fix Me: description needs review!
  }
  function TPdfGaniTalbert.invCumulative( prob: double ): integer;
    begin
      if( ( 0.0 > prob ) or ( 1.0 < prob ) ) then
        begin
          raise exception.create( 'Invalid probability (' + usFloatToStr( prob ) + ') in TPdfDiscrete.invCumulative()' );
          result := 0;
        end
      // Remember the special case!
      else if( 0 = k ) then
        result := 0
      else
        result := inherited invCumulative( prob )
      ;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// DLL handling
//-----------------------------------------------------------------------------
  /// Sets the function pointer to nil
  procedure freeLoadedPointers();
    begin
      aphi_gani_talbert_probability := nil;
    end
  ;


  {*
     Attempts to load the library libaphi.dll
     @return true if the library is loaded sucessfully, else false
     @comment if the return value is false check _libAphiGTLoadErrors for error message
     @commment _libAphiGTLoaded is set by calling this function during unit initialization
  }
  function loadAPHIPointers(): boolean;
    var
      rngHandle: THandle; //Handle used to open the DLL.  Defined in unit Windows.

      lib_version: function(): pchar; cdecl;
      libVersion: string;
      lastError: integer;
    begin
      _libAphiGTLoadErrors := '';

      try
        setLastError( 0 );
        rngHandle := loadLibrary( 'libaphi.dll' );
      except
        result := false;
        lastError := getLastError();
        appendToString( _libAphiGTLoadErrors, 'Could not load library libaphi.dll: system error code ' + intToStr( lastError ) );
        freeLoadedPointers();
        exit;
      end;

      if(  rngHandle >= 32 ) then // libAPHI was successfully loaded.  Assign function pointers now.
        begin
          result := true;

          // Check the version of libAPHI
          //-----------------------------
          lib_version := getProcAddress( rngHandle, 'lib_version' );
          if( nil = @lib_version ) then
            begin
              appendToString( _libAphiGTLoadErrors, 'MISSING FUNCTION lib_version' );
              result := false;
              exit;
            end
          else
            begin
              libVersion := lib_version();
              if( REQUIRED_APHI_DLL_VERSION <> libVersion ) then
                begin
                  appendToString( _libAphiGTLoadErrors, 'INCOMPATIBLE DLL VERSION: ' + libVersion );
                  result := false;
                end
              ;
            end
          ;

          // GT probability function
          //------------------------
          aphi_gani_talbert_probability := getProcAddress( rngHandle, 'aphi_gani_talbert_probability' );
          if( nil = @aphi_gani_talbert_probability ) then
            begin
              appendToString( _libAphiGTLoadErrors, 'MISSING FUNCTION aphi_gani_talbert_probability' );
              result := false;
            end
          ;
        end
      else
        begin
          result := false;
          lastError := getLastError();
          appendToString( _libAphiGTLoadErrors, 'Could not load library libaphi.dll: system error code ' + intToStr( lastError ) );
        end
      ;

      dbcout( 'APHI rng loaded', DBSHOWMSG );
      dbcout( result, DBSHOWMSG );
    end
  ;
//-----------------------------------------------------------------------------

initialization
  _libAphiGTLoaded := loadAphiPointers();

end.



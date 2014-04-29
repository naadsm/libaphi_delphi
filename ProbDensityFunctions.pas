{*
ProbDensityFunctions.pas
-------------------------
Begin: 2005/03/15
Last revision: $Date: 2012-08-14 19:14:07 $ $Author: areeves $
Version number: $Revision: 1.54 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/delphi
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
Author: Shaun Case <Shaun.Case@colostate.edu>
--------------------------------------------------
Copyright (C) 2005 - 2012 Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
---------------------------------------------------
This unit defines the base classes for discreet and continuous PDFs and
a derived class for each supported type of distribution. The derived
class instance objects interact with the APHI and GSL libraries to fill each
with data that meet the shape and other characteristics of the distribution
type for the input parameters.
}

 (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator

    Base classes have been documented pretty well, the derived classes not so much
  *)

unit ProbDensityFunctions;

{$INCLUDE Defs.inc}

(*
  TODO:
    - All of the calculations should probably check for function validity before allowing use.
      (Consider setting a flag, so the validation function doesn't have to be rerun for
      every calculation.)
*)

interface

  uses
    // Standard Delphi units
    Math,
    SysUtils,
    Dialogs,

    // General purpose units
    {$IFDEF DATABASE_ENABLED}
    SqlClasses,
    {$ENDIF}

    Sdew,
    Points,
    MyDelphiArrayUtils,

    // QClasses for Delphi
    QIntegerMaps,
    QVectors,

    // APHI Modeling Library
    AphiRng,
    ChartFunction
  ;

  // REMEMBER: When new types are added, do all of the following:
  //   - Update the functions for returning sets of PDF types (continuousPdfs(), continuousBoundedPdfs(), etc.).
  //   - Update XML import capabilities:
  //     - Write an appropriate importPdfAbcXml() function.
  //     - Update NumXmlImportedPdfTypes and XmlImportedPdfTypes below.
  //   - Deal with PDF editing capabilities, especially in FrameFunctionParams2.

  type
  /// Enumeration of the types of Probability Density Functions supported
  TPdfType = (
    PdfUndefined,

    // Continuous types
    //-----------------
    PdfBeta,
    PdfBetaPERT,
    PdfExponential,
    PdfGamma,
    PdfGaussian,
    PdfHistogram,
    PdfInverseGaussian,     // needs CDF, inverse CDF
    PdfLogistic,
    PdfLogLogistic,
    PdfLognormal,
    PdfPareto,
    PdfPearson5,
    PdfPoint,
    PdfPiecewise,
    PdfTriangular,
    PdfUniform,
    PdfWeibull,

    // Discrete types
    //---------------
    PdfBernoulli,
    PdfBinomial,
    PdfDiscreteUniform,
    PdfHypergeometric,
    PdfNegativeBinomial,
    PdfPoisson,

    // Specialized types
    //------------------
    PdfSpGaniTalbert // A discrete type: see unit PdfGaniTalbert
  );

  type
  /// A set of objects of type TPdfType, often the set associated with a particular model
  TPdfTypeSet = set of TPdfType;

  /// Returns the set of continuous PDF types
  function continuousPdfs(): TPdfTypeSet;
  /// Returns the set of continuous bounded PDF types
  function continuousBoundedPdfs(): TPdfTypeSet;
  /// Returns the set of continuous and discrete bounded PDF types
  function boundedPdfs(): TPdfTypeSet;
  /// Returns the set of discrete PDF types
  function discretePdfs(): TPdfTypeSet;
  /// Returns the set of all PDF types
  function allPdfs(): TPdfTypeSet;

  /// Returns a string description for the specified type.
  function pdfTypeDescr( val: TPdfType ): string;

  /// Returns true if all the required DLLs are loaded successfully.
  function pdfFnsLoaded(): boolean;
  /// Returns true if the GSL library sucessfully loaded
  function pdfGslFnsLoaded(): boolean;
   /// Returns true if the APHI library sucessfully loaded
  function pdfAphiFnsLoaded(): boolean;

  /// Returns the error message, if any, generated during DLL loading.
  function pdfFnsLoadedErr(): string;

  {*
    For those who are inclined to forget how these functions work,
    the following R code may be useful:

      # Fun things to do with a standard Normal distribution (mean = 0, standard deviation = 1)
      x = seq( -1.96, 1.96, length = 10 );
      cumulProbDens = pnorm( x );
      cumulProbDens;
      # [1] 0.02499790 0.06369886 0.13810144 0.25677070 0.41380113 0.58619887 0.74322930 0.86189856 0.93630114 0.97500210
      probDens = dnorm( x );
      probDens;
      # [1] 0.05844094 0.12481733 0.22051754 0.32227155 0.38959322 0.38959322 0.32227155 0.22051754 0.12481733 0.05844094
      p = seq( 0.1, 0.9, length = 9 );
      invCumul = qnorm( p );
      invCumul;
      # [1] -1.2815516 -0.8416212 -0.5244005 -0.2533471  0.0000000  0.2533471  0.5244005  0.8416212  1.2815516
  *}

  // FIX ME: Consider making all of the calculations check for function validity before allowing use.
  type
  /// Base class for implementing discrete Probability Density Functions, inherits from TChartFunction
  TPdf = class( TChartFunction )
    protected
      _pdfType: TPdfType;        /// identifies the PDF type

      _infLeftTail: boolean;     /// whether the distribution has an infinite left tail
      _infRightTail: boolean;    /// whether the distribution has an infinite right tail

      _discreteVals: TQIntegerDoubleMap;    /// structure for the probabilities of a set of discrete values
      _posDiscreteVals: TQIntegerDoubleMap; /// structure for the probabilities of a set of positive discrete values, for distributions truncated at 0

      /// Returns the minimum of the distribution rounded to a discrete number
      function getMinD(): integer; virtual;
      /// Returns the maximum of the distribution rounded to a discrete number
      function getMaxD(): integer; virtual;

      /// Returns an object reference to _discreteVals
      function getDiscreteVals(): TQIntegerDoubleMap;
      /// Returns an object reference to _posDiscreteVals
      function getPosDiscreteVals(): TQIntegerDoubleMap;
      /// To be used by derived classes to initialize distribution parameters
      procedure initialize(); virtual; abstract;
      /// frees _discreteVals and _posDiscreteVals from memory
      procedure freeDiscreteValues();
       /// To be used by derived classes to initialize _discreteVals and _posDiscreteVals
      procedure initializeDiscreteVals(); virtual; abstract;
      /// Populates _posDiscreteVals from a subset of  _discreteVals values and ensures that the total probability will be 1.0
      procedure standardizePosDiscreteVals();

      {$IFDEF DATABASE_ENABLED }
      function createDBFieldList( db: TSqlDatabase ): TQueryDictionary;
      {$ENDIF}

      /// Get function for property pdfType
      function getPdfType() : TPdfType;
      /// Set function for property pdfType, setting it to val
      procedure setPdfType( val: TPdfType );

      /// always returns true
      function getYStartsAtZero(): boolean; override;

      /// Get function for property isDiscrete, to be implemented by derived class
      function getIsDiscrete(): boolean; virtual; abstract;
      /// Get function for property isContinuous
      function getIsContinuous(): boolean;

      /// Get function for property hasMean, to be implemented by derived class
      function getHasMean(): boolean; virtual; abstract;
      /// Get function for property mean, to be implemented by derived class
      function getMean(): double; virtual; abstract;

      /// Get function for property hasMin
      function getHasMin(): boolean;
      /// Get function for property hasMax
      function getHasMax(): boolean;
      /// Get function for property min, to be implemented by derived class
      function getMin(): double; virtual; abstract;
      /// Get function for property max, to be implemented by derived class
      function getMax(): double; virtual; abstract;

    public
      /// Creates and empty PDF and initializes it with safe but unuseful values (nil, false, 0 etc.)
      constructor create(); overload; override;
      /// Creates a PDF object by copying src
      constructor create( const src: TPdf ); overload;
      /// frees the object form memory
      destructor destroy(); override;
      /// Reads the unit element of the XML but derived class will have to handle their own parameters
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      /// Defines behavior for when err is not nil, but derived class must provide validation specifics for err
      function validate( err: pstring = nil ): boolean; override;
      /// Chart type, PDF type, chart name, and axes units of Chart are compared to self
      function compare( Chart:TChartFunction ):boolean; override;

      /// Returns the probability density at x for the distribution.  Same as GSL function gsl_ran_xxx_pdf or R function dxxx.
      function probDensity( x: double ): double; virtual; abstract;
      /// Returns the cumulative density at x for the distribution.  Same as GSL function gsl_cdf_xxx_P or R function pxxx.
      function cumulativeDistr( x: double ): double; virtual; abstract;
      /// Returns the inverse cumulative density at p for the distribution.  Same as GSL function gsl_cdf_xxx_Pinv or R function qxxx.
      function invCumulative( p: double ): double; virtual; abstract;
      /// Returns a random variate from the distribution.  Same as GSL function gsl_ran_xxx or R function rxxx.
      function rand(): double; virtual; abstract;

      /// Rounds the result of rand() to an integer.
      function randInt(): integer;

      /// Truncates the distribution at 0 returning a positive random number.
      function randNonNeg(): double;

      /// Truncates the distribution at 0, and rounds the result of rand() to an integer.
      function randNonNegInt(): integer;

      /// Truncates the distribution at 1, and rounds the result of rand() to an integer.
      function randPosInt(): integer;

      /// Returns the probability of x occurring
      function discreteProb( x: integer ): double;
      /// Returns the cummulative probability for x
      function discreteCumulativeProb( x: integer ): double;
      //function discreteInvCumulative( p: double ): integer; // See discrete PDFs for how to do this.

      /// Returns the probability of x occurring from the positive discrete distribution
      function posDiscreteProb( x: integer ): double;

      /// Outputs to the debug window the value returned from rand(), nSamples times
      procedure randTest( const nSamples: integer = 10 );

      property pdfType: TPdfType read getPdfType write setPdfType;  /// access to _pdfType

      property infiniteLeftTail: boolean read _infLeftTail;    /// read-only access to _infLeftTail
      property infiniteRightTail: boolean read _infRightTail;  /// read-only access to _infRightTail

      property isDiscrete: boolean read getIsDiscrete;       /// indicates whether PDF is a discrete type of function
      property isContinuous: boolean read getIsContinuous;   /// indicates whether PDF is a continuous type of function

      property mean: double read getMean;        /// provides the mean value for the PDF
      property hasMean: boolean read getHasMean; /// indicates whether the PDF has a mean value (points don't)

      property hasMin: boolean read getHasMin;   /// indicates whether the PDF has a minimum
      property hasMax: boolean read getHasMax;   /// indicates whether the PDF has a maximum

      property min: double read getMin;  /// provides the minimum value of the PDF
      property max: double read getMax;  /// provides the maximum value of the PDF

      property discreteVals: TQIntegerDoubleMap read getDiscreteVals;       /// read-only access to _discreteVals
      property posDiscreteVals: TQIntegerDoubleMap read getPosDiscreteVals; /// read-only access to _posDiscreetVals
      property minD: integer read getMinD;  /// provides the minimum discreet value
      property maxD: integer read getMaxD;  /// provides the maximum discreet value
    end
  ;


//-----------------------------------------------------------------------------
// Continuous PDF types
//-----------------------------------------------------------------------------
  type
  /// Base class for implementing continuous Probability Density Functions, inherits from TPdf
  TPdfContinuous = class( TPdf )
    protected
      /// Initialize _discreteVals and standardizes _posDiscreteVals
      procedure initializeDiscreteVals(); override;

      /// Get function for property isDiscrete of the base class, the result is false
      function getIsDiscrete(): boolean; override;

      /// Currently a stub that always returns false
      function getCanConvertToPiecewise(): boolean;

    public

      /// Creates a list of x,y coordinates of the PDF in 0.04 probability increments
      function createPointArray(): T2DPointList; virtual;

      /// Generates an APPROXIMATE (uses createPointArray) plain-text representation of the function.
      function asCsv(): string; override;

      /// Indicates whether the PDF can be converted to type PdfPiecewise
      property canConvertToPiecewise: boolean read getCanConvertToPiecewise;
    end
  ;


  type
  /// Used to decribe uncertainty of a probability or prevalence, especially about a binomial probability
  TPdfBeta = class( TPdfContinuous )
    protected
      _alpha1: double;
      _alpha2: double;
      _min: double;
      _max: double;

      procedure setAlpha1( val: double );
      procedure setAlpha2( val: double );
      procedure setMin( val: double );
      procedure setMax( val: double );
      function getAlpha1(): double;
      function getAlpha2(): double;
      function getMin(): double; override;
      function getMax(): double; override;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;
      
      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfBeta ); overload;
      constructor create( alpha1, alpha2, min, max: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;

      property alpha1: double read getAlpha1;
      property alpha2: double read getAlpha2;
    end
  ;


  {*
    FIX ME: is this comment still correct??

    NAADSM uses a set of approximations (based on Davis) for the BetaPERT PDF
    which are slightly different from those described by Vose and Palisade.
    1) Davis, R.E. http://www.cob.sjsu.edu/facstaff/davis_r/courses/QBAreader/QBAtoc.html
    2) Vose, D. 1996.  Quantitative risk analysis: a guide to Monte Carlo Simulation Modelling.
       John Wiley and Sons, New York.
    3) Palisade Corporation. 2002.  A concise summary of @RISK probability distribution functions.
       Available (probably illegally) from http://project.zf.jcu.cz/risk/data/distfunc.pdf
  }
  type
  /// For modeling expert opinion, having a minimum, most likely and maximum estimates.
  TPdfBetaPERT = class( TPdfContinuous )
    protected
      _a: double;
      _c: double;
      _b: double;

      procedure setA( val: double );
      procedure setC( val: double );
      procedure setB( val: double );
      function getA(): double;
      function getC(): double;
      function getB(): double;

      function getMin(): double; override;
      function getMax(): double; override;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;
      
      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfBetaPERT ); overload;
      constructor create( min, mode, max: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;
      
      procedure debug(); override;

      property mode: double read getC;
    end
  ;


  type
  /// Models the time until the occurrence of a first event in a Poisson process
  TPdfExponential = class( TPdfContinuous )
    protected
      _mean: double;

      function getMean(): double; override;
      procedure setMean( val: double );
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

      function getDescr(): string; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfExponential ); overload;
      constructor create( mean: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;
    end
  ;


  type
  /// Used in modeling aspects of poisson processes
  TPdfGamma = class( TPdfContinuous )
    protected
      _alpha: double;
      _beta: double;

      procedure setAlpha( val: double );
      procedure setBeta( val: double );
      function getAlpha(): double;
      function getBeta(): double;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfGamma ); overload;
      constructor create( alpha, beta: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;

      property alpha: double read getAlpha;
      property beta: double read getBeta;
    end
  ;


  type
  /// The Normal distribution
  TPdfGaussian = class( TPdfContinuous )
    protected
      _mean: double;
      _stddev: double;

      function getMean(): double; override;
      function getStdDev(): double;
      procedure setMean( val: double );
      procedure setStdDev( val: double );

      function getMin(): double; override;
      function getMax(): double; override;

      function getHasMean(): boolean; override;
      function getDescr(): string; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfGaussian ); overload;
      constructor create( mean, stddev: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; overload; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; overload; override;

      class function probDensity( const x, mean, stdev: double ): double; reintroduce; overload;
      class function rand( const mean, stdev: double ): double; reintroduce; overload;

      procedure debug(); override;

      property stddev: double read getStdDev;
    end
  ;


  type
  /// Defines the characteristics of a point of a histogram
  RHistogramPoint = record
    range: double;   /// the bin length
    count: double;   /// the number of occurrences in range
    density: double; /// the proportion of occurrences in range
  end;

  type
  /// Structure to delimit a histogram
  RHistogramPointArray = array of RHistogramPoint;

  type
  /// Technique for replicating the distribution (shape) of a large set of data
  TPdfHistogram = class( TPdfContinuous )
    protected
      _cStructPtr: pointer;
      _ranges: TARDoubleArray;
      _counts: TARDoubleArray;
      _densities: TARDoubleArray;


      procedure standardizeArea();

      function getMean(): double; override;

      function getMin(): double; override;
      function getMax(): double; override;

      function getHasMean(): boolean; override;
      function getDescr(): string; override;

      function getXCount(): integer;
      function getYCount(): integer;
      function getNBins(): integer;

      function xyValsOK( msg: Pstring = nil ): boolean;

      procedure initialize(); override;
      procedure setRangesAndCounts( ranges, counts: TARDoubleArray );
      procedure processRangesAndDensities();

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfHistogram ); overload;
      constructor create( ranges, counts: TARDoubleArray; u: TChartUnitType = UUnknown ); overload;
      constructor create(  rcd: RHistogramPointArray; u: TChartUnitType = UUnknown ); overload;

      {$IFDEF DATABASE_ENABLED}
      constructor create( db: TSqlDatabase; chartID: integer; u: TChartUnitType = UUnknown ); overload;
      {$ENDIF}
      destructor destroy(); override;

      class function standardize( var rcd: RHistogramPointArray; msg: PString = nil ): boolean;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      { Generates a plain-text representation of the function. }
      function asCsv(): string; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createPointArray(): T2DPointList; override;
      function createHistogramPointArray(): RHistogramPointArray;

      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; overload; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; overload; override;

      procedure debug(); override;

      property mean: double read getMean;
      property nBins: integer read getNBins;
      property xCount: integer read getXCount;
      property yCount: integer read getYCount;
    end
  ;


  type
  /// Has been used in place of Lognormal distribution having too heavy a right tail
  TPdfInverseGaussian = class( TPdfContinuous )
    protected
      _mean: double;
      _shape: double;

      function getMean(): double; override;
      function getShape(): double;
      procedure setMean( val: double );
      procedure setShape( val: double );

      function getMin(): double; override;
      function getMax(): double; override;

      function getHasMean(): boolean; override;
      function getDescr(): string; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfInverseGaussian ); overload;
      constructor create( mean, shape: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;

      property shape: double read getShape;
    end
  ;


  type
  /// Looks similar to the Gaussian distribution but has a higher kurtosis
  TPdfLogistic = class( TPdfContinuous )
    protected
      _location: double;
      _scale: double;

      procedure setLocation( val: double );
      procedure setScale( val: double );
      function getLocation(): double;
      function getScale(): double;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfLogistic ); overload;
      constructor create( location, scale: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;

      property location: double read getLocation;
      property scale: double read getScale;
    end
  ;


  type
  /// When Log(X) takes a Logistic distribution then X takes a LogLogistic distribution
  TPdfLogLogistic = class( TPdfContinuous )
    protected
      _location: double;
      _scale: double;
      _shape: double;

      procedure setLocation( val: double );
      procedure setScale( val: double );
      procedure setShape( val: double );
      function getLocation(): double;
      function getScale(): double;
      function getShape(): double;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfLogLogistic ); overload;
      constructor create( location, scale, shape: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;

      property location: double read getLocation;
      property scale: double read getScale;
      property shape: double read getShape;
    end
  ;


  type
  /// Often a good representation for a physical quantity that extends from zero to + infinity and is positively skewed
  TPdfLognormal = class( TPdfContinuous )
    protected
      _mean: double;
      _stddev: double;

      _zeta: double;
      _sigma: double;

      procedure setMean( val: double );
      procedure setStddev( val: double );
      procedure setZeta( val: double );
      procedure setSigma( val: double );
      function getZeta(): double;
      function getSigma(): double;
      function getMean(): double; override;
      function getStddev(): double;

      function getMin(): double; override;
      function getMax(): double; override;

      function getHasMean(): boolean; override;
      function getDescr(): string; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfLognormal ); overload;
      constructor create( mean, stddev: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;

      class function calculateZeta( mean, stddev: double ): double;
      class function calculateSigmaPrime( mean, stddev: double ): double;
      class function calculateMean( zeta, sigma: double ): double;
      class function calculateStddev( zeta, sigma: double ): double;

      property zeta: double read getZeta write setZeta;
      property sigma: double read getSigma write setSigma;
      property mean: double read getMean;
      property stddev: double read getStddev;
    end
  ;


  type
  /// An exponential shape, right skewed where mode and minimum are equal
  TPdfPareto = class( TPdfContinuous )
    protected
      _theta: double;
      _a: double;

      procedure setTheta( val: double );
      procedure setA( val: double );
      function getTheta(): double;
      function getA(): double;

      function getMin(): double; override;
      function getMax(): double; override;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfPareto ); overload;
      constructor create( theta, a: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;

      property theta: double read getTheta;
      property a: double read getA;
    end
  ;


  type
  /// Domain: 0 <= x < +infinity; Restrictions: alpha > 0, beta > 0 
  TPdfPearson5 = class( TPdfContinuous )
    protected
      _alpha: double;
      _beta: double;

      procedure setAlpha( val: double );
      procedure setBeta( val: double );
      function getAlpha(): double;
      function getBeta(): double;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfPearson5 ); overload;
      constructor create( alpha, beta: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;

      property alpha: double read getAlpha;
      property beta: double read getBeta;

    end
  ;


  type
  /// User-defined
  TPdfPiecewise = class( TPdfContinuous )
    protected
      _points: RPointArray;
      _slope: TARDoubleArray;
      _cumul: TARDoubleArray;

      function getPointCount(): integer;
      function getPoint( i: integer ): RPoint;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

      procedure initialize(); override;

      procedure calculateSlopeAndCumul();

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfPiecewise ); overload;
      constructor create( pnts: RPointArray; u: TChartUnitType = UUnknown ); overload;
      constructor create( pnts: T2DPointList; u: TChartUnitType = UUnknown ); overload;
      {$IFDEF DATABASE_ENABLED}
      constructor create( db: TSqlDatabase; chartID: integer; u: TChartUnitType = UUnknown ); reintroduce; overload;
      {$ENDIF}
      destructor destroy(); override;

      //procedure addPoint( x, y: double ); // AR 10/30/09 This function isn't really safe.  Make it go away.
      procedure clearPoints();

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function createCopy(): TChartFunction; override;

      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function probDensity( x: double ): double; override;
      function rand(): double; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createRecordPointArray(): RPointArray;
      function createPointArray(): T2DPointList; override;

      class function standardize( var points: RPointArray; msg: PString = nil ): boolean;

      property pointCount: integer read getPointCount;
      property points[i: integer]: RPoint read getPoint;

      //property pointArray: RPointArray read _points;

      procedure debug(); override;
    end
  ;


  type
  /// Point estimate, rather than a distribution
  TPdfPoint = class( TPdfContinuous )
    protected
      _point: double;

      procedure setPoint( val: double );
      function getPoint(): double;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

      function getMinD(): integer; override;
      function getMaxD(): integer; override;

      procedure initialize(); override;

      procedure initializeDiscreteVals(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfPoint ); overload;
      constructor create( c: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function invCumulative( p: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function probDensity( x: double ): double; override;
      function rand(): double; override;

      function createPointArray(): T2DPointList; override;

      procedure debug(); override;

      property mode: double read getPoint;
      property value: double read getPoint;
    end
  ;


  type
  /// For modeling expert opinion, a triangular shape having a minimum, most likely and maximum estimates.
  TPdfTriangular = class( TPdfContinuous )
    protected
      _a: double;
      _c: double;
      _b: double;

      procedure setA( val: double );
      procedure setC( val: double );
      procedure setB( val: double );
      function getA(): double;
      function getC(): double;
      function getB(): double;
      function getMin(): double; override;
      function getMax(): double; override;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfTriangular ); overload;
      constructor create( a, c, b: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      function createPointArray(): T2DPointList; override;

      procedure debug(); override;

      property a: double read getA;
      property c: double read getC;
      property b: double read getB;
      property mode: double read getC;
    end
  ;


  type
  /// Assigns equal probability to all values between min and max
  TPdfUniform = class( TPdfContinuous )
    protected
      _a: double;
      _b: double;

      procedure setA( val: double );
      procedure setB( val: double );
      function getA(): double;
      function getB(): double;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfUniform ); overload;
      constructor create( a, b: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; overload; override;
      class function rand( const min, max: double ): double; reintroduce; overload;

      function createPointArray(): T2DPointList; override;

      procedure debug(); override;

      property a: double read getA;
      property b: double read getB;
    end
  ;


  type
  /// Domain: -infinity < x < +infinity; Restrictions: alpha > 0, beta > 0
  TPdfWeibull = class( TPdfContinuous )
    protected
      _alpha: double;
      _beta: double;

      procedure setAlpha( val: double );
      procedure setBeta( val: double );
      function getAlpha(): double;
      function getBeta(): double;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

      procedure initialize(); override;

    public
      constructor create(); overload; override;
      constructor create( const src: TPdfWeibull ); overload;
      constructor create( alpha, beta: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( Chart:TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( x: double ): double; override;
      function cumulativeDistr( x: double ): double; override;
      function invCumulative( p: double ): double; override;
      function rand(): double; override;

      procedure debug(); override;

      property alpha: double read getAlpha;
      property beta: double read getBeta;

    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Discrete PDF types
//-----------------------------------------------------------------------------
  type
  /// Point data for discete PDF types
  TMassHistogramValue = class( T2DPoint );


  type
  /// Structure to hold the point data of discete PDF types
  TMassHistogramValues = class( T2DPointList )
    public
      function at( const idx: integer ): TMassHistogramValue;
    end
  ;


  type
  /// Base class for discrete probability density functions.
  TPdfDiscrete = class( TPdf )
    protected
      _histValues: TMassHistogramValues;       /// container for point values of the histogram
      _cumulHistValues: TMassHistogramValues;  /// container point values of the cummulative histogram

      /// Uses massHistogramValues for the initialization
      procedure initializeDiscreteVals(); override;

      /// frees private member histogram data data structures
      procedure freeHistogramValues();

      /// Get function for property isDiscrete, always returns true
      function getIsDiscrete(): boolean; override;
      /// Get function for property massHistogramVaalues, upto the derived class to implement
      function getHistValues(): TMassHistogramValues; virtual; abstract;
      /// Get function for property cumulHistogramValues
      function getCumulHistValues(): TMassHistogramValues;
      /// Rebuilds  _cumulHistValues based on the current data of _histValues
      procedure buildCumulHistValues();
      
    public
      /// Creates a PdfDiscrete object with private members set to nil
      constructor create(); override;
      /// Calls freeHistogramValues and then frees itself from memory
      destructor destroy(); override;

      /// Generates an APPROXIMATE plain-text representation of this function.
      function asCsv(): string; override;

      /// Returns the cummulative probabily value for k
      function cumulativeDistr( k: integer; raiseExceptionOnFailure: boolean = true ): double; reintroduce;

      /// Returns inverse cumulative function for probability prob.
      function invCumulative( prob: double ): integer; reintroduce; virtual;

      /// Generates a random variate using the "inverse cumulative" function.
      function randInvCumulative(): integer;

      /// Read-only access to _histValues
      property massHistogramValues: TMassHistogramValues read getHistValues;
      /// Read-only access to _cumulHistValues
      property cumulHistogramValues: TMassHistogramValues read getCumulHistValues;
    end
  ;


  type
  /// Used to model an event occurring or not.
  TPdfBernoulli = class( TPdfDiscrete )
    protected
      _p: double;

      procedure setP( val: double );

      procedure initialize(); override;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

    public
      constructor create(); overload; override;
      constructor create( p: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      constructor create( const src: TPdfBernoulli ); overload;
      destructor destroy(); override;

      procedure setParams( p: double; u: TChartUnitType = UUnknown );

      function validate( err: pstring = nil ): boolean; override;
      function compare( chart: TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( k: integer ): double; reintroduce;
      function rand(): double; overload; override;

      // FIX ME: If this is protected (as it should be), I sometimes get 'abstract' errors.
      function getHistValues(): TMassHistogramValues; override;

      class function rand( const p: double ): double; reintroduce; overload;

      procedure debug(); override;

      property p: double read _p write setP;
    end
  ;


  type
  /// Models the number of successes from n independent trials where there is a probability p of success in each trial.
  TPdfBinomial = class( TPdfDiscrete )
    protected
      _n: integer;
      _p: double;

      procedure setN( val: integer );
      procedure setP( val: double );

      procedure initialize(); override;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

    public
      constructor create(); overload; override;
      constructor create( n: integer; p: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      constructor create( const src: TPdfBinomial ); overload;
      destructor destroy(); override;

      procedure setParams( n: integer; p: double; u: TChartUnitType = UUnknown );

      function validate( err: pstring = nil ): boolean; override;
      function compare( chart: TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( k: integer ): double; reintroduce;
      function rand(): double; overload; override;

      // FIX ME: If this is protected (as it should be), I sometimes get 'abstract' errors.
      function getHistValues(): TMassHistogramValues; override;

      class function rand( const n: integer; p: double ): double; reintroduce; overload;

      // For lack of a better place at the moment, multinomial functions are included here
      class procedure multinomialRand( n: integer; const p: TQDoubleVector; output: TQIntegerVector );

      procedure debug(); override;

      property n: integer read _n write setN;
      property p: double read _p write setP;
    end
  ;


  type
  /// Describes a variable that can take one of several explicit discrete values with equal probabilities of taking any particular value.
  TPdfDiscreteUniform = class( TPdfDiscrete )
    protected
      _min: integer;
      _max: integer;

      procedure initialize(); override;

      procedure setMin( val: integer );
      procedure setMax( val: integer );
      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

    public
      constructor create(); overload; override;
      constructor create( min, max: integer; u: TChartUnitType = UUnknown ); reintroduce; overload;
      constructor create( const src: TPdfDiscreteUniform ); overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( chart: TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( k: integer ): double; reintroduce;
      function rand(): double; override;

      // FIX ME: If this is protected (as it should be), I sometimes get 'abstract' errors.
      function getHistValues(): TMassHistogramValues; override;

      procedure debug(); override;
    end
  ;


 type
 {*
    Models the number of items of a particular type that there are in a sample of size n
    where that sample is drawn from a population of size m of which d are also of that particular type.
  }
 TPdfHypergeometric = class( TPdfDiscrete )
    protected
      _m: integer; // Population size
      _d: integer; // Size of subpopulation of interest
      _n: integer; // Sample size

      procedure initialize(); override;

      procedure setM( val: integer );
      procedure setD( val: integer );
      procedure setN( val: integer );

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

    public
      constructor create(); overload; override;
      constructor create( n, d, m: integer; u: TChartUnitType = UUnknown ); reintroduce; overload;
      constructor create( const src: TPdfHypergeometric ); overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( chart: TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( k: integer ): double; reintroduce;
      function rand(): double; overload; override;

      class function rand( const n, d, m: integer ): double; reintroduce; overload;

      // For lack of a better place at the moment, multivariate hypergeometric functions are included here
      class procedure multiHypergeometricRand( n: integer; const d: TQIntegerVector; output: TQIntegerVector );

      // FIX ME: If this is protected (as it should be), I sometimes get 'abstract' errors.
      function getHistValues(): TMassHistogramValues; override;

      procedure debug(); override;

      property m: integer read _m;
      property d: integer read _d;
      property n: integer read _n;
    end
  ;


  type
  /// Estimates the number of failures there will be before s successes are achieved where there is a probability p of success with each trial.
  TPdfNegativeBinomial = class( TPdfDiscrete )
    protected
      _s: integer;
      _p: double;

      procedure setS( val: integer );
      procedure setP( val: double );

      procedure initialize(); override;

      function getDescr(): string; override;
      function getMean(): double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

    public
      constructor create(); overload; override;
      constructor create( s: integer; p: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      constructor create( const src: TPdfNegativeBinomial ); overload;
      destructor destroy(); override;

      procedure setParams( s: integer; p: double; u: TChartUnitType = UUnknown );

      function validate( err: pstring = nil ): boolean; override;
      function compare( chart: TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( k: integer ): double; reintroduce;
      function rand(): double; override;

      // FIX ME: If this is protected (as it should be), I sometimes get 'abstract' errors.
      function getHistValues(): TMassHistogramValues; override;

      procedure debug(); override;

      property s: integer read _s write setS;
      property p: double read _p write setP;
    end
  ;


  type
  {*
    Models the number of occurrences of an event in a given time with an expected rate
    of 1/mean events when the time between successive events follows a Poisson process.
  }
  TPdfPoisson = class( TPdfDiscrete )
    protected
      _mean: double;

      procedure initialize(); override;

      procedure setMean( val: double );
      function getDescr(): string; override;
      function getMean: double; override;
      function getHasMean(): boolean; override;

      function getMin(): double; override;
      function getMax(): double; override;

    public
      constructor create(); overload; override;
      constructor create( mean: double; u: TChartUnitType = UUnknown ); reintroduce; overload;
      constructor create( const src: TPdfPoisson ); overload;
      destructor destroy(); override;

      function validate( err: pstring = nil ): boolean; override;
      function compare( chart: TChartFunction ):boolean; override;

      function ssXml( const indent: integer ): string; override;
      function importXml( element: Pointer; sdew: TSdew ): boolean; override;
      {$IFDEF DATABASE_ENABLED}
      function populateDatabase( db: TSqlDatabase; update: boolean = false ): integer; override;
      {$ENDIF}
      function createCopy(): TChartFunction; override;

      function probDensity( k: integer ): double; reintroduce;
      function rand(): double; override;

      // FIX ME: If this is protected (as it should be), I sometimes get 'abstract' errors.
      function getHistValues(): TMassHistogramValues; override;

      procedure debug(); override;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// XML import of probability density functions
//-----------------------------------------------------------------------------
  // The function that does all of the heavy lifting.
  /// The primary XML parsing function, sdew has already consumed the XML and will now try to parse out the PDF
  function createPdfFromXml( element: pointer; sdew: TSdew; unused: pointer ): TObject; overload;
  /// Overloaded version to save you from having to pass in a nil for unused pointer
  function createPdfFromXml( element: pointer; sdew: TSdew ): TPdf; overload;
  /// Creates the PDF from string xml that holds the contents of a probability-density-function XML element
  function createPdfFromXml( xml: string ): TPdf; overload;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Global helper functions
//-----------------------------------------------------------------------------
  {*
   Converts a SHARCSpread XML/database string description into the appropriate ordinal type.
  }
  function pdfType( descr: string ): TPdfType;

  /// FIX ME: I'm not sure if this function is still necessary.
  function createPdfFromBaseObject( obj: TPdf ): TChartFunction;

  (**
  // Two local functions are used in this unit, and are called in the initialization section.
  // These functions handle the assignment of function pointers to functions available in
  // dynamically loaded libraries.  The function pointers themselves are declared below
  // in the implementation section.

  { Initializes function pointers to nil.  }
  procedure freeLoadedPointers();

  { Assigns function pointers from the DLLs.  If successful, sets libAphiLoaded and gslLoaded to true. }
  function loadGslPointers(): boolean;
  **)
//-----------------------------------------------------------------------------


implementation

  uses
    // Standard Delphi units
    Windows, // Defines THandle
    TypInfo,
    Types,

    // General purpose units
    FunctionPointers,
    myStrUtils,
    DebugWindow,
    I88n,
    ARMath,
    ARMathAdvanced,
    RoundToXReplacement_3c

    {$IFDEF DATABASE_ENABLED}
    , FunctionEnums
    {$ENDIF}

    , XMLReader
  ;


  const
    REQUIRED_APHI_DLL_VERSION: string = '0.5';  /// Expected version number of the APHI library

    DBSHOWMSG: boolean = false; /// Set to true to enable debugging messages for this unit
    DBSHOWMSGSELECTED: boolean = false; /// Set to true to enable selected debugging messages for this unit.

  var
    dllLoadErrMsgs: string;  /// If there were library loading errors, find them here
    libAphiLoaded: boolean;  /// True if APHI library loaded sucessfully, else false
    gslLoaded: boolean;      /// True if GSL library loaded sucessfully, else false

    //-------------------------------------------------------------------------
    // Continuous types
    //-------------------------------------------------------------------------
    // Beta PDF
    //---------
    aphi_beta_pdf: TCFnDouble_5_Double_Double_Double_Double_Double;
    aphi_beta_cdf: TCFnDouble_5_Double_Double_Double_Double_Double;
    aphi_beta_inverse_cdf: TCFnDouble_5_Double_Double_Double_Double_Double;
    aphi_beta_rand: TCFnDouble_5_Cardinal_Double_Double_Double_Double;

    // BetaPERT pdf
    //-------------
    aphi_beta_pert_pdf: TCFnDouble_4_Double_Double_Double_Double;
    aphi_beta_pert_cdf: TCFnDouble_4_Double_Double_Double_Double;
    aphi_beta_pert_inverse_cdf: TCFnDouble_4_Double_Double_Double_Double;
    aphi_beta_pert_rand: TCFnDouble_4_Cardinal_Double_Double_Double;

    // Exponential PDF
    //----------------
    gslRanExponentialPdf: TCFnDouble_2_Double_Double;
    gslCdfExponentialP: TCFnDouble_2_Double_Double;
    gslCdfExponentialPinv: TCFnDouble_2_Double_Double;
    gslRanExponential: TCFnDouble_2_Cardinal_Double;

    // Gamma PDF
    //----------
    gslRanGammaPdf: TCFnDouble_3_Double_Double_Double;
    gslCdfGammaP: TCFnDouble_3_Double_Double_Double;
    gslCdfGammaPinv: TCFnDouble_3_Double_Double_Double;
    gslRanGamma: TCFnDouble_3_Cardinal_Double_Double;

    // Gaussian PDF
    //--------------
    gslRanGaussianPdf: TCFnDouble_2_Double_Double;
    gslCdfGaussianP: TCFnDouble_2_Double_Double;
    gslCdfGaussianPinv: TCFnDouble_2_Double_Double;
    gslRanGaussian: TCFnDouble_2_Cardinal_Double;

    // Histogram PDF
    //--------------
    aphi_create_histogram: TCFnPointer_3_Int_DArr_DArr;
    aphi_free_histogram: TCFnVoid_1_Pointer;
    aphi_histogram_mean: TCFnDouble_1_Pointer;
    aphi_histogram_pdf: TCFnDouble_2_Double_Pointer;
    aphi_histogram_cdf: TCFnDouble_2_Double_Pointer;
    aphi_histogram_inverse_cdf: TCFnDouble_2_Double_Pointer;
    aphi_histogram_rand: TCFnDouble_2_Cardinal_Pointer;

    // Inverse Gaussian PDF
    //---------------------
    aphi_inverse_gaussian_pdf: TCFnDouble_3_Double_Double_Double;
    aphi_inverse_gaussian_cdf: TCFnDouble_3_Double_Double_Double;
    aphi_inverse_gaussian_inverse_cdf: TCFnDouble_3_Double_Double_Double;

    // Logistic PDF
    //-------------
    gslRanLogisticPDF: TCFnDouble_2_Double_Double;
    gslCdfLogisticP: TCFnDouble_2_Double_Double;
    gslCdfLogisticPinv: TCFnDouble_2_Double_Double;
    gslRanLogistic: TCFnDouble_2_Cardinal_Double;

    // Loglogistic PDF
    //----------------
    aphi_loglogistic_pdf: TCFnDouble_4_Double_Double_Double_Double;
    aphi_loglogistic_cdf: TCFnDouble_4_Double_Double_Double_Double;
    aphi_loglogistic_inverse_cdf: TCFnDouble_4_Double_Double_Double_Double;
    aphi_loglogistic_rand: TCFnDouble_4_Cardinal_Double_Double_Double;

    // Lognormal PDF
    //--------------
    gslRanLognormalPDF: TCFnDouble_3_Double_Double_Double;
    gslCdfLognormalP: TCFnDouble_3_Double_Double_Double;
    gslCdfLognormalPinv: TCFnDouble_3_Double_Double_Double;
    gslRanLognormal: TCFnDouble_3_Cardinal_Double_Double;
    
    // Pareto PDF
    //-------------
    gslRanParetoPdf: TCFnDouble_3_Double_Double_Double;
    gslCdfParetoP: TCFnDouble_3_Double_Double_Double;
    gslCdfParetoPinv: TCFnDouble_3_Double_Double_Double;
    gslRanPareto: TCFnDouble_3_Cardinal_Double_Double;

    // Pearson5 PDF
    //-------------
    aphi_pearson5_pdf: TCFnDouble_3_Double_Double_Double;
    aphi_pearson5_cdf: TCFnDouble_3_Double_Double_Double;
    aphi_pearson5_inverse_cdf: TCFnDouble_3_Double_Double_Double;
    aphi_pearson5_rand: TCFnDouble_3_Cardinal_Double_Double;

    // Piecewise PDF
    //--------------
    gsl_poly_solve_quadratic: TCFnInt_5_Double_Double_Double_PDouble_PDouble;

    // Triangular PDF
    //---------------
    aphi_triangular_cdf: TCFnDouble_4_Double_Double_Double_Double;
    aphi_triangular_inverse_cdf: TCFnDouble_4_Double_Double_Double_Double;
    aphi_triangular_rand: TCFnDouble_4_Cardinal_Double_Double_Double;

    // Uniform PDF
    //------------
    gslCdfUniformP: TCFnDouble_3_Double_Double_Double;
    gslCdfUniformPinv: TCFnDouble_3_Double_Double_Double;
    gslRanUniform: TCFnDouble_3_Cardinal_Double_Double;

    // Weibull PDF
    //------------
    gslRanWeibullPDF: TCFnDouble_3_Double_Double_Double;
    gslCdfWeibullP: TCFnDouble_3_Double_Double_Double;
    gslCdfWeibullPinv: TCFnDouble_3_Double_Double_Double;
    gslRanWeibull: TCFnDouble_3_Cardinal_Double_Double;

    //-------------------------------------------------------------------------
    // Discrete types
    //-------------------------------------------------------------------------
    // Binomial PDF
    //-------------
    gslRanBinomialPdf: TCFnDouble_3_Int_Double_Int;
    gslRanBinomial: TCFnInt_3_Cardinal_Double_Int;

    // Hypergeometric PDF
    //-------------------
    gslRanHypergeometricPdf: TCFnDouble_4_Int_Int_Int_Int;
    gslRanHypergeometric: TCFnInt_4_Cardinal_Int_Int_Int;

    // NegativeBinomial PDF
    //-------------
    gslRanNegativeBinomialPdf: TCFnDouble_3_Int_Double_Double;
    gslRanNegativeBinomial: TCFnInt_3_Cardinal_Double_Double;

    // Poisson PDF
    //------------
    gslRanPoissonPdf: TCFnDouble_2_Int_Double;
    gslRanPoisson: TCFnInt_2_Cardinal_Double;
    //-------------------------------------------------------------------------


  
//-----------------------------------------------------------------------------
// Global DLL management functions
//-----------------------------------------------------------------------------
  function pdfFnsLoaded(): boolean;
    begin
      result := libAphiLoaded and gslLoaded;
    end
  ;

  function pdfGslFnsLoaded(): boolean;
    begin
      result := gslLoaded;
    end
  ;

  function pdfAphiFnsLoaded(): boolean;
    begin
      result := libAphiLoaded;
    end
  ;

  function pdfFnsLoadedErr(): string;
    begin
      result := trim( dllLoadErrMsgs );
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Global helper functions
//-----------------------------------------------------------------------------
  function continuousPdfs(): TPdfTypeSet;
    begin
      result := [ PdfBeta..PdfWeibull ];
    end
  ;

  function continuousBoundedPdfs(): TPdfTypeSet;
    begin
      result := [
        PdfBeta,
        PdfBetaPERT,
        PdfHistogram,
        PdfPiecewise,
        PdfPoint,
        PdfTriangular,
        PdfUniform,
        PdfHistogram
      ];
    end
  ;

  function boundedPdfs(): TPdfTypeSet;
    begin
      result := continuousBoundedPdfs() +
        [
          PdfBernoulli,
          PdfBinomial,
          PdfHypergeometric,
          PdfDiscreteUniform
        ]
      ;
    end
  ;

  function discretePdfs(): TPdfTypeSet;
    begin
      result := [ PdfBernoulli..PdfPoisson ];
    end
  ;

  function allPdfs(): TPdfTypeSet;
    begin
      result := [ PdfBeta..PdfPoisson ];
    end
  ;


  function createPdfFromBaseObject( obj: TPdf ): TChartFunction;
    var
      newobj: TChartFunction;
    begin
      dbcout( '*** PDF type ' + intToStr( ord( obj.pdfType ) ), DBSHOWMSG );
      dbcout( pdfTypeDescr( obj.pdfType ), DBSHOWMSG );

      case obj.pdfType of
        // Continuous types
        //-----------------
        PdfBeta:
          newobj := TPdfBeta.create( (obj as TPdfBeta) )
        ;
        PdfBetaPERT:
          newobj := TPdfBetaPERT.create( (obj as TPdfBetaPERT) )
        ;
        PdfExponential:
          newobj := TPdfExponential.create( (obj as TPdfExponential) )
        ;
        PdfGamma:
          newobj := TPdfGamma.create( (obj as TPdfGamma) )
        ;
        PdfGaussian:
          newobj := TPdfGaussian.create( (obj as TPdfGaussian) )
        ;
        PdfHistogram:
          newObj := TPdfHistogram.create( (obj as TPdfHistogram) )
        ;
        PdfInverseGaussian:
          newObj := TPdfInverseGaussian.create( (obj as TPdfInverseGaussian) )
        ;
        PdfLogistic:
          newobj := TPdfLogistic.create( (obj as TPdfLogistic) )
        ;
        PdfLogLogistic:
          newobj := TPdfLoglogistic.create( (obj as TPdfLoglogistic) )
        ;
        PdfLognormal:
          newobj := TPdfLognormal.create( (obj as TPdfLognormal) )
        ;
        PdfPareto:
          newobj := TPdfPareto.create( (obj as TPdfPareto) )
        ;
        PdfPearson5:
          newobj := TPdfPearson5.create( (obj as TPdfPearson5) )
        ;
        PdfPiecewise:
          newobj := TPdfPiecewise.create( (obj as TPdfPiecewise) )
        ;
        PdfPoint:
          newobj := TPdfPoint.create( (obj as TPdfPoint) )
        ;
        PdfTriangular:
          newobj :=  TPdfTriangular.create( (obj as TPdfTriangular) )
        ;

        PdfUniform:
          newobj := TPdfUniform.create( (obj as TPdfUniform) )
        ;
        PdfWeibull:
          newobj := TPdfWeibull.create( (obj as TPdfWeibull) )
        ;

        // Discrete types
        //---------------
        PdfBernoulli:
          newObj := TPdfBernoulli.create( obj as TPdfBernoulli )
        ;
        PdfBinomial:
          newObj := TPdfBinomial.create( (obj as TPdfBinomial) )
        ;
        PdfDiscreteUniform:
          newObj := TPdfDiscreteUniform.create( (obj as TPdfDiscreteUniform) )
        ;
        PdfHypergeometric:
          newObj := TPdfHypergeometric.create( (obj as TPdfHypergeometric) )
        ;
        PdfNegativeBinomial:
          newObj := TPdfNegativeBinomial.create( (obj as TPdfNegativeBinomial) )
        ;
        PdfPoisson:
          newObj := TPdfPoisson.create( (obj as TPdfPoisson) )
        ;

        else
          begin
            raise exception.create( 'Unrecognized PDF type in createPdfFromBaseObject' );
            newobj := nil;
          end
        ;
      end;

      newobj.name := obj.name;
      newobj.id := obj.id;
      newobj.dbField := obj.dbField;

      result := newobj;
    end
  ;

  
  function pdfType( descr: string ): TPdfType;
    begin
      descr := fixup( descr );

      // Continuous types
      //-----------------
      if( 'beta' = descr ) then result := PdfBeta
      else if( 'beta-pert' = descr ) then result := PdfBetaPERT
      else if( 'exponential'= descr ) then result := PdfExponential
      else if( 'gamma' = descr ) then result := PdfGamma
      else if( 'gaussian' = descr ) then result := PdfGaussian
      else if( 'histogram' = descr ) then result := PdfHistogram
      else if( ( 'inverse-gaussian' = descr ) or ( 'inversegaussian' = descr ) ) then result := PdfInverseGaussian
      else if( 'logistic' = descr ) then result := PdfLogistic
      else if( 'loglogistic' = descr ) then result := PdfLogLogistic
      else if( 'lognormal' = descr ) then result := PdfLognormal
      else if( 'pareto' = descr ) then result := PdfPareto
      else if( 'pearson5' = descr ) then result := PdfPearson5
      else if( 'piecewise' = descr ) then result := PdfPiecewise
      else if( 'point' = descr ) then result := PdfPoint
      else if( 'triangular'= descr ) then result := PdfTriangular
      else if( 'uniform' = descr ) then result := PdfUniform
      else if( 'weibull' = descr ) then result := PdfWeibull

      // Discrete types
      //---------------
      else if( 'bernoulli' = descr ) then result := PdfBernoulli
      else if( 'binomial' = descr ) then result := PdfBinomial
      else if( 'discrete-uniform' = descr ) then result := PdfDiscreteUniform
      else if( 'hypergeometric' = descr ) then result := PdfHyperGeometric
      else if( 'negative-binomial' = descr ) then result := PdfNegativeBinomial
      else if( 'poisson' = descr ) then result := PdfPoisson

      else
        begin
          raise exception.create( 'Unrecognized PDF type (' + descr + ') in function pdfType' );
          result := PdfUndefined;
        end
      ;
    end
  ;


  function pdfTypeDescr( val: TPdfType ): string;
    begin
      case val of
        PdfUndefined: result := tr( '(Unspecified)' );

        // Continuous types
        //-----------------
        PdfBeta: result := tr( 'Beta' );
        PdfBetaPERT: result := tr( 'BetaPERT' );
        PdfExponential: result := tr( 'Exponential' );
        PdfGamma: result := tr( 'Gamma' );
        PdfGaussian: result := tr( 'Gaussian (normal)' );
        PdfHistogram: result := tr( 'Histogram' );
        PdfInverseGaussian: result := tr( 'Inverse Gaussian' );
        PdfLogistic: result := tr( 'Logistic' );
        PdfLogLogistic: result := tr( 'Loglogistic' );
        PdfLognormal: result := tr( 'Lognormal' );
        PdfPareto: result := tr( 'Pareto' );
        PdfPearson5: result := tr( 'Pearson 5' );
        PdfPiecewise: result := tr( 'Piecewise (general)' );
        PdfPoint: result := tr( 'Fixed value' );
        PdfTriangular: result := tr( 'Triangular' );
        PdfUniform: result := tr( 'Uniform' );
        PdfWeibull: result := tr( 'Weibull' );

        // Discrete types
        //---------------
        PdfBernoulli: result := tr( 'Bernoulli' );
        PdfBinomial: result := tr( 'Binomial' );
        PdfDiscreteUniform: result := tr( 'Discrete uniform' );
        PdfHypergeometric: result := tr( 'Hypergeometric' );
        PdfNegativeBinomial: result := tr( 'Negative binomial' );
        PdfPoisson: result := tr( 'Poisson' );

      else
        raise exception.Create( 'Unrecognized PdfType in pdfTypeDescr' );
      end;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Base class for all probability density functions
//-----------------------------------------------------------------------------
  constructor TPdf.create();
    begin
      inherited create();

      setChartType( CTPdf );
      yUnits := UUnitless;
      xUnits := UUnknown;

      _discreteVals := nil;
      _posDiscreteVals := nil;

      initialize();
    end
  ;


  constructor TPdf.create( const src: TPdf );
    begin
      dbcout( '*** Calling TPdf copy constructor.', DBSHOWMSG );
      inherited create( src );
      initialize();

      _chartType := CTPdf;

      _pdfType := src._pdfType;
      _yUnits := src._yUnits;
      _xUnits := src._xUnits;

      _infLeftTail := src._infLeftTail;
      _infRightTail := src._infRightTail;

      if( nil <> src._discreteVals ) then
        _discreteVals := TQIntegerDoubleMap.create( src._discreteVals )
      else
        _discreteVals := nil
      ;

      if( nil <> src._posDiscreteVals ) then
        _posDiscreteVals := TQIntegerDoubleMap.create( src._posDiscreteVals )
      else
        _posDiscreteVals := nil
      ;
    end
  ;


  destructor TPdf.destroy();
    begin
      //dbcout( 'Destroying TProbDensFunction ' + self.name, DBSHOWMSG );

      freeDiscreteValues();

      inherited destroy();
    end
  ;


  procedure TPdf.freeDiscreteValues();
    begin
      freeAndNil( _discreteVals );
      freeAndNil( _posDiscreteVals );
    end
  ;


  function TPdf.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      unitElement: pointer;
      nameAttribute: pointer;
    begin
      result := true; // This won't matter.
      
      // Each of the derived classes will have to handle its own parameters.
      //--------------------------------------------------------------------
      initialize();

      // Get the function name
      //----------------------
      nameAttribute := sdew.getAttributeByName( element, 'name' );
      if( nil <> nameAttribute ) then
        name := sdew.GetElementAttribute( element, 'name' )
      ;

      // Take care of the units
      //-----------------------
      //dbcout2( sdew.GetElementName( element ) );
      unitElement := Sdew.GetElementByName( element, 'units' );
      if ( nil <> unitElement ) then
        begin
          if ( nil <> Sdew.GetElementByName( unitElement, 'xdf:unit' ) ) then
            xUnits := chartUnitTypeFromXml( Sdew.GetElementContents( Sdew.GetElementByName( unitElement, 'xdf:unit' ) ) )
          else if ( nil <> Sdew.GetElementByName( unitElement, 'xdf:unitless' ) ) then
            xUnits := UUnitless
          ;
        end
      else
        xUnits := UUnknown
      ;
    end
  ;

  
  function TPdf.validate( err: pstring = nil ): boolean;
    begin
      if( err <> nil ) then err^ := tr( 'Function type is unspecified.' );
      result := false;
    end
  ;

  
  function TPdf.compare( Chart: TChartFunction ):boolean;
    begin
      result := inherited compare( Chart );

      if( result ) then
        result := ( self.pdfType = (Chart as TPdf).pdfType )
      ;
    end
  ;


  function TPdf.randNonNeg(): double;
    var
      r: double;
    begin
      r := -1.0;

      while( 0 > r ) do
        r := rand()
      ;

      result := r;
    end
  ;


  function TPdf.randInt(): integer;
    begin
      result := roundDbl( rand() );
    end
  ;


  function TPdf.randNonNegInt(): integer;
    var
      r: double;
    begin
      r := -1.0;

      // Anything between -0.5 and 0.5 will be rounded to 0.
      // Also watch for numbers too big to fit in an int.
      // (These should be extremely unlikely, but it doesn't
      // hurt to be prepared.)
      while( ( -0.5 > r ) or ( 2147483647 < r ) ) do
        r := rand()
      ;

      result := roundDbl( r );
    end
  ;


  function TPdf.randPosInt(): integer;
    var
      r: double;
    begin
      r := 0.0;

      // Anything between 0.5 and 1 will be rounded to 1.
      // Also watch for numbers too big to fit in an int.
      // (These should be extremely unlikely, but it doesn't
      // hurt to be prepared.)
      while( ( 0.5 > r ) or ( 2147483647 < r ) ) do
        r := rand()
      ;

      result := roundDbl( r );
    end
  ;


  procedure TPdf.randTest( const nSamples: integer = 10 );
    var
      i: integer;
      msg: string;
    begin
      dbcout( self.descr + ':', true );

      if( not( validate( @msg ) ) ) then
        begin
          dbcout( tr( 'Function is not valid:' ), true );
          dbcout( msg, true );
          dbcout( '', true );
        end
      else
        begin
          for i := 0 to nSamples - 1 do
            dbcout( self.rand(), true )
          ;
          dbcout( '', true );
        end
      ;
    end
  ;


  function TPdf.getMinD(): integer;
    begin
      if( self.isContinuous ) then
        begin
          if( self.hasMin ) then
            result := floor( self.min )
          else
            result := floor( self.invCumulative( 0.00001 ) )
          ;
        end
      else
        begin
          if( self.hasMin ) then
            result := floor( self.min )
          else
            result := (self as TPdfDiscrete).invCumulative( 0.00001 )
          ;
        end
      ;
    end
  ;


  function TPdf.getMaxD(): integer;
    begin
      if( self.isContinuous ) then
        begin
          if( self.hasMax ) then
            result := ceil( self.max )
          else
            result := ceil( self.invCumulative( 0.99999 ) )
          ;
        end
      else
        begin
          if( self.hasMax ) then
            result := ceil( self.max )
          else
            result := (self as TPdfDiscrete).invCumulative( 0.99999 )
          ;
        end
      ;
    end
  ;


  procedure TPdf.standardizePosDiscreteVals();
    var
      i: integer;
      total: double;
      min: integer;
    begin
      _posDiscreteVals := TQIntegerDoubleMap.create();

      total := 0.0;

      if( 0 <= minD ) then
        begin
          for i := minD to maxD do
            begin
              if( _discreteVals.contains( i ) ) then
                begin
                  total := total + _discreteVals.value( i );
                  _posDiscreteVals.insert( i, _discreteVals.value( i ) );
                end
              ;
            end
          ;

        end
      else
        begin
          min := math.Max( 0, minD );

          for i := min to maxD do
            begin
              if( _discreteVals.contains( i ) ) then
                total := total + _discreteVals.value( i )
              ;
            end
          ;

          for i := min to maxD do
            begin
              if( _discreteVals.contains( i ) ) then
                _posDiscreteVals.insert( i, ( _discreteVals.value( i ) / total ) )
              ;
            end
          ;
        end
      ;

      // Distributions without an upper bound won't reach a total probability of 1.
      // This ensures that, for our purposes, every distribution does.
      if( 1.0 > total ) then
        _posDiscreteVals.Add( maxD + 1, 1.0 - total )
      ;

      //_posDiscreteVals.debug();
    end
  ;

  
  function TPdf.getDiscreteVals(): TQIntegerDoubleMap;
    begin
      if( nil = _discreteVals ) then
        initializeDiscreteVals()
      ;
      result := _discreteVals;
    end
  ;


  function TPdf.getPosDiscreteVals(): TQIntegerDoubleMap;
    begin
      if( nil = _posDiscreteVals ) then
        initializeDiscreteVals()
      ;
      result := _posDiscreteVals;
    end
  ;


  function TPdf.discreteProb( x: integer ): double;
    begin
      if( nil = _discreteVals ) then
        initializeDiscreteVals()
      ;

      if( x < minD ) then
        result := 0.0
      else if( x > maxD ) then
        result := 0.0
      else
        result := _discreteVals.value( x )
      ;
    end
  ;


  function TPdf.discreteCumulativeProb( x: integer ): double;
    var
      max: integer;
      i: integer;
    begin
      if( nil = _discreteVals ) then
        initializeDiscreteVals()
      ;

      if( x < minD ) then
        result := 0.0
      else if( x > maxD ) then
        result := 1.0
      else
        begin
          max := math.Min( x, maxD );
          result := 0.0;
          for i := minD to max do
            result := result + _discreteVals.value( i )
          ;
        end
      ;
    end
  ;


  function TPdf.posDiscreteProb( x: integer ): double;
    begin
      if( nil = _discreteVals ) then
        initializeDiscreteVals()
      ;

      if( x < 0 ) then
        result := 0.0
      else if( x > maxD ) then
        result := 0.0
      else if( _posDiscreteVals.contains( x ) ) then
        result := _posDiscreteVals.value( x )
      else
        result := 0.0
      ;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdf.createDBFieldList( db: TSqlDatabase ): TQueryDictionary;
    var
      dict: TQueryDictionary;
    begin
      dict := TQueryDictionary.create();

      dict['fieldName']   := db.sqlQuote( smChartStr( TSMChart( self.dbField ) ) );
      dict['chartName']   := db.sqlQuote( self.name );
      dict['isPdf']       := usBoolToText( true );
      dict['chartType']   := '"Unspecified"';
      dict['mean']        := DATABASE_NULL_VALUE;
      dict['stddev']      := DATABASE_NULL_VALUE;
      dict['min']         := DATABASE_NULL_VALUE;
      dict['max']         := DATABASE_NULL_VALUE;
      dict['mode']        := DATABASE_NULL_VALUE;
      dict['location']    := DATABASE_NULL_VALUE;
      dict['scale']       := DATABASE_NULL_VALUE;
      dict['shape']       := DATABASE_NULL_VALUE;
      dict['alpha']       := DATABASE_NULL_VALUE;
      dict['alpha2']      := DATABASE_NULL_VALUE;
      dict['beta']        := DATABASE_NULL_VALUE;
      dict['n']           := DATABASE_NULL_VALUE;
      dict['p']           := DATABASE_NULL_VALUE;
      dict['m']           := DATABASE_NULL_VALUE;
      dict['d']           := DATABASE_NULL_VALUE;
      dict['dMin']        := DATABASE_NULL_VALUE;
      dict['dMax']        := DATABASE_NULL_VALUE;
      dict['theta']       := DATABASE_NULL_VALUE;
      dict['a']           := DATABASE_NULL_VALUE;
      dict['s']           := DATABASE_NULL_VALUE;
      dict['xAxisUnits']  := db.sqlQuote( chartUnitTypeAsSqlString( xUnits ) );
      dict['notes']       := db.sqlQuote( notes );

      result := dict;
    end
  ;
  {$ENDIF}

  function TPdf.getPdfType() : TPdfType; begin Result := _pdfType; end;
  procedure TPdf.setPdfType( val: TPdfType ); begin _pdfType := val; end;

  function TPdf.getYStartsAtZero(): boolean; begin result := true; end;

  function TPdf.getIsContinuous(): boolean; begin result := not( isDiscrete ); end;

  function TPdf.getHasMin(): boolean; begin result := not( _infLeftTail ); end;
  function TPdf.getHasMax(): boolean; begin result := not( _infRightTail ); end;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Base class for all continuous PDFs
//-----------------------------------------------------------------------------
  procedure fillLeftTail(
        arr: T2DPointList;
        tailPoints: integer;
        startPosition: double;
        desiredIntervalSize: double;
        maxY: double;
        recurseCounter: integer;
        fn: TPdf
      );
    var
      cumul: double;
      i: integer;
      stopNow: boolean;
      stopYLow, stopYHigh: double;
      x, y: double;
    begin
      stopNow := false;
      cumul := -1;
      stopYLow := maxY / 500;
      stopYHigh := maxY * 10;
      y := stopYLow;

      dbcout( '---LEFT TAIL', DBSHOWMSG );
      dbcout( 'Start position: ' + usFloatToStr( startPosition ), DBSHOWMSG );
      dbcout( 'desired interval size: ' + usFloatToStr( desiredIntervalSize ), DBSHOWMSG );
      dbcout( 'Tail points: ' + intToStr( tailPoints ), DBSHOWMSG );

      // Work backwards, to allow points to be prepended until 0 is reached
      for i := tailPoints - 1 downto 1 do
        begin
          cumul := startPosition - ( ( tailPoints - i ) * desiredIntervalSize );

          if
            ( fn.infiniteLeftTail and ( 0 = cumul ) )
          or
            ( 0 > cumul )
          then
            begin
              dbcout( '*** cumul has reached ' + usFloatToStr( cumul ), DBSHOWMSG );
              stopNow := true;
              break;
            end
          ;

          x := fn.invCumulative( cumul );
          y := fn.probDensity( x );

          if( fn.hasMin and ( x < fn.min ) ) then
            begin
              arr.append( T2DPoint.Create( fn.min, 0.0 ) );
              dbcout( '*** min has been reached.', DBSHOWMSG );
              stopNow := true;
              break;
            end
          else
            arr.prepend( T2DPoint.Create( x, y ) )
          ;

          dbcout( usFloatToStr( x ) + ', ' + usFloatToStr( y ), DBSHOWMSG );
        end
      ;

      if( y < stopYLow ) then
        begin
          dbcout( '*** min Y reached', DBSHOWMSG );
          stopNow := true;
        end
      ;

      if( y > stopYHigh ) then
        begin
          dbcout( '*** max Y reached', DBSHOWMSG );
          stopNow := true;
        end
      ;

      if( 25 = recurseCounter ) then
        begin
          dbcout( '*** Max recursion depth reached', DBSHOWMSG );
          stopNow := true;
        end
      ;

      if( stopNow ) then
        dbcout( 'Recursion stopped', DBSHOWMSG )
      else
        begin
          dbcout( 'Making recursive call', DBSHOWMSG );
          fillLeftTail( arr, tailPoints, cumul, cumul/5, maxY, recurseCounter + 1, fn );
        end
      ;

    end
  ;


  procedure fillRightTail(
        arr: T2DPointList;
        tailPoints: integer;
        startPosition: double;
        desiredIntervalSize: double;
        maxY: double;
        recurseCounter: integer;
        fn: TPdf
      );
    var
      cumul: double;
      i: integer;
      stopNow: boolean;
      stopYLow, stopYHigh: double;
      x, y: double;
    begin
      stopNow := false;
      cumul := -1;
      stopYLow := maxY / 500;
      stopYHigh := maxY * 10;
      y := stopYLow;

      dbcout( '--- RIGHT TAIL', DBSHOWMSG );
      dbcout( 'Start position: ' + usFloatToStr( startPosition ), DBSHOWMSG );
      dbcout( 'desired interval size: ' + usFloatToStr( desiredIntervalSize ), DBSHOWMSG );
      dbcout( 'Tail points: ' + intToStr( tailPoints ), DBSHOWMSG );

      for i := 1 to tailPoints - 1 do
        begin
          cumul := startPosition + ( i * desiredIntervalSize );

          if
            ( fn.infiniteRightTail and ( 1 = cumul ) )
          or
            ( 1 < cumul )
          then
            begin
              dbcout( '*** cumul has reached ' + usFloatToStr( cumul ), DBSHOWMSG );
              stopNow := true;
              break;
            end
          ;


          x := fn.invCumulative( cumul );
          y := fn.probDensity( x );

          if( fn.hasMax and ( x > fn.max ) ) then
            begin
              arr.append( T2DPoint.Create( fn.max, 0.0 ) );
              dbcout( '*** max has been reached.', DBSHOWMSG );
              stopNow := true;
              break;
            end
          else
            arr.append( T2DPoint.Create( x, y ) )
          ;

          dbcout( usFloatToStr( cumul ) + ', ' + usFloatToStr( x ) + ', ' + usFloatToStr( y ), DBSHOWMSG );
        end
      ;

      if( y < stopYLow ) then
        begin
          dbcout( '*** min Y reached', DBSHOWMSG );
          stopNow := true;
        end
      ;

      if( y > stopYHigh ) then
        begin
          dbcout( '*** max Y reached', DBSHOWMSG );
          stopNow := true;
        end
      ;

      if( 25 = recurseCounter ) then
        begin
          dbcout( '*** Max recursion depth reached', DBSHOWMSG );
          stopNow := true;
        end
      ;


      if( stopNow ) then
        dbcout( 'Recursion stopped', DBSHOWMSG )
      else
        begin
          dbcout( 'Making recursive call', DBSHOWMSG );
          fillRightTail( arr, tailPoints, cumul, desiredIntervalSize/5, maxY, recurseCounter + 1, fn );
        end
      ;

    end
  ;


  procedure trimLeftTail( arr: T2DPointList );
    var
      max: double;
    begin
      max := arr.maxY;

      dbcout( 'Checking for left trim:', DBSHOWMSG );
      dbcout( usFloatToStr( max ) + ', ' + usFloatToStr( arr.at(0).y ), DBSHOWMSG );

      if
        ( arr.at(0).y = max )
      and
        ( 10 < max - arr.minY )
      then
        begin
          dbcout( '*** trimming left tail', DBSHOWMSG );
          arr.removeAt(0);
          trimLeftTail( arr );
        end
      ;
    end
  ;


  procedure trimRightTail( arr: T2DPointList );
    var
      max: double;
    begin
      max := arr.maxY;
      if
        ( arr.at(arr.Count-1).y = max )
      and
        ( 10 < max - arr.minY )
      then
        begin
          dbcout( '*** trimming right tail', DBSHOWMSG );
          arr.removeAt(arr.Count-1);
          trimRightTail( arr );
        end
      ;
    end
  ;


  function TPdfContinuous.createPointArray(): T2DPointList;
    var
      arr: T2DPointList;

      i: integer;
      x, y: double;
      cumul: double;

      maxY, minY: double;
      xInterval, avgXInterval: double;
      desiredIntervalSize: double;

      tailPoints: integer;
    begin
      arr := T2DPointList.Create();

      // Create the bulk of the points first,
      // from p = 0.04 to p = 0.96
      //-------------------------------------
      dbcout( endl + '***', DBSHOWMSG );
      dbcout( descr, DBSHOWMSG );
      dbcout( 'i, cumul, x, y', DBSHOWMSG );
      for i := 1 to 24 do
        begin
          cumul := 0.04 * (i);
          x := invCumulative( cumul );
          y := probDensity( x );
          dbcout(
            intToStr( i )
              + ', ' + usFloatToStr( cumul )
              + ', ' + usFloatToStr( x )
              + ', ' + usFloatToStr( y ),
            DBSHOWMSG
          );

          arr.Append( T2DPoint.create( x, y ) );
        end
      ;

      
      // Determine some characteristics based on the bulk of the points
      //---------------------------------------------------------------
      maxY := arr.maxY();
      minY := arr.minY();
      avgXInterval := arr.avgXInterval();

      dbcout( '------------', DBSHOWMSG );
      dbcout( usFloatToStr( minY ) + ', ' + usFloatToStr( maxY ), DBSHOWMSG );
      dbcout( usFloatToStr( maxY - minY ), DBSHOWMSG );
      dbcout( usBoolToText( 10 < maxY - minY ), DBSHOWMSG );
      dbcout( '------------', DBSHOWMSG );


      // Deal with the left tail
      //------------------------

      // What is the smallest x value?
      if( infiniteLeftTail ) then
        cumul := 0.001
      else
        cumul := 0.0
      ;

      x := invCumulative( cumul );
      y := probDensity( x );

      // What is the interval from this smallest x to the next one that we've calculated?
      xInterval := arr.at(0).x - x;

      // If this interval is smaller than or equal to the average interval,
      // PREPEND the smallest x and call it a day.
      // Otherwise, divide the remaining tail area into pieces and PREPEND several points.
      if( xInterval <= avgXInterval ) then
        arr.prepend( T2DPoint.Create( x, y ) )
      else
        begin
          tailPoints := 5;
          desiredIntervalSize := 0.04 / 5;

          dbcout( 'MAXY: ' + usFloatToStr( maxY ), DBSHOWMSG );
          dbcout( 'STOPY: ' + usFloatToStr( maxY / 500 ), DBSHOWMSG );

          fillLeftTail( arr, tailPoints, 0.04, desiredIntervalSize, maxY, 0, self );
        end
      ;


      // Deal with the right tail
      //-------------------------

      // What is the largest X value?
      if( infiniteRightTail ) then
        cumul := 0.999
      else
        cumul := 1.0
      ;

      x := invCumulative( cumul );
      y := probDensity( x );

      // What is the interval from this largest x to the next one that we've calculated?
      xInterval := x - arr.at(arr.Count-1).x;

      // If this interval is larger than than or equal to the average interval,
      // APPEND the smallest x and call it a day.
      // Otherwise, divide the remaining tail area into pieces and APPEND several points.
      if( xInterval <= avgXInterval ) then
        arr.append( T2DPoint.Create( x, y ) )
      else
        begin
          tailPoints := 5;
          desiredIntervalSize := 0.04 / 5;

          dbcout( 'MAXY: ' + usFloatToStr( maxY ), DBSHOWMSG );
          dbcout( 'STOPY: ' + usFloatToStr( maxY / 500 ), DBSHOWMSG );

          fillRightTail( arr, tailPoints, 0.96, desiredIntervalSize, maxY, 0, self );
        end
      ;

      (*
      // AR 4/15/10: I don't think this tail trimming is necessary or desirable any more (if it ever really was).
      // Try it for a while without, and see what happens.
      
      // Check min and max again to see if tails should be trimmed.
      //-----------------------------------------------------------
      maxY := arr.maxY();
      minY := arr.minY();

      if( 10 < maxY - minY ) then
        begin
          dbcout( 'Should be trimming...', DBSHOWMSG );
          // Eliminate points from each end until this condition is no longer met
          trimRightTail( arr );
          trimLeftTail( arr );
        end
      ;
      *)

      //arr.debug();
      result := arr;
    end
  ;


  function TPdfContinuous.asCsv(): string;
    var
      i: integer;
      points: T2DPointList;
    begin
      result := '';

      if( self.validate() ) then
        begin
          points := createPointArray();

          // Set the first and last Y values 0.
          // This will make conversion to piecewise easier.
          points[0].y := 0.0;
          points[points.Count - 1].y := 0.0;

          try
            result := 'x' + csvListSep() + ' y' + endl;

            for i := 0 to points.Count - 1 do
              begin
                result :=
                  result
                  + csvFloatToStr( points[i].x )
                  + csvListSep() + ' '
                  + csvFloatToStr( points[i].y )
                  + endl
                ;
              end
            ;
          except
            raise exception.Create( 'An exception occurred in TPdfContinuous.asCsv()' );
            result := '';
          end;

          freeAndNil( points );
        end
      ;
    end
  ;


  procedure TPdfContinuous.initializeDiscreteVals();
    var
      i: integer;
      d: double;
      m: double;
      total: double;
    begin
      _discreteVals := TQIntegerDoubleMap.create();

      total := 0.0;

      for i := minD to maxD do
        begin
          m := math.max( self.min, i - 0.5 );
          d := self.cumulativeDistr( i + 0.5 ) - self.cumulativeDistr( m );
          _discreteVals.Add( i, d );
          total := total + d;
        end
      ;

      // Distributions without an upper bound won't reach a total probability of 1.
      // This ensures that, for our purposes, every distribution does.
      if( 1.0 > total ) then
        _discreteVals.Add( maxD + 1, 1.0 - total )
      ;

      //_discreteVals.debug();

      standardizePosDiscreteVals();
    end
  ;


  // FIX ME: write this function!
  function TPdfContinuous.getCanConvertToPiecewise(): boolean;
    begin
      result := false;
    end
  ;


  function TPdfContinuous.getIsDiscrete(): boolean;
    begin
      result := false;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Beta PDF
//-----------------------------------------------------------------------------
  constructor TPdfBeta.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;

  constructor TPdfBeta.create( alpha1, alpha2, min, max: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor
      
      setAlpha1( alpha1 );
      setAlpha2( alpha2 );
      setMin( min );
      setMax( max );
      xUnits := u;
    end
  ;


  constructor TPdfBeta.create( const src: TPdfBeta );
    begin
      dbcout( '*** Calling TPdfBeta copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setAlpha1( src.alpha1 );
      setAlpha2( src.alpha2 );
      setMin( src.min );
      setMax( src.max );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfBeta.initialize();
    begin
      setPdfType( PdfBeta );

      setAlpha1( NaN );
      setAlpha2( NaN );
      setMin( NaN );
      setMax( NaN );

      _infLeftTail := false;
      _infRightTail := false;
    end
  ;


  destructor TPdfBeta.destroy();
    begin
      //dbcout( 'Destroying TPdfBeta ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;
  
  function TPdfBeta.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;

      if( inherited compare( Chart ) ) then
        begin
          if ( self._alpha1 = (Chart as TPdfBeta)._alpha1 ) and
             ( self._alpha2 = (Chart as TPdfBeta)._alpha2 ) and
             ( self._min = (Chart as TPdfBeta)._min ) and
             ( self._max = (Chart as TPdfBeta)._max )
          then
            result := true
          ;
        end
      ;
    end;

  procedure TPdfBeta.setAlpha1( val: double ); begin _alpha1 := val; freeDiscreteValues(); end;
  procedure TPdfBeta.setAlpha2( val: double ); begin _alpha2 := val; freeDiscreteValues(); end;
  procedure TPdfBeta.setMin( val: double ); begin _min := val; freeDiscreteValues(); end;
  procedure TPdfBeta.setMax( val: double ); begin _max := val; freeDiscreteValues(); end;
  function TPdfBeta.getAlpha1(): double; begin Result := _alpha1; end;
  function TPdfBeta.getAlpha2(): double; begin Result := _alpha2; end;
  function TPdfBeta.getMin(): double; begin result := _min; end;
  function TPdfBeta.getMax(): double; begin result := _max; end;


  function TPdfBeta.getDescr(): string;
    begin
      result := tr( 'Beta' ) + format( ' ( %f, %f, %f, %f )', [alpha1, alpha2, min, max] );
    end
  ;


  function TPdfBeta.getMean(): double;
    begin
      result := min + ( max - min )*( alpha1 / ( alpha1 + alpha2 ) );
    end
  ;


  function TPdfBeta.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfBeta.validate( err: pstring = nil ): boolean;
    begin
      result := true;

      if( isNan(alpha1) or isNaN(alpha2) or isNaN(min) or isNan(max) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if ( 0.0 >= alpha1 ) or ( 0.0 >= alpha2 )
      then
        begin
          if( err <> nil ) then err^ := tr( 'Alpha1 and alpha2 must be greater than 0.' );
          result := false;
        end
      else if( max <= min ) then
        begin
          if( err <> nil ) then err^ := tr( 'Min must be less than max.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( ( 0.0 > min ) or ( 1.0 < max ) ) then
            begin
              if( err <> nil ) then err^ := tr( 'Proportions and probabilities must be between 0 and 1, inclusive.' );
              result := false;
            end
          ;
        end
      ;
    end
  ;


  function TPdfBeta.probDensity( x: double ): double;
    begin
      if( nil <> @aphi_beta_pdf ) then
        result := aphi_beta_pdf( x, self.alpha1, self.alpha2, self.min, self.max )
      else
        result := NaN
      ;
    end
  ;


  function TPdfBeta.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @aphi_beta_cdf ) then
        result := aphi_beta_cdf( x, self.alpha1, self.alpha2, self.min, self.max )
      else
        result := NaN
      ;
    end
  ;


  function TPdfBeta.invCumulative( p: double ): double;
    begin
      if( nil <> @aphi_beta_inverse_cdf ) then
        result := aphi_beta_inverse_cdf( p, self.alpha1, self.alpha2, self.min, self.max )
      else
        result := NaN
      ;
    end
  ;


  function TPdfBeta.rand(): double;
    begin
      if( nil <> @aphi_beta_rand ) then
        result := aphi_beta_rand( rngPtr(), self.alpha1, self.alpha2, self.min, self.max )
      else
        result := NaN
      ;
    end
  ;


  procedure TPdfBeta.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str
        + 'Beta PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 'alpha1=' + usFloatToStr(alpha1)
        + ', alpha2=' + usFloatToStr(alpha2)
        + ', min=' + usFloatToStr(min)
        + ', max=' + usFloatToStr(max)
      ;

      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfBeta.createCopy(): TChartFunction;
    begin
      result := TPdfBeta.create( self );
    end
  ;


  function TPdfBeta.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <beta>' + endl;
      result := result + _indent + '    <alpha>' + usFloatToStr( alpha1 ) + '</alpha>' + endl;
      result := result + _indent + '    <beta>' + usFloatToStr( alpha2 ) + '</beta>' + endl;
      result := result + _indent + '    <location>' + usFloatToStr( min ) + '</location>' + endl;
      result := result + _indent + '    <scale>' + usFloatToStr( max ) + '</scale>' + endl;
      result := result + _indent + '  </beta>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfBeta.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'beta' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setAlpha1( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'alpha' )), NaN ) );
          setAlpha2( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'beta' )), NaN ) );
          setMin( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'location' )), NaN ) );
          setMax( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'scale' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfBeta.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Beta"';
      dict['min']         := usFloatToStr( min );
      dict['max']         := usFloatToStr( max );
      dict['alpha']       := usFloatToStr( alpha1 );
      dict['alpha2']      := usFloatToStr( alpha2 );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// BetaPERT PDF
//-----------------------------------------------------------------------------
  constructor TPdfBetaPERT.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;

  constructor TPdfBetaPERT.create( min, mode, max: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setA( min );
      setC( mode );
      setB( max );
      xUnits := u;
    end
  ;


  constructor TPdfBetaPERT.create( const src: TPdfBetaPERT );
    begin
      dbcout( '*** Calling TPdfBetaPERT copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setA( src._a );
      setC( src._c );
      setB( src._b );
      xUnits := src._xUnits;
    end
  ;


  procedure TPdfBetaPERT.initialize();
    begin
      setPdfType( PdfBetaPERT );

      setA( NaN );
      setC( NaN );
      setB( NaN );

      _infLeftTail := false;
      _infRightTail := false;
    end
  ;


  destructor TPdfBetaPERT.destroy();
    begin
      //dbcout( 'Destroying TPdfTriangular ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfBetaPERT.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._a = (Chart as TpdfBetaPERT)._a ) and
             ( self._b = (Chart as TpdfBetaPERT)._b ) and
             ( self._c = (Chart as TpdfBetaPERT)._c ) then
            result := true;
        end;
    end;

  procedure TPdfBetaPERT.setA( val: double ); begin _a := val; freeDiscreteValues(); end;
  procedure TPdfBetaPERT.setC( val: double ); begin _c := val; freeDiscreteValues(); end;
  procedure TPdfBetaPERT.setB( val: double ); begin _b := val; freeDiscreteValues(); end;
  function TPdfBetaPERT.getA(): double; begin result := _a; end;
  function TPdfBetaPERT.getC(): double; begin result := _c; end;
  function TPdfBetaPERT.getB(): double; begin result := _b; end;
  function TPdfBetaPERT.getMin(): double; begin result := _a; end;
  function TPdfBetaPERT.getMax(): double; begin result := _b; end;

  function TPdfBetaPERT.getDescr(): string;
    begin
      result := tr( 'BetaPERT' ) + format( ' ( %f, %f, %f )', [_a, _c, _b] );
    end
  ;



  function TPdfBetaPERT.getMean(): double;
    begin
      result := ( min + 4 * mode + max ) / 6;
    end
  ;


  function TPdfBetaPERT.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfBetaPERT.validate( err: pstring = nil ): boolean;
    begin
      result := true;
      
      if( isNaN( _a ) or isNaN( _b ) or isNan( _c ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( _a < 0 ) or ( _b < 0 ) or ( _c < 0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Values must be greater than or equal to 0.' );
          result := false;
        end
      else if not( (_a < _c) and (_c < _b) ) then
        begin
          if( err <> nil ) then err^ := tr( 'The minimum value must be less than the mode, which must be less than the maximum.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( ( 0.0 > min ) or ( 1.0 < max ) ) then
            begin
              if( err <> nil ) then err^ := tr( 'Proportions and probabilities must be between 0 and 1, inclusive.' );
              result := false;
            end
          ;
        end
      ;
    end
  ;


  function TPdfBetaPERT.probDensity( x: double ): double;
    begin
      if( nil <> @aphi_beta_pert_pdf ) then
        result := aphi_beta_pert_pdf( x, self.min, self.mode, self.max )
      else
        result := NaN
      ;
    end
  ;

  function TPdfBetaPERT.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @aphi_beta_pert_cdf ) then
        result :=aphi_beta_pert_cdf( x, self.min, self.mode, self.max )
      else
        result := NaN
      ;
    end
  ;

  function TPdfBetaPERT.invCumulative( p: double ): double;
    begin
      if( nil <> @aphi_beta_pert_inverse_cdf ) then
        result := aphi_beta_pert_inverse_cdf( p, self.min, self.mode, self.max )
      else
        result := NaN
      ;
    end
  ;


  function TPdfBetaPERT.rand(): double;
    begin
      if( nil <> @aphi_beta_pert_rand ) then
        result := aphi_beta_pert_rand( rngPtr(), self.min, self.mode, self.max )
      else
        result := NaN
      ;
    end
  ;


  procedure TPdfBetaPERT.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'BetaPERT PDF (id ' + intToStr( self.id ) + ') :' + endl + 'min=' + usFloatToStr(_a) + ', max=' + usFloatToStr(_b) + ', mode=' + usFloatToStr(_c);
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfBetaPERT.createCopy(): TChartFunction;
    begin
      result := TPdfBetaPERT.create( self );
    end
  ;


  function TPdfBetaPERT.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <beta-pert>' + endl;
      // FIX ME: Is this comment still correct?
      result := result + _indent + '    <!-- NAADSM uses a set of approximations (from Davis) for the BetaPERT PDF -->' + endl;
      result := result + _indent + '    <!-- which are slightly different from those described by Vose and Palisade. -->' + endl;
      result := result + _indent + '    <!-- 1) Davis, R.E. http://www.cob.sjsu.edu/facstaff/davis_r/courses/QBAreader/QBAtoc.html -->' + endl;
      result := result + _indent + '    <!-- 2) Vose, D. 1996.  Quantitative risk analysis: a guide to Monte Carlo Simulation Modelling. -->' + endl;
      result := result + _indent + '    <!--    John Wiley and Sons, New York. -->' + endl;
      result := result + _indent + '    <!-- 3) Palisade Corporation. 2002.  A concise summary of @RISK probability distribution functions. -->' + endl;
      result := result + _indent + '    <!--    Available (probably illegally) from http://project.zf.jcu.cz/risk/data/distfunc.pdf -->' + endl;

      result := result + _indent + '    <min>' + usFloatToStr( min ) + '</min>' + endl;
      result := result + _indent + '    <mode>' + usFloatToStr( mode ) + '</mode>' + endl;
      result := result + _indent + '    <max>' + usFloatToStr( max ) + '</max>' + endl;
      result := result + _indent + '  </beta-pert>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfBetaPERT.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'beta-pert' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setA( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'min' )), NaN ) );
          setB( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'max' )), NaN ) );
          setC( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'mode' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfBetaPERT.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Beta-Pert"';
      dict['min']         := usFloatToStr( min );
      dict['mode']        := usFloatToStr( mode );
      dict['max']         := usFloatToStr( max );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Exponential PDF
//-----------------------------------------------------------------------------
  constructor TPdfExponential.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfExponential.create( mean: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setMean( mean );
      xUnits := u;
    end
  ;

  constructor TPdfExponential.create( const src: TPdfExponential );
    begin
      dbcout( '*** Calling TPdfExponential copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setMean( src.mean );
      xUnits := src._xUnits;
    end
  ;


  procedure TPdfExponential.initialize();
    begin
      setPdfType( PdfExponential );

      setMean( NaN );

      _infLeftTail := false;
      _infRightTail := true;
    end
  ;


  destructor TPdfExponential.destroy();
    begin
      //dbcout( 'Destroying TPdfExponential ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfExponential.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._mean = (Chart as TPdfExponential)._mean ) then
            result := true;
        end;
    end;

  function TPdfExponential.getMean(): double; begin result := _mean; end;
  procedure TPdfExponential.setMean( val: double ); begin _mean := val; freeDiscreteValues(); end;

  function TPdfExponential.getMin(): double;
    begin
      result := 0.0;
    end
  ;

  function TPdfExponential.getMax(): double;
    begin
      result := Infinity;
    end
  ;

  
  function TPdfExponential.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfExponential.validate( err: pstring = nil ): boolean;
    begin
      if( isNan( mean ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'Mean must be specified.' );
          result := false;
        end
      else if( mean <= 0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Mean must be greater than 0.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfExponential.probDensity( x: double ): double;
    begin
      if( nil <> @gslRanExponentialPdf ) then
        result := gslRanExponentialPdf( x, self.mean )
      else
        result := NaN
      ;
    end
  ;


  function TPdfExponential.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @gslCdfExponentialP ) then
        result := gslCdfExponentialP( x, self.mean )
      else
        result := NaN
      ;
    end
  ;


  function TPdfExponential.invCumulative( p: double ): double;
    begin
      if( nil <> @gslCdfExponentialPinv ) then
        result := gslCdfExponentialPInv( p, self.mean )
      else
        result := NaN
      ;
    end
  ;


  function TPdfExponential.rand(): double;
    begin
      if( nil <> @gslRanExponential ) then
        result := gslRanExponential( gslRngPtr(), self.mean )
      else
        result := NaN
      ;
    end
  ;


  function TPdfExponential.getDescr(): string;
    begin
      result := tr( 'Exponential' ) + format( ' ( %f )', [_mean] );
    end
  ;


  procedure TPdfExponential.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Exponential PDF (id ' + intToStr( self.id ) + ') : mean=' + usFloatToStr(self.mean);
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfExponential.createCopy(): TChartFunction;
    begin
      result := TPdfExponential.create( self );
    end
  ;


  function TPdfExponential.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <exponential>' + endl;
      result := result + _indent + '    <mean>' + usFloatToStr( mean ) + '</mean>' + endl;
      result := result + _indent + '  </exponential>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfExponential.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'exponential' );


      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;
                  
          setMean( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'mean' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfExponential.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Exponential"';
      dict['mean']        := usFloatToStr( mean );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Gamma PDF
//-----------------------------------------------------------------------------
  constructor TPdfGamma.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;

  constructor TPdfGamma.create( alpha, beta: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setAlpha( alpha );
      setBeta( beta );
      xUnits := u;
    end
  ;


  constructor TPdfGamma.create( const src: TPdfGamma );
    begin
      dbcout( '*** Calling TPdfGamma copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setAlpha( src.alpha );
      setBeta( src.beta );
      xUnits := src._xUnits;
    end
  ;


  procedure TPdfGamma.initialize();
    begin
      setPdfType( PdfGamma );

      setAlpha( NaN );
      setBeta( NaN );

      _infLeftTail := false;
      _infRightTail := true;
    end
  ;


  destructor TPdfGamma.destroy();
    begin
      //dbcout( 'Destroying TPdfGamma ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfGamma.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._alpha = (Chart as TPdfGamma)._alpha ) and
             ( self._beta = (Chart as TPdfGamma)._beta ) then
            result := true;
        end;
    end;

  procedure TPdfGamma.setAlpha( val: double ); begin _alpha := val; freeDiscreteValues(); end;
  procedure TPdfGamma.setBeta( val: double ); begin _beta := val; freeDiscreteValues(); end;
  function TPdfGamma.getAlpha(): double; begin Result := _alpha; end;
  function TPdfGamma.getBeta(): double; begin Result := _beta; end;


  function TPdfGamma.validate( err: pstring = nil ): boolean;
    begin
      if( isNan( alpha ) or isNaN( beta ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( beta <= 0.1 ) then
        begin
          if( err <> nil ) then err^ := tr( 'The specified beta value is too small for practical calculation.' );
          result := false;
        end
      else if( alpha <= 0.1 ) then
        begin
          if( err <> nil ) then err^ := tr( 'The specified alpha value is too small for practical calculation.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfGamma.probDensity( x: double ): double;
    begin
      if( nil <> @gslRanGammaPdf ) then
        result := gslRanGammaPdf( x, self.alpha, self.beta )
      else
        result := NaN
      ;
    end
  ;


  function TPdfGamma.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @gslCdfGammaP ) then
        result := gslCdfGammaP( x, self.alpha, self.beta )
      else
        result := NaN
      ;
    end
  ;


  function TPdfGamma.invCumulative( p: double ): double;
    begin
      if( nil <> @gslCdfGammaPinv ) then
        result := gslCdfGammaPInv( p, self.alpha, self.beta )
      else
        result := NaN
      ;
    end
  ;


  function TPdfGamma.rand(): double;
    begin
      if( nil <> @gslRanGamma ) then
        result := gslRanGamma( gslRngPtr(), self.alpha, self.beta )
      else
        result := NaN
      ;
    end
  ;


  function TPdfGamma.getDescr(): string;
    begin
      result := tr( 'Gamma' ) + format( ' ( %f, %f )', [_alpha, _beta] );
    end
  ;


  function TPdfGamma.getMean(): double;
    begin
      result := alpha * beta;
    end
  ;


  function TPdfGamma.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfGamma.getMin(): double;
    begin
      result := 0.0;
    end
  ;

  function TPdfGamma.getMax(): double;
    begin
      result := infinity;
    end
  ;


  procedure TPdfGamma.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Gamma PDF (id ' + intToStr( self.id ) + ') :' + endl + 'alpha=' + usFloatToStr(alpha) + ', beta=' + usFloatToStr(beta);
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfGamma.createCopy(): TChartFunction;
    begin
      result := TPdfGamma.create( self );
    end
  ;


  function TPdfGamma.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <gamma>' + endl;
      result := result + _indent + '    <alpha>' + usFloatToStr( alpha ) + '</alpha>' + endl;
      result := result + _indent + '    <beta>' + usFloatToStr( beta ) + '</beta>' + endl;
      result := result + _indent + '  </gamma>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfGamma.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'gamma' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setAlpha( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'alpha' )), NaN ) );
          setBeta( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'beta' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfGamma.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Gamma"';
      dict['alpha']       := usFloatToStr( alpha );
      dict['beta']        := usFloatToStr( beta );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Gaussian (normal) PDF
//-----------------------------------------------------------------------------
  constructor TPdfGaussian.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfGaussian.create( mean, stddev: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setMean( mean );
      setStddev( stddev );
      xUnits := u;
    end
  ;


  constructor TPdfGaussian.create( const src: TPdfGaussian );
    begin
      dbcout( '*** Calling TPdfGaussian copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setMean( src.mean );
      setStddev( src.stddev );
      xUnits := src._xUnits;
    end
  ;


  procedure TPdfGaussian.initialize();
    begin
      setPdfType( PdfGaussian );

      setMean( NaN );
      setStddev( NaN );

      _infLeftTail := true;
      _infRightTail := true;
    end
  ;


  destructor TPdfGaussian.destroy();
    begin
      //dbcout( 'Destroying TPdfGaussian ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfGaussian.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._mean = (Chart as TPdfGaussian)._mean ) and
             ( self._stddev = (Chart as TPdfGaussian)._stddev ) then
            result := true;
        end;
    end;

  function TPdfGaussian.getMean(): double; begin Result := _mean; end;
  function TPdfGaussian.getStdDev(): double; begin Result := _stddev; end;
  procedure TPdfGaussian.setMean( val: double ); begin _mean := val; freeDiscreteValues(); end;
  procedure TPdfGaussian.setStdDev( val: double ); begin _stddev := val; freeDiscreteValues(); end;

   function TPdfGaussian.getMin(): double;
    begin
      result := negInfinity;
    end
  ;

  function TPdfGaussian.getMax(): double;
    begin
      result := infinity;
    end
  ;

  function TPdfGaussian.getHasMean(): boolean;
    begin
      result := true;
    end
  ;

      
  function TPdfGaussian.validate( err: pstring = nil ): boolean;
    begin
      if( isNaN( mean ) or isNan( stddev ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( mean < 0 ) or ( stddev < 0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Values must be greater than or equal to 0.' );
          result := false;
        end
      else if( stddev <= 0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Standard deviation must be greater than 0.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfGaussian.getDescr(): string;
    begin
      result := tr( 'Gaussian' ) + format( ' ( %f, %f )', [_mean, _stddev] );
    end
  ;


  function TPdfGaussian.probDensity( x: double ): double;
    begin
      if( nil <> @gslRanGaussianPdf ) then
        result := gslRanGaussianPdf( x - self.mean, self.stddev )
      else
        result := NaN
      ;
    end
  ;


  class function TPdfGaussian.probDensity( const x, mean, stdev: double ): double;
    begin
      if( nil <> @gslRanGaussianPdf ) then
        result := gslRanGaussianPdf( x - mean, stdev )
      else
        result := NaN
      ;
    end
  ;

  function TPdfGaussian.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @gslCdfGaussianP ) then
        result := gslCdfGaussianP( x - self.mean, self.stddev )
      else
        result := NaN
      ;
    end
  ;


  function TPdfGaussian.invCumulative( p: double ): double;
    begin
      if( nil <> @gslCdfGaussianPinv ) then
        result := gslCdfGaussianPInv( p, self.stddev ) + self.mean
      else
        result := NaN
      ;
    end
  ;


  function TPdfGaussian.rand(): double;
    begin
      if( nil <> @gslRanGaussian ) then
        result := gslRanGaussian( gslRngPtr(), self.stddev ) + self.mean
      else
        result := NaN
      ;
    end
  ;


  class function TPdfGaussian.rand( const mean, stdev: double ): double;
    begin
      if( nil <> @gslRanGaussian ) then
        result := gslRanGaussian( gslRngPtr(), stdev ) + mean
      else
        result := NaN
      ;
    end
  ;

  
  procedure TPdfGaussian.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Gaussian PDF (id ' + intToStr( self.id ) + ') : mean=' + usFloatToStr(self.mean) + ', ' + endl + 'stddev=' + usFloatToStr(stddev);
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfGaussian.createCopy(): TChartFunction;
    begin
      result := TPdfGaussian.create( self );
    end
  ;


  function TPdfGaussian.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <gaussian>' + endl;
      result := result + _indent + '    <mean>' + usFloatToStr( mean ) + '</mean>' + endl;
      result := result + _indent + '    <stddev>' + usFloatToStr( stddev ) + '</stddev>' + endl;
      result := result + _indent + '  </gaussian>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfGaussian.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'gaussian' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setMean( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'mean' )), NaN ) );
          setStdDev( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'stddev' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfGaussian.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Gaussian"';
      dict['mean']        := usFloatToStr( mean );
      dict['stddev']      := usFloatToStr( stddev );

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

      // KEEP THIS CODE: it might be useful when I write the chart conversion function.
      (*
      // FIX ME: This junk is only here until the interface is updated.
      // Fill the points in the database, according to Mark's procedure
      xMin, xMax: real;
      interval : real;
      i : integer;
      dbX, dbY: real;

      xMin := mean - (3*stddev);
      xMax := mean + (3*stddev);
      interval := (xmax - xmin) / 19;

      q := 'INSERT INTO [inChartDetail] ([chartID], [pointOrder], [x], [y])';


      for i := 0 to 19 do
        begin
          dbX := xMin + (interval*i);

          if i in [0,19] then
            dbY := 0
          else
            dbY := GetNormalFunction( dbX )
          ;

          q2 := q + ' VALUES (' + usFloatToStr( id ) + ', ' + intToStr(i) + ', ' + usFloatToStr( dbX ) + ', ' + usFloatToStr( dbY) + ')';
          db.execute( q2 );

        end
      ;
      *)
    end
  ;
  {$ENDIF}
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Histogram PDF
//-----------------------------------------------------------------------------
  // Must be used with setRangesAndBinValues to be useful.
  constructor TPdfHistogram.create();
    begin
      dbcout( '-- Begin TPdfHistogram.create():1' , DBSHOWMSG );
      inherited create();
      // initialize() is called by the inherited constructor

      dbcout( '--Done TPdfHistogram.create():1' , DBSHOWMSG );
    end
  ;


  constructor TPdfHistogram.create( ranges, counts: TARDoubleArray; u: TChartUnitType = UUnknown );
    begin
      dbcout( '-- Begin TPdfHistogram.create():2' , DBSHOWMSG );
      inherited create();
      // initialize() is called by the inherited constructor

      xUnits := u;

      setRangesAndCounts( ranges, counts );
      dbcout( '-- Done TPdfHistogram.create():2' , DBSHOWMSG );
    end
  ;


  constructor TPdfHistogram.create( rcd: RHistogramPointArray; u: TChartUnitType = UUnknown );
    var
      i: integer;
    begin
      dbcout( '-- Begin TPdfHistogram.create():fromHistogramPointArray' , DBSHOWMSG );
      inherited create();
      // initialize() is called by the inherited constructor

      xUnits := u;

      // There must be at least 2 points (upper and lower range and one bin value)
      // to define a histogram distribution.
      if( 1 < length( rcd ) ) then
        begin
          setLength( _ranges, length( rcd ) );
          setLength( _densities, length( rcd ) - 1 );
          setLength( _counts, length( rcd ) - 1 );

          for i := 0 to length( rcd ) - 1 do
            begin
              _ranges[i] := rcd[i].range;

              if( i < ( length( rcd ) - 1 ) ) then
                _counts[i] := rcd[i].count
              ;
            end
          ;

          processRangesAndDensities();
        end
      ;
      dbcout( '-- Done TPdfHistogram.create():fromHistogramPointArray' , DBSHOWMSG );
    end
  ;


  constructor TPdfHistogram.create( const src: TPdfHistogram );
    begin
      dbcout( '-- Begin TPdfHistogram.create():CopyConstructor' , DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      xUnits := src.xUnits;

      setRangesAndCounts( src._ranges, src._counts );

      dbcout( '-- Done TPdfHistogram.create():CopyConstructor' , DBSHOWMSG );
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  constructor TPdfHistogram.create( db: TSqlDatabase; chartID: integer; u: TChartUnitType = UUnknown );
    var
      q: string;
      res: TSqlResult;
      row: TSqlRow;
      i: integer;
      nRows: integer;
    begin
      dbcout( '-- Begin TPdfHistogram.create():Database' , DBSHOWMSG );

      inherited create();
      // initialize() is called by the inherited constructor

      xUnits := u;

      q := 'SELECT [x], [y] FROM [inChartDetail] WHERE [chartID] = ' + intToStr( chartID ) + ' ORDER BY [pointOrder]';

      res := TSqlResult.create( q, db );

      // There must be at least 2 points (upper and lower range and one bin value)
      // to define a histogram distribution.
      nRows := res.numRows;
      if( 1 < nRows ) then
        begin
          setLength( _ranges, res.numRows );
          setLength( _counts, res.numRows - 1 );
          setLength( _densities, res.numRows - 1 );

          row := res.fetchArrayFirst();
          i := 0;

          while( row <> nil ) do
            begin
              _ranges[i] := double( row.field(0) );

              if( i < ( nRows - 1 ) ) then
                _counts[i] := double( row.field(1) )
              ;

              inc( i );
              row := res.fetchArrayNext();
            end
          ;

          processRangesAndDensities();
        end
    ;

      freeAndNil( res );

      dbcout( '-- Done TPdfHistogram.create():Database' , DBSHOWMSG );
    end
  ;
  {$ENDIF}


  procedure TPdfHistogram.initialize();
    begin
      dbcout( '-- Begin TPdfHistogram.initialize' , DBSHOWMSG );

      setPdfType( PdfHistogram );

      setLength( _ranges, 0 );
      setLength( _counts, 0 );
      setLength( _densities, 0 );

      _infLeftTail := false;
      _infRightTail := false;

      _cStructPtr := nil;

      dbcout( '-- Done TPdfHistogram.initialize' , DBSHOWMSG );
    end
  ;


  procedure TPdfHistogram.setRangesAndCounts( ranges, counts: TARDoubleArray );
    begin
      dbcout( '-- Begin TPdfHistogram.setRangesAndCounts' , DBSHOWMSG );

      freeDiscreteValues();

      setLength( _ranges, length( ranges ) );
      setLength( _counts, length( counts ) );
      setLength( _densities, length( counts ) );

      _ranges := copy( ranges, 0, length( ranges ) );
      _counts := copy( counts, 0, length( counts ) );

      processRangesAndDensities();

      dbcout( '-- Done TPdfHistogram.setRangesAndCounts' , DBSHOWMSG );
    end
  ;


  procedure TPdfHistogram.processRangesAndDensities();
    var
      s: string;
    begin
      dbcout( '-- Begin TPdfHistogram.processRangesAndDensities' , DBSHOWMSG );

      if( libAphiLoaded and ( nil <> _cStructPtr ) ) then
        aphi_free_histogram( _cStructPtr )
      ;
      _cStructPtr := nil;

      if( validate ) then
        begin
          standardizeArea();

          if( libAphiLoaded and validate( @s ) ) then
            _cStructPtr := aphi_create_histogram( length( _ranges ), _ranges, _counts )
          else
            dbcout( s, true )
          ;
        end
      ;

      dbcout( '-- Done TPdfHistogram.processRangesAndDensities' , DBSHOWMSG );
    end
  ;


  destructor TPdfHistogram.destroy();
    begin
      if( nil <> _cStructPtr ) then
        begin
          aphi_free_histogram( _cStructPtr );
          _cStructPtr := nil;
        end
      ;

      //dbcout( 'Destroying TPdfHistogram ' + self.name, DBSHOWMSG );
      setLength( _ranges, 0 );
      setLength( _counts, 0 );
      setLength( _densities, 0 );
      
      inherited destroy();
    end
  ;

  function TPdfHistogram.compare( Chart:TChartFunction ):boolean;
    var
      i, count: integer;
      pdf: TPdfHistogram;
    begin
      result := false;
      
      if ( inherited compare( Chart ) ) then
        begin
          pdf := Chart as TPdfHistogram; 
        
          if ( ( self.xCount = pdf.xCount ) and ( self.yCount = pdf.yCount ) )  then
            begin
              count := length( self._ranges );

              result := true;
               
              for i := 0 to count - 1 do
                begin
                  if( self._ranges[i] <> pdf._ranges[i] ) then
                    begin
                      result := false;
                      break;
                    end
                  ;
                end
              ;

              if( true = result ) then
                begin
                  count := length( self._counts );
                  
                  for i := 0 to count - 1 do
                    begin
                      if( self._counts[i] <> pdf._counts[i] ) then
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
    end
  ;


  function TPdfHistogram.getMin(): double;
    begin
      result := _ranges[0];
    end
  ;


  function TPdfHistogram.getMax(): double;
    begin
      result := _ranges[ length( _ranges ) - 1 ];
    end
  ;
  

  function TPdfHistogram.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfHistogram.xyValsOK( msg: Pstring = nil ): boolean;
    var
      i: integer;
      prev: double;
      sum: double;
    begin
      // Each x value must be greater than the preceding one.
      // All y (probability) values must be greater than or equal to 0.
      // The sum of the Y values must be greater than 0 (otherwise, it's a flat line).
      
      result := true;
      
      // Verify that each x is greater than the preceding one
      //-----------------------------------------------------
      prev := _ranges[0];
      
      for i := 1 to high( _ranges ) do
        begin
          //dbcout2( _ranges[i] );
          if( _ranges[i] <= prev ) then
            begin
              result := false;
              if ( nil <> msg ) then msg^ := tr( 'Each X value must be greater than the preceding one.' );
              exit;
            end
          else
            prev := _ranges[i]
          ;
        end
      ;

      // Verify that all counts are greater than or equal to 0
      //-------------------------------------------------------
      sum := 0.0;
      
      for i := 0 to length( _counts ) - 1 do
        begin
          if( 0.0 > _counts[i] ) then
            begin
              result := false;
              if( nil <> msg ) then msg^ := tr( 'Every histogram bin count must be greater than or equal to 0.' );
              exit;
            end
          else
            sum := sum + _counts[i]
          ;
        end
      ;

      // Verify that the sum of the bin counts is greater than 0.
      //---------------------------------------------------------
      if( 0.0 >= sum ) then
        begin
          result := false;
          if( nil <> msg ) then msg^ := tr( 'The sum of the bin values of a histogram PDF must be greater than 0.' );
        end
      ; 
    end
  ;



  function TPdfHistogram.createPointArray(): T2DPointList;
    var
      arr: T2DPointList;
      i: integer;
    begin
      arr := T2DPointList.Create();

      arr.append( T2DPoint.create( _ranges[0], 0 ) );

      for i := 0 to length( _counts ) - 1 do
        begin
          arr.append( T2DPoint.create( _ranges[i], _densities[i] ) );
          arr.append( T2DPoint.create( _ranges[i+1], _densities[i] ) )
        end
      ;

      arr.append( T2DPoint.create( _ranges[ length( _ranges ) - 1 ], 0 ) );

      result := arr;
    end
  ;


  function TPdfHistogram.createHistogramPointArray(): RHistogramPointArray;
    var
      rcd: RHistogramPointArray;
      i: integer;
    begin
      setLength( rcd, length( _ranges ) );

      for i := low( _ranges ) to high( _ranges ) - 1 do
        begin
          rcd[i].range := _ranges[i];
          rcd[i].count := _counts[i];
          rcd[i].density := _densities[i];
        end
      ;

      rcd[ high( rcd ) ].range := _ranges[ high( _ranges ) ];
      rcd[ high( rcd ) ].count := 0.0;
      rcd[ high( rcd ) ].density := 0.0;

      result := rcd;
    end
  ;

      
  function TPdfHistogram.validate( err: pstring = nil ): boolean;
    begin
      // xCount must have a minimum of 2
      // xCount must be equal to yCount + 1
      // Each x value must be greater than the preceding one.
      // All y (probability) values must be greater than or equal to 0.
      // The sum of the Y values must be greater than 0 (otherwise, it's a flat line).
      // If a probability or proportion, min must be >= 0 and max must be <= 1

      result := true;

      if( 2 > xCount ) then
        begin
          if( err <> nil ) then err^ := tr( 'At least one histogram bin must be specified.' );
          result := false;
        end
      else if( yCount <> xCount - 1 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Histogram bins are not properly specified.' );
          result := false;
        end
      else if( not( xyValsOK( err ) ) ) then
        result := false  
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( ( 0.0 > min ) or ( 1.0 < max ) ) then
            begin
              if( err <> nil ) then err^ := tr( 'Proportions and probabilities must be between 0 and 1, inclusive.' );
              result := false;
            end
          ;
        end
      ;
    end
  ;


  procedure TPdfHistogram.standardizeArea();
    var
      i: integer;
      sum: double;
    begin
      dbcout( '--- Begin TPdfHistogram.standardizeArea' , DBSHOWMSG );

      sum := 0.0;

      for i := 0 to length( _counts ) - 1 do
        sum := sum + _counts[i]
      ;

      for i := 0 to length( _counts ) - 1 do
        _densities[i] := _counts[i] / sum / ( _ranges[i+1] - _ranges[i] )
      ;

      dbcout( '--- Done TPdfHistogram.standardizeArea' , DBSHOWMSG );
    end
  ;

  
  class function TPdfHistogram.standardize( var rcd: RHistogramPointArray; msg: PString = nil ): boolean;
    var
      i: integer;
      ranges: TARDoubleArray;
      counts: TARDoubleArray;
      sum: double;
    begin
      dbcout( '--- Begin TPdfHistogram.standardize (class function)' , DBSHOWMSG );

      result := true;
      
      // Check that there are enough points to work with
      //------------------------------------------------
      // There must be at least 2 points (upper and lower range and one bin value)
      // to define a histogram distribution.
      if( 2 > length( rcd ) ) then
        begin
          //raise exception.create( 'There are not enough points to define a histogram in TPdfHistogram.standardize()' );
          result := false;
          exit;
        end
      ;

      // Create ranges and bin vals from the points array
      //-------------------------------------------------
      setLength( ranges, length( rcd ) );
      setLength( counts, length( rcd ) - 1 );

      for i := 0 to length( rcd ) - 1 do
        begin
          ranges[i] := rcd[i].range;

          if( i < ( length( rcd ) - 1 ) ) then
            counts[i] := rcd[i].count
          ;
        end
      ;

      // Standardize the area under the curve to 1
      //------------------------------------------
      sum := 0.0;

      for i := 0 to length( counts ) - 1 do
        sum := sum + counts[i]
      ;

      for i := 0 to length( counts ) -1 do
        rcd[i].density := counts[i] / sum / ( ranges[i+1] - ranges[i] )
      ;

      // Clean up
      //---------
      setLength( ranges, 0 );
      setLength( counts, 0 );

      dbcout( '--- Done TPdfHistogram.standardize (class function)' , DBSHOWMSG );
    end
  ;

  
  function TPdfHistogram.getMean(): double;
    begin
      if( libAphiLoaded and ( nil <> _cStructPtr ) ) then
        result := aphi_histogram_mean( _cStructPtr )
      else
        begin
          //raise exception.create( 'DLL problem in TPdfHistogram.getMean()' );
          result := NaN;
        end
      ;
    end
  ;


  function TPdfHistogram.getDescr(): string;
    begin
      result := tr( 'Histogram' );
    end
  ;


  function TPdfHistogram.getXCount(): integer;
    begin
      result := length( _ranges );
    end
  ;


  function TPdfHistogram.getYCount(): integer;
    begin
      result := length( _counts );
    end
  ;

  function TPdfHistogram.getNBins(): integer;
    begin
      result := length( _counts );
    end
  ;


  function TPdfHistogram.probDensity( x: double ): double;
    begin
      if( nil <> @aphi_histogram_pdf ) then
        result := aphi_histogram_pdf( x, _cStructPtr )
      else
        result := NaN
      ;
    end
  ;


  function TPdfHistogram.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @aphi_histogram_cdf ) then
        result := aphi_histogram_cdf( x, _cStructPtr )
      else
        result := NaN
      ;
    end
  ;


  function TPdfHistogram.invCumulative( p: double ): double;
    begin
      if( nil <> @aphi_histogram_inverse_cdf ) then
        result := aphi_histogram_inverse_cdf( p, _cStructPtr )
      else
        result := NaN
      ;
    end
  ;


  function TPdfHistogram.rand(): double;
    begin
      if( nil <> @aphi_histogram_rand ) then
        result := aphi_histogram_rand( rngPtr(), _cStructPtr )
      else
        result := NaN
      ;
    end
  ;

  
  procedure TPdfHistogram.debug();
    var
      str: string;
      i: integer;
      errMsg: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Histogram PDF (id ' + intToStr( self.id ) + ')' + endl;

      if( not( validate( @errMsg ) ) ) then
        str := str + 'INVALID FUNCTION:' + endl + errMsg
      else
        begin
          str := str + 'Ranges, counts' + endl;

          for i := 0 to length( _ranges ) - 1 do
            begin
              if( i < length( _ranges ) - 1 ) then
                str := str +  '(' + usFloatToStr(_ranges[i]) + ', ' + usFloatToStr(_counts[i]) + ')' + endl
              else
                str := str +  '(' + usFloatToStr(_ranges[i]) + ')' + endl
              ;
            end
          ;
        end
      ;

      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfHistogram.createCopy(): TChartFunction;
    begin
      result := TPdfHistogram.create( self );
    end
  ;


  function TPdfHistogram.asCsv(): string;
    var
      i: integer;
    begin
      result := '';

      if( validate() ) then
        begin
          try
            result := 'x' + csvListSep() + ' y' + endl;

            for i := 0 to length( _counts ) - 1 do
              begin
                result :=
                  result
                  + csvFloatToStr( _ranges[i] )
                  + csvListSep() + ' '
                  + csvFloatToStr( _counts[i] )
                  + endl
                ;
              end
            ;

            result := result
              + csvFloatToStr( _ranges[high(_ranges)] )
              + csvListSep() + ' '
              + csvFloatToStr( 0.0 )
              + endl
            ;
          except
            raise exception.Create( 'An exception occurred in TPdfContinuous.asCsv()' );
            result := '';
          end;
        end
      ;
    end
  ;


  function TPdfHistogram.ssXml( const indent: integer ): string;
    var
      i: integer;
      str: string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <histogram>' + endl;

      for i := 0 to length( _counts ) - 1 do
        begin
          str := '    <value><x0>' + usFloatToStr( _ranges[i] ) +  '</x0>    <x1>' + usFloatToStr( _ranges[i+1] ) + '</x1>';
          str := str + '    <p>' + usFloatToStr( _counts[i] ) + '</p></value>';
          result := result + _indent + str + endl;
        end
      ;

      result := result + _indent + '  </histogram>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfHistogram.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
      newStyle: boolean;

      loRange, hiRange: TARDoubleArray;
      xmlCounts: TARDoubleArray;

      qRanges, qCounts: TQDoubleVector;
      ranges, counts: TARDoubleArray;

      success: boolean;
      i: integer;

      // Returns true on success
      function fillValuesFromOldXml(): boolean;
        var
          i: integer;
          count: integer;
          x0, x1, p: string;
          index: integer;
          triples: integer;
          ex0, ex1, ep: pointer;
        begin
          count := Sdew.GetElementCount( fnElement );

          result := false; // until shown otherwise

          if( 0 = ( count mod 3 ) ) then
            begin
              triples := count div 3;

              if( 0 < triples ) then
                begin
                  setLength( xmlCounts, triples );
                  setLength( loRange, triples );
                  setLength( hiRange, triples );

                  result := true;
                  index := 0;

                  for i := 1 to triples do
                    begin
                      ex0 := sdew.getElementByIndex( fnElement, index );
                      if( 'x0' = sdew.GetElementName( ex0 ) ) then
                        begin
                          x0 := sdew.getElementContents( ex0 );

                          ex1 := sdew.getElementByIndex( fnElement, index + 1 );
                          if( 'x1' = sdew.getElementName( ex1 ) ) then
                            begin
                              x1 := sdew.getElementContents( ex1 );

                              ep := sdew.GetElementByIndex( fnElement, index + 2 );
                              if( 'p' = sdew.GetElementName( ep ) ) then
                                begin
                                  p := sdew.getElementContents( ep );
                                  index := index + 3;
                                  loRange[i-1] := usStrToFloat( x0 );
                                  hiRange[i-1] := usStrToFloat( x1 );
                                  xmlCounts[i-1] := usStrToFloat( p );
                                end
                              else
                                begin
                                  dbcout( 'Bad name found were "p" was expected: ' + sdew.GetElementName( ep ), true );
                                  result := false;
                                  break;
                                  //Application.MessageBox('The histogram XML import method in this file has invalid parameters','ERROR');
                                end
                              ;
                            end
                          else
                            begin
                              dbcout( 'Bad name found were "x1" was expected: ' + sdew.GetElementName( ex1 ), true );
                              result := false;
                              break;
                              //Application.MessageBox('The histogram XML import method in this file has invalid parameters','ERROR');
                            end
                          ;
                        end
                      else
                        begin
                          dbcout( 'Bad name found were "x0" was expected: ' + sdew.GetElementName( ex0 ), true );
                          result := false;
                          break;
                          //Application.MessageBox('The histogram XML import method in this file has invalid parameters','ERROR');
                        end
                      ;
                    end
                  ; //  End For Loop
                end
              ;
            end
          ;
        end
      ;

      // Returns success on true
      function fillValuesFromNewXml(): boolean;
        var
          x0, x1, p: double;
          ev, ex0, ex1, ep: pointer;
          i, count: integer;
        begin
          result := true;

          count := Sdew.GetElementCount( fnElement );

          setLength( xmlCounts, count );
          setLength( loRange, count );
          setLength( hiRange, count );

          for i := 0 to count - 1 do
            begin
              ev := sdew.GetElementByIndex( fnElement, i );
              if( 'value' <> sdew.GetElementName( ev ) ) then
                begin
                  result := false;
                  break;
                end
              else
                begin
                  ex0 := sdew.GetElementByName( ev, 'x0' );
                  ex1 := sdew.GetElementByName( ev, 'x1' );
                  ep := sdew.GetElementByName( ev, 'p' );

                  if( nil <> ex0 ) then
                    x0 := usStrToFloat( sdew.GetElementContents( ex0 ), NaN )
                  else
                    x0 := NaN
                  ;

                  if( nil <> ex1 ) then
                    x1 := usStrToFloat( sdew.GetElementContents( ex1 ), NaN )
                  else
                    x1 := NaN
                  ;

                  if( nil <> ep ) then
                    p := usStrToFloat( sdew.GetElementContents( ep ), NaN )
                  else
                    p := NaN
                  ;

                  if( isNaN( x0 ) or isNaN( x1 ) or isNaN( p ) ) then
                    begin
                      result := false;
                      break;
                    end
                  else
                    begin
                      loRange[i] := x0;
                      hiRange[i] := x1;
                      xmlCounts[i] := p;
                    end
                  ;
                end
              ;
            end
          ;
        end
      ;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'histogram' );

      if( nil = fnElement ) then
        begin
          result := false;
          exit;
        end
      ;

      // Old-style XML may have the name attached to the type tag
      if( strIsEmpty( self.name ) ) then
        begin
          nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
          if( nil <> nameAttribute ) then
            name := sdew.GetElementAttribute( fnElement, 'name' )
          ;
        end
      ;

      newStyle := ( nil <> sdew.GetElementByName( fnElement, 'value' ) );

      if( newStyle ) then
        success := fillValuesFromNewXml()
      else
        success := fillValuesFromOldXml()
      ;

      if( success ) then
        begin
          // Once all of the XML has been parsed, the range borders need to be assembled to give
          // ranges that the rest of the class can use.  If ranges are completely consecutive,
          // this is pretty easy.  Unfortunately, there is no guarantee that this will be the case.
          qRanges := TQDoubleVector.create();
          qCounts := TQDoubleVector.create();

          qRanges.append( loRange[0] );
          qCounts.append( xmlCounts[0] );

          for i := 1 to length( xmlCounts ) - 1 do
            begin
              if( loRange[i] = hiRange[i-1] ) then
                begin
                  qRanges.append( loRange[i] );
                  qCounts.append( xmlCounts[i] );
                end
              else
                begin
                  qRanges.append( hiRange[i-1] );
                  qCounts.append( 0.0 );
                  qRanges.append( loRange[i] );
                  qCounts.append( xmlCounts[i] );
                end
              ;

            end
          ; // For loop

          // Don't forget the top range border!
          qRanges.append( hiRange[ length( hiRange ) - 1 ] );

          ranges := qRanges.createDoubleArray();
          counts := qCounts.createDoubleArray();

          setRangesAndCounts( ranges, counts );

          // Don't forget to clean up!
          qRanges.Free();
          qCounts.Free();

          setLength( ranges, 0 );
          setLength( counts, 0 );
        end
      ;

      setLength( xmlCounts, 0 );
      setLength( loRange, 0 );
      setLength( hiRange, 0 );

      if( success ) then
        result := true
      else
        begin
          initialize();
          result := false;
        end
      ;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfHistogram.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
    var
      q, q2: string;
      dict: TQueryDictionary;
      i: integer;
    begin
      dict := createDBFieldList( db );
      
      if( update ) then
        begin
          q := 'DELETE FROM `inChartDetail` WHERE `chartID` = ' + intToStr( id );
          db.execute( q );
        end
      ;

      dict['chartType']   := '"Histogram"';

      if( not update ) then
        q := writeQuery( 'inChart', QInsert, dict )
      else
        q := writeQuery( 'inChart', QUpdate, dict, 'WHERE [chartID] = ' + intToStr( id ) );
      ;

      dict.Clear();
      freeAndNil( dict );

      db.execute( q );

      if( not update ) then id := db.lastInsertID();

      // Populate the inChartDetail with the individual points
      q := 'INSERT INTO [inChartDetail] ([chartID], [pointOrder], [x], [y])';

      for i := 0 to length( _ranges ) - 2 do
        begin
          q2 := q + ' VALUES (' + usFloatToStr( id ) + ', '+ intToStr( i ) + ', ' + usFloatToStr( _ranges[i] ) + ', ' + usFloatToStr( _counts[i] ) + ')';
          db.execute( q2 );
        end
      ;

      // There is one more X value than there are Y values.
      // The last Y value in the database will be null.
      q2 := q + ' VALUES (' + usFloatToStr( id ) + ', '+ intToStr( length( _ranges ) - 1 ) + ', ' + usFloatToStr( _ranges[length( _ranges ) - 1] ) + ', NULL )';
      db.execute( q2 );

      result := id;
    end
  ;
  {$ENDIF}
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Inverse Gaussian PDF
//-----------------------------------------------------------------------------
  constructor TPdfInverseGaussian.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfInverseGaussian.create( const src: TPdfInverseGaussian );
    begin
      dbcout( '*** Calling TPdfInverseGaussian copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setMean( src.mean );
      setShape( src.shape );
      xUnits := src._xUnits;
    end
  ;


  constructor TPdfInverseGaussian.create( mean, shape: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setMean( mean );
      setShape( shape );
      xUnits := u;
    end
  ;


  procedure TPdfInverseGaussian.initialize();
    begin
      setPdfType( PdfInverseGaussian );

      setMean( NaN );
      setShape( NaN );

      _infLeftTail := false;
      _infRightTail := true;
    end
  ;


  destructor TPdfInverseGaussian.destroy();
    begin
      //dbcout( 'Destroying TPdfInverseGaussian ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;


  function TPdfInverseGaussian.compare( Chart:TChartFunction ): boolean;
    begin
      result := false;

      if ( inherited compare( Chart ) ) then
        begin
          if
            ( self._mean = (Chart as TPdfInverseGaussian)._mean )
          and
             ( self._shape = (Chart as TPdfInverseGaussian)._shape )
          then
            result := true
          ;
        end
      ;
    end
  ;


  function TPdfInverseGaussian.getMean(): double; begin Result := _mean; end;
  function TPdfInverseGaussian.getShape(): double; begin Result := _shape; end;
  procedure TPdfInverseGaussian.setMean( val: double ); begin _mean := val; freeDiscreteValues(); end;
  procedure TPdfInverseGaussian.setShape( val: double ); begin _shape := val; freeDiscreteValues(); end;


  function TPdfInverseGaussian.getMin(): double;
    begin
      result := 0;
    end
  ;


  function TPdfInverseGaussian.getMax(): double;
    begin
      result := infinity;
    end
  ;


  function TPdfInverseGaussian.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfInverseGaussian.validate( err: pstring = nil ): boolean;
    begin
      if( isNaN( mean ) or isNan( shape ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( mean <= 0 ) or ( shape <= 0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Values must be greater than 0.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfInverseGaussian.getDescr(): string;
    begin
      result := tr( 'Inverse Gaussian' ) + format( ' ( %f, %f )', [_mean, _shape] );
    end
  ;


  function TPdfInverseGaussian.probDensity( x: double ): double;
    (*
    var
      t1, t2: double;
    *)
    begin
      if( nil <> @aphi_inverse_gaussian_pdf ) then
        result := aphi_inverse_gaussian_pdf( x, self.mean, self.shape )
      else
        result := NaN
      ;
      (*
      t1 := shape / ( 2 * pi * power( x, 3 ) );

      t2 := x - mean;
      t2 := t2 * t2;
      t2 := shape * t2;
      t2 := exp( t2 / ( 2 * mean * mean * x ) );

      result := sqrt( t1 * t2 );
      *)
    end
  ;


  function TPdfInverseGaussian.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @aphi_inverse_gaussian_cdf ) then
        result := aphi_inverse_gaussian_cdf( x, self.mean, self.shape )
      else
        result := NaN
      ;
    end
  ;


  function TPdfInverseGaussian.invCumulative( p: double ): double;
    begin
      if( nil <> @aphi_inverse_gaussian_inverse_cdf ) then
        result := aphi_inverse_gaussian_inverse_cdf( p, self.mean, self.shape )
      else
        result := NaN
      ;
    end
  ;


  function TPdfInverseGaussian.rand(): double;
    var
      x, y, n, ymean, u: double;
    begin
      // Algorithm is from:
      // Devroye, L.  1986.  Non-Uniform Random Variate Generation. New York: Springer-Verlag.
      // See http://cg.scs.carleton.ca/~luc/rnbookindex.html
      n := TPdfGaussian.rand( 0.0, 1.0 );
      y := n * n;

      ymean := y * mean;

      x := (
        mean
        + ((mean * ymean) / (2.0 * shape))
        - ((mean / (2.0 * shape)) * sqrt (4.0 * ymean * shape + ymean * ymean ))
      );

      u := rngRand();

      if( ( mean/(mean + x) ) >= u ) then
        result := x
      else
        result := (mean*mean)/x
      ;
    end
  ;


  procedure TPdfInverseGaussian.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Inverse Gaussian PDF (id ' + intToStr( self.id ) + ') : mean=' + usFloatToStr(self.mean) + ', ' + endl + 'shape=' + usFloatToStr(shape);
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfInverseGaussian.createCopy(): TChartFunction;
    begin
      result := TPdfInverseGaussian.create( self );
    end
  ;


  function TPdfInverseGaussian.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <inverse-gaussian>' + endl;
      result := result + _indent + '    <mu>' + usFloatToStr( mean ) + '</mu>' + endl;
      result := result + _indent + '    <lambda>' + usFloatToStr( shape ) + '</lambda>' + endl;
      result := result + _indent + '  </inverse-gaussian>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfInverseGaussian.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'inverse-gaussian' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setMean( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'mu' )), NaN ) );
          setShape( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'lambda' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfInverseGaussian.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"InverseGaussian"';
      dict['mean']        := usFloatToStr( mean );
      dict['shape']      := usFloatToStr( shape );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Logistic PDF
//-----------------------------------------------------------------------------
  constructor TPdfLogistic.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfLogistic.create( location, scale: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setLocation( location );
      setScale( scale );
      xUnits := u;
    end
  ;


  constructor TPdfLogistic.create( const src: TPdfLogistic );
    begin
      dbcout( '*** Calling TPdfLogistic copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setLocation( src.location );
      setScale( src.scale );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfLogistic.initialize();
    begin
      setPdfType( PdfLogistic );

      setLocation( NaN );
      setScale( NaN );

      _infLeftTail := true;
      _infRightTail := true;
    end
  ;


  destructor TPdfLogistic.destroy();
    begin
      //dbcout( 'Destroying TPdfLogistic ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfLogistic.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._location = (Chart as TPdfLogistic)._location ) and
             ( self._scale = (Chart as TPdfLogistic)._scale ) then
            result := true;
        end;
    end;

  procedure TPdfLogistic.setLocation( val: double ); begin _location := val; freeDiscreteValues(); end;
  procedure TPdfLogistic.setScale( val: double ); begin _scale := val; freeDiscreteValues(); end;
  function TPdfLogistic.getLocation(): double; begin result := _location; end;
  function TPdfLogistic.getScale(): double; begin result := _scale; end;


  function TPdfLogistic.validate( err: pstring = nil ): boolean;
    begin
      if( isNaN( location ) or isNaN( scale ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( (* ( location <= 0.0 ) or *) ( scale <= 0.0 ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'Scale must be greater than 0.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfLogistic.probDensity( x: double ): double;
    begin
      if( nil <> @gslRanLogisticPdf ) then
        result := gslRanLogisticPdf( x - self.location, self.scale )
      else
        result := NaN
      ;
    end
  ;


  function TPdfLogistic.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @gslCdfLogisticP ) then
        result := gslCdfLogisticP( x - self.location, self.scale )
      else
        result := NaN
      ;
    end
  ;


  function TPdfLogistic.invCumulative( p: double ): double;
    begin
      if( nil <> @gslCdfLogisticPinv ) then
        result := gslCdfLogisticPInv( p, self.scale ) + self.location
      else
        result := NaN
      ;
    end
  ;


  function TPdfLogistic.rand(): double;
    begin
      if( nil <> @gslRanLogistic ) then
        result := gslRanLogistic( gslRngPtr(), self.scale ) + self.location
      else
        result := NaN
      ;
    end
  ;


  function TPdfLogistic.getDescr(): string;
    begin
      result := tr( 'Logistic' ) + format( ' ( %f, %f )', [location, scale] );
    end
  ;


  function TPdfLogistic.getMean(): double;
    begin
      result := location;
    end
  ;


  function TPdfLogistic.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfLogistic.getMin(): double;
    begin
      result := negInfinity;
    end
  ;

  function TPdfLogistic.getMax(): double;
    begin
      result := infinity;
    end
  ;


  procedure TPdfLogistic.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str
        + 'Logistic PDF (id ' + intToStr( self.id ) + ') :' + endl
        + ', location=' + usFloatToStr(location)
        + ', scale=' + usFloatToStr(scale)
      ;

      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfLogistic.createCopy(): TChartFunction;
    begin
      result := TPdfLogistic.create( self );
    end
  ;


 function TPdfLogistic.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <logistic>' + endl;
      result := result + _indent + '    <location>' + usFloatToStr( location ) + '</location>' + endl;
      result := result + _indent + '    <scale>' + usFloatToStr( scale ) + '</scale>' + endl;
      result := result + _indent + '  </logistic>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfLogistic.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'logistic' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ; 

          setLocation( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'location' )), NaN ) );
          setScale( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'scale' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfLogistic.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Logistic"';
      dict['location']    := usFloatToStr( location );
      dict['scale']       := usFloatToStr( scale );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Loglogistic PDF
//-----------------------------------------------------------------------------
  constructor TPdfLogLogistic.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;

  constructor TPdfLogLogistic.create( location, scale, shape: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setLocation( location );
      setScale( scale );
      setShape( shape );
      xUnits := u;
    end
  ;


  constructor TPdfLogLogistic.create( const src: TPdfLogLogistic );
    begin
      dbcout( '*** Calling TPdfLogLogistic copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setLocation( src.location );
      setScale( src.scale );
      setShape( src.shape );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfLogLogistic.initialize();
    begin
      setPdfType( PdfLogLogistic );

      setLocation( NaN );
      setScale( NaN );
      setShape( NaN );

      _infLeftTail := false;
      _infRightTail := true;
    end
  ;


  destructor TPdfLogLogistic.destroy();
    begin
      //dbcout( 'Destroying TPdfLogLogistic ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfLogLogistic.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._location = (Chart as TPdfLogLogistic)._location ) and
             ( self._scale = (Chart as TPdfLogLogistic)._scale ) and
             ( self._shape = (Chart as TPdfLogLogistic)._shape ) then
            result := true;
        end;
    end;

  procedure TPdfLogLogistic.setLocation( val: double ); begin _location := val; freeDiscreteValues(); end;
  procedure TPdfLogLogistic.setScale( val: double ); begin _scale := val; freeDiscreteValues(); end;
  procedure TPdfLogLogistic.setShape( val: double ); begin _shape := val; freeDiscreteValues(); end;
  function TPdfLogLogistic.getLocation(): double; begin result := _location; end;
  function TPdfLogLogistic.getScale(): double; begin result := _scale; end;
  function TPdfLogLogistic.getShape(): double; begin result := _shape; end;


  function TPdfLogLogistic.validate( err: pstring = nil ): boolean;
    begin
      if( isNaN( location ) or isNaN( scale ) or isNaN( shape ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( ( scale <= 0.0 ) or ( shape <= 0.0 ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'Scale and shape must be greater than 0.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfLogLogistic.getDescr(): string;
    begin
      result := tr( 'LogLogistic' ) + format( ' ( %f, %f, %f )', [location, scale, shape] );
    end
  ;


  function TPdfLogLogistic.getMean(): double;
    var
      theta: double;
    begin
      if( shape <= 1 ) then
        begin
          raise exception.create( 'Mean cannot be calculated for shape <= 1 in TPdfLogLogistic.getMean()' );
          result := NaN;
        end
      else
        begin
          theta := pi / shape;
          result := location + ( scale * theta * cosecant( theta ) );
        end
      ;
    end
  ;


  function TPdfLogLogistic.getHasMean(): boolean;
    begin
      result := ( shape > 1 );
    end
  ;

  function TPdfLogLogistic.getMin(): double;
    begin
      result := location;
    end
  ;

  function TPdfLogLogistic.getMax(): double;
    begin
      result := infinity;
    end
  ;


  function TPdfLogLogistic.probDensity( x: double ): double;
    begin
        if( nil <> @aphi_loglogistic_pdf ) then
          result := aphi_loglogistic_pdf( x, self.location, self.scale, self.shape )
        else
          result := NaN
        ;
    end
  ;

  function TPdfLogLogistic.cumulativeDistr( x: double ): double;
    begin
        if( nil <> @aphi_loglogistic_cdf ) then
          result := aphi_loglogistic_cdf( x, self.location, self.scale, self.shape )
        else
          result := NaN
        ;
    end
  ;


  function TPdfLogLogistic.invCumulative( p: double ): double;
    begin
        if( nil <> @aphi_loglogistic_inverse_cdf ) then
          result := aphi_loglogistic_inverse_cdf( p, self.location, self.scale, self.shape )
        else
          result := NaN
        ;
    end
  ;


  function TPdfLogLogistic.rand(): double;
    begin
      if( nil <> @aphi_loglogistic_rand ) then
        result := aphi_loglogistic_rand( rngPtr(), self.location, self.scale, self.shape )
      else
        result := NaN
      ;
    end
  ;


  procedure TPdfLogLogistic.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str
        + 'LogLogistic PDF (id ' + intToStr( self.id ) + ') : '
        + descr
      ;

      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfLogLogistic.createCopy(): TChartFunction;
    begin
      result := TPdfLogLogistic.create( self );
    end
  ;


  function TPdfLogLogistic.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <loglogistic>' + endl;
      result := result + _indent + '    <location>' + usFloatToStr( location ) + '</location>' + endl;
      result := result + _indent + '    <scale>' + usFloatToStr( scale ) + '</scale>' + endl;
      result := result + _indent + '    <shape>' + usFloatToStr( shape ) + '</shape>' + endl;
      result := result + _indent + '  </loglogistic>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfLogLogistic.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'loglogistic' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setLocation( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'location' )), NaN ) );
          setScale( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'scale' )), NaN ) );
          setShape( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'shape' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfLogLogistic.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Loglogistic"';
      dict['location']    := usFloatToStr( location );
      dict['scale']       := usFloatToStr( scale );
      dict['shape']       := usFloatToStr( shape );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Lognormal PDF
//-----------------------------------------------------------------------------
  constructor TPdfLognormal.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfLognormal.create( mean, stddev: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setMean( mean );
      setStddev( stddev );
      xUnits := u;
    end
  ;


  constructor TPdfLognormal.create( const src: TPdfLognormal );
    begin
      dbcout( '*** Calling TPdfLognormal copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setMean( src.mean );
      setStddev( src.stddev );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfLognormal.initialize();
    begin
      setPdfType( PdfLognormal );

      setMean( NaN );
      setStddev( NaN );

      _infLeftTail := false;
      _infRightTail := true;
    end
  ;


  destructor TPdfLognormal.destroy();
    begin
      //dbcout( 'Destroying TPdfLognormal ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfLognormal.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._zeta = (Chart as TPdfLognormal)._zeta ) and
             ( self._mean = (Chart as TPdfLognormal)._mean ) and
             ( self._stddev = (Chart as TPdfLognormal)._stddev ) and
             ( self._sigma = (Chart as TPdfLognormal)._sigma ) then
            result := true;
        end;
    end;

  procedure TPdfLognormal.setZeta( val: double );
    begin
      _zeta := val;
      _mean := calculateMean( _zeta, _sigma );
      _stddev := calculateStddev( _zeta, _sigma );
      freeDiscreteValues(); 
    end
  ;


  procedure TPdfLognormal.setSigma( val: double );
    begin
      _sigma := val;
      _mean := calculateMean( _zeta, _sigma );
      _stddev := calculateStddev( _zeta, _sigma );
      freeDiscreteValues();
    end
  ;


  procedure TPdfLognormal.setMean( val: double );
    begin
      _mean := val;
      _zeta := calculateZeta( _mean, _stddev );
      _sigma := calculateSigmaPrime( _mean, _stddev );
      freeDiscreteValues();
    end
  ;


  procedure TPdfLognormal.setStddev( val: double );
    begin
      _stddev := val;
      _zeta := calculateZeta( _mean, _stddev );
      _sigma := calculateSigmaPrime( _mean, _stddev );
      freeDiscreteValues(); 
    end
  ;


  class function TPdfLognormal.calculateMean( zeta, sigma: double ): double;
    begin
      if( isNaN( zeta ) or isNaN( sigma ) ) then
        result := NaN
      else if( 0 >= sigma ) then
        result := NaN
      else
        result := exp( zeta + (sqr(sigma)/2) )
      ;
    end
  ;


  class function TPdfLognormal.calculateStddev( zeta, sigma: double ): double;
    var
      omega: double;
    begin
      if( isNaN( zeta ) or isNaN( sigma ) ) then
        result := NaN
      else if( 0 >= sigma ) then
        result := NaN
      else
        begin
          omega := exp( sqr(sigma) );
          result := sqrt( exp(2*zeta) * omega * (omega - 1) );
        end
      ;
    end
  ;


  class function TPdfLognormal.calculateZeta( mean, stddev: double ): double;
    begin
      if( isNaN( mean ) or isNan( stddev ) ) then
        result := NaN
      else if( ( 0 >= mean ) or ( 0 >= stddev ) ) then
        result := NaN
      else
        result := ln( sqr(mean) / sqrt( sqr(stddev) + sqr(mean) ) )
      ;
    end
  ;

  class function TPdfLognormal.calculateSigmaPrime( mean, stddev: double ): double;
    begin
      if( isNaN( mean ) or isNan( stddev ) ) then
        result := NaN
      else if( ( 0 >= mean ) or ( 0 >= stddev ) ) then
        result := NaN
      else
        result := sqrt( ln( 1 + sqr(stddev/mean) ) )
      ;
    end
  ;



  function TPdfLognormal.getZeta(): double; begin Result := _zeta; end;
  function TPdfLognormal.getSigma(): double; begin Result := _sigma; end;
  function TPdfLognormal.getMean(): double; begin Result := _mean; end;
  function TPdfLognormal.getStddev(): double; begin Result := _stddev; end;

  function TPdfLognormal.getHasMean(): boolean;
    begin
      result := true;
    end
  ;

  function TPdfLognormal.getMin(): double;
    begin
      result := 0.0;
    end
  ;

  function TPdfLognormal.getMax(): double;
    begin
      result := infinity;
    end
  ;

  function TPdfLognormal.validate( err: pstring = nil ): boolean;
    begin
      if
        ( isNaN( getMean() ) )
      or
        ( isNaN( getStddev() ) )
      or
        ( isNaN( getSigma() ) )
      or
        ( isNaN( getZeta() ) )
      then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
          exit;
        end
      ;


      if
        ( 0.0 >= getSigma() )
      or
        ( 0.0 >= getMean() )
      or
        ( 0.0 >= getStddev() )
      then
        begin
          if( err <> nil ) then err^ := tr( 'Mean, standard deviation, and sigma'' must be greater than 0.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfLognormal.probDensity( x: double ): double;
    begin
      if( nil <> @gslRanLognormalPdf ) then
        result := gslRanLognormalPdf( x, self.zeta, self.sigma )
      else
        result := NaN
      ;
    end
  ;


  function TPdfLognormal.cumulativeDistr( x: double ): double;
    begin
      //dbcout2( usFloatToStr( x ) + ', ' + usFloatToStr( zeta ) + ', ' + usFloatToStr( sigma ) );
      if( self.min = x ) then
        result := 0.0
      else if( nil <> @gslCdfLognormalP ) then
        result := gslCdfLognormalP( x, self.zeta, self.sigma )
      else
        result := NaN
      ;
    end
  ;


  function TPdfLognormal.invCumulative( p: double ): double;
    begin
      if( nil <> @gslCdfLognormalPinv ) then
        result := gslCdfLognormalPInv( p, self.zeta, self.sigma )
      else
        result := NaN
      ;
    end
  ;


  function TPdfLognormal.rand(): double;
    begin
      if( nil <> @gslRanLognormal ) then
        result := gslRanLognormal( gslRngPtr(), self.zeta, self.sigma )
      else
        result := NaN
      ;
    end
  ;


  function TPdfLognormal.getDescr(): string;
    begin
      result := tr( 'Lognormal' ) + format( ' ( %f, %f )', [mean, stddev] );
    end
  ;


  procedure TPdfLognormal.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Lognormal PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 'mean = ' + usFloatToStr( self.mean )
        + ', stddev = ' + usFloatToStr( self.stddev )
        + ', zeta = ' + usFloatToStr( self.zeta )
        + ', sigma = ' + usFloatToStr( self.sigma )
      ;

      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfLognormal.createCopy(): TChartFunction;
    begin
      result := TPdfLognormal.create( self );
    end
  ;


  function TPdfLognormal.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;

      result := result + _indent + '  <lognormal>' + endl;
      result := result + _indent + '    <!-- The NAADSM interface calculates this zeta and sigma from the -->' + endl;
      result := result + _indent + '    <!-- stored mean and standard deviation, as described in these references: -->' + endl;
      result := result + _indent + '    <!-- 1) Vose, D. 1996.  Quantitative risk analysis: a guide to Monte Carlo Simulation Modelling. -->' + endl;
      result := result + _indent + '    <!--    John Wiley and Sons, New York. -->' + endl;
      result := result + _indent + '    <!-- 2) Palisade Corporation. 2002.  A concise summary of @RISK probability distribution functions. -->' + endl;
      result := result + _indent + '    <!--    Available (probably illegally) from http://project.zf.jcu.cz/risk/data/distfunc.pdf -->' + endl;
      result := result + _indent + '    <zeta>' + usFloatToStr( self.zeta ) + '</zeta>' + endl;
      result := result + _indent + '    <sigma>' + usFloatToStr( self.sigma ) + '</sigma>' + endl;
      result := result + _indent + '  </lognormal>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfLognormal.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'lognormal' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setZeta( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'zeta' )), NaN ) );
          setSigma( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'sigma' )), NaN ) );

          if ( (  not IsNaN( _zeta )  ) and ( not IsNaN( _sigma ) ) ) then
            begin
              setMean( calculateMean( _zeta, _sigma )  );
              setStddev( calculateStddev( _zeta, _sigma ) );
            end
          else
            initialize();
          ;
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfLognormal.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Lognormal"';
      dict['mean']        := usFloatToStr( self.mean );
      dict['stddev']      := usFloatToStr( self.stddev );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Pareto PDF
//-----------------------------------------------------------------------------
  constructor TPdfPareto.create();
    begin
      inherited create();
      // initialize() is called by the inherited contructor
    end
  ;


  constructor TPdfPareto.create( theta, a: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setTheta( theta );
      setA( a );
      xUnits := u;
    end
  ;


  constructor TPdfPareto.create( const src: TPdfPareto );
    begin
      dbcout( '*** Calling TPdfPareto copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setTheta( src.theta );
      setA( src.a );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfPareto.initialize();
    begin
      setPdfType( PdfPareto );

      setTheta( NaN );
      setA( NaN );

      _infLeftTail := true;
      _infRightTail := true;
    end
  ;


  destructor TPdfPareto.destroy();
    begin
      //dbcout( 'Destroying TPdfPareto ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfPareto.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._theta = (Chart as TPdfPareto)._theta ) and
             ( self._a = (Chart as TPdfPareto)._a ) then
            result := true;
        end;
    end;

  procedure TPdfPareto.setTheta( val: double ); begin _theta := val; freeDiscreteValues(); end;
  procedure TPdfPareto.setA( val: double ); begin _a := val; freeDiscreteValues(); end;
  function TPdfPareto.getTheta(): double; begin result := _theta; end;
  function TPdfPareto.getA(): double; begin result := _a; end;


  function TPdfPareto.validate( err: pstring = nil ): boolean;
    begin
      if( isNaN( theta ) or isNaN( a ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( theta <= 0.0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Theta must be greater than 0.' );
          result := false;
        end
      else if( a <= 0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'a must be greater than 0.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfPareto.getDescr(): string;
    begin
      result := tr( 'Pareto' ) + format( ' ( %f, %f )', [theta, a] );
    end
  ;


  function TPdfPareto.getMean(): double;
    begin
      if( 1 >= theta ) then
        begin
          raise exception.create( 'Mean cannot be calculated with theta <= 1 in TPdfPareto.getMean()' );
          result := NaN;
        end
      else
        result := theta * a / (theta - 1.0)
      ;
    end
  ;

  function TPdfPareto.getHasMean(): boolean;
    begin
      result := ( 1 < theta );
    end
  ;


  function TPdfPareto.getMin(): double;
    begin
      result := a;
    end
  ;

  function TPdfPareto.getMax(): double;
    begin
      result := infinity;
    end
  ;


  function TPdfPareto.probDensity( x: double ): double;
    begin
      if( nil <> @gslRanParetoPdf ) then
        result := gslRanParetoPdf( x, self.theta, self.a )
      else
        result := NaN
      ;
    end
  ;


  function TPdfPareto.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @gslCdfParetoP ) then
        result := gslCdfParetoP( x, self.theta, self.a )
      else
        result := NaN
      ;
    end
  ;


  function TPdfPareto.invCumulative( p: double ): double;
    begin
      if( nil <> @gslCdfParetoPinv ) then
        result := gslCdfParetoPinv( p, self.theta, self.a )
      else
        result := NaN
      ;
    end
  ;


  function TPdfPareto.rand(): double;
    begin
      if( nil <> @gslRanPareto ) then
        result := gslRanPareto( gslRngPtr(), self.theta, self.a )
      else
        result := NaN
      ;
    end
  ;


  procedure TPdfPareto.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Pareto PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 'shape=' + usFloatToStr( theta )
        + ', location=' + usFloatToStr( a )
      ;
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfPareto.createCopy(): TChartFunction;
    begin
      result := TPdfPareto.create( self );
    end
  ;


  function TPdfPareto.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <pareto>' + endl;
      result := result + _indent + '    <theta>' + usFloatToStr( theta ) + '</theta>' + endl;
      result := result + _indent + '    <a>' + usFloatToStr( a ) + '</a>' + endl;
      result := result + _indent + '  </pareto>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfPareto.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'pareto' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setTheta( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'theta' )), NaN ) );
          setA( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'a' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfPareto.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Pareto"';
      dict['a']       := usFloatToStr( a );
      dict['theta']        := usFloatToStr( theta );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Pearson5 PDF
//-----------------------------------------------------------------------------
  constructor TPdfPearson5.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfPearson5.create( alpha, beta: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setAlpha( alpha );
      setBeta( beta );
      xUnits := u;
    end
  ;


  constructor TPdfPearson5.create( const src: TPdfPearson5 );
    begin
      dbcout( '*** Calling TPdfPearson5 copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setAlpha( src.alpha );
      setBeta( src.beta );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfPearson5.initialize();
    begin
      setPdfType( PdfPearson5 );

      setAlpha( NaN );
      setBeta( NaN );

      _infLeftTail := true;
      _infRightTail := true;
    end
  ;


  destructor TPdfPearson5.destroy();
    begin
      //dbcout( 'Destroying TPdfPearson5 ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfPearson5.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._alpha = (Chart as TPdfPearson5)._alpha ) and
             ( self._beta = (Chart as TPdfPearson5)._beta ) then
            result := true;
        end;
    end;

  procedure TPdfPearson5.setAlpha( val: double ); begin _alpha := val; freeDiscreteValues(); end;
  procedure TPdfPearson5.setBeta( val: double ); begin _beta := val; freeDiscreteValues(); end;
  function TPdfPearson5.getAlpha(): double; begin Result := _alpha; end;
  function TPdfPearson5.getBeta(): double; begin Result := _beta; end;


  function TPdfPearson5.validate( err: pstring = nil ): boolean;
    begin
      if( isNaN( alpha ) or isNaN( beta ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( beta <= 0.0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Beta must be greater than 0.' );
          result := false;
        end
      else if( alpha <= 1 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Alpha must be greater than 1.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfPearson5.getDescr(): string;
    begin
      result := tr( 'Pearson 5' ) + format( ' ( %f, %f )', [alpha, beta] );
    end
  ;


  function TPdfPearson5.getMean(): double;
    begin
      if( 1 >= alpha ) then
        begin
          raise exception.create( 'Mean cannot be calculated with alpha <= 1 in TPdfPearson5.getMean()' );
          result := NaN;
        end
      else
        result := beta / ( alpha - 1 )
      ;
    end
  ;

  function TPdfPearson5.getHasMean(): boolean;
    begin
      result := ( 1 < alpha );
    end
  ;


  function TPdfPearson5.getMin(): double;
    begin
      result := negInfinity;
    end
  ;

  function TPdfPearson5.getMax(): double;
    begin
      result := infinity;
    end
  ;


  function TPdfPearson5.probDensity( x: double ): double;
    begin
      if( nil <> @aphi_pearson5_pdf ) then
        result := aphi_pearson5_pdf( x, self.alpha, self.beta )
      else
        result := NaN
      ;
    end
  ;


  function TPdfPearson5.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @aphi_pearson5_cdf ) then
        result := aphi_pearson5_cdf( x, self.alpha, self.beta )
      else
        result := NaN
      ;
    end
  ;


  function TPdfPearson5.invCumulative( p: double ): double;
    begin
      if( nil <> @aphi_pearson5_inverse_cdf  ) then
        result :=aphi_pearson5_inverse_cdf( p, self.alpha, self.beta )
      else
        result := NaN
      ;
    end
  ;


  function TPdfPearson5.rand(): double;
    begin
      if( nil <> @aphi_pearson5_rand ) then
        result := aphi_pearson5_rand( rngPtr(), self.alpha, self.beta )
      else
        result := NaN
      ;
    end
  ;


  procedure TPdfPearson5.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Pearson5 PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 'alpha=' + usFloatToStr(alpha)
        + ', beta=' + usFloatToStr(beta)
      ;
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfPearson5.createCopy(): TChartFunction;
    begin
      result := TPdfPearson5.create( self );
    end
  ;


  function TPdfPearson5.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <pearson5>' + endl;
      result := result + _indent + '    <alpha>' + usFloatToStr( alpha ) + '</alpha>' + endl;
      result := result + _indent + '    <beta>' + usFloatToStr( beta ) + '</beta>' + endl;
      result := result + _indent + '  </pearson5>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfPearson5.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'pearson5' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setAlpha( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'alpha' )), NaN ) );
          setBeta( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'beta' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfPearson5.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Pearson5"';
      dict['alpha']       := usFloatToStr( alpha );
      dict['beta']        := usFloatToStr( beta );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Piecewise PDF
//-----------------------------------------------------------------------------
  constructor TPdfPiecewise.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfPiecewise.create( pnts: RPointArray; u: TChartUnitType = UUnknown );
    var
      i: integer;
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setLength( _points, length( pnts ) );

      for i := 0 to length( pnts ) - 1 do
        begin
          _points[i].x := pnts[i].x;
          _points[i].y := pnts[i].y;
        end
      ;

      xUnits := u;
      
      // Standardize the points to give an area of 1
      standardize( _points );
      calculateSlopeAndCumul();
    end
  ;

  
  constructor TPdfPiecewise.create( pnts: T2DPointList; u: TChartUnitType = UUnknown );
    var
      i: integer;
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      xUnits := u;

      setLength( _points, pnts.Count );

      for i := 0 to pnts.Count - 1 do
        begin
          _points[i].x := pnts[i].x;
          _points[i].y := pnts[i].y;
        end
      ;

      standardize( _points );
      calculateSlopeAndCumul();
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  constructor TPdfPiecewise.create( db: TSqlDatabase; chartID: integer; u: TChartUnitType = UUnknown );
    var
      q: string;
      res: TSqlResult;
      row: TSqlRow;
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      q := 'SELECT [x], [y] FROM [inChartDetail] WHERE [chartID] = ' + intToStr( chartID ) + ' ORDER BY [pointOrder]';

      res := TSqlResult.create( q, db );

      row := res.fetchArrayFirst();

      while( row <> nil ) do
        begin
          setLength( _points, length(_points)+1 );
          _points[high(_points)].x := double( row.field(0) );
          _points[high(_points)].y := double( row.field(1) );
          row := res.fetchArrayNext();
        end
      ;

      setXUnits( u );


      freeAndNil( res );

      // Standardize the points to give an area of 1
      standardize( _points );
      calculateSlopeAndCumul();
    end
  ;
  {$ENDIF}
  

  constructor TPdfPiecewise.create( const src: TPdfPiecewise );
    begin
      dbcout( '*** Calling TPdfPiecewise copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setLength( _points, length( src._points ) );
      _points := copy( src._points, 0, length( src._points ) );

      xUnits := src.xUnits;
    end
  ;


  procedure TPdfPiecewise.initialize();
    begin
      setPdfType( PdfPiecewise );

      _infLeftTail := false;
      _infRightTail := false;

      clearPoints();
    end
  ;


  destructor TPdfPiecewise.destroy();
    begin
      //dbcout( 'Destroying TPdfPiecewise ' + self.name, DBSHOWMSG );
      clearPoints();
      inherited destroy();
    end
  ;


  function TPdfPiecewise.compare( Chart: TChartFunction ):boolean;
    var
      I, count: integer;
    begin
      result := false;

      if ( inherited compare( Chart ) ) then
        begin
          if ( length( self._points ) = (Chart as TPdfPiecewise ).getPointCount() )  then
            begin
               Count := length( self._points );

               result := true;
               for I := 0 to Count - 1 do
                 begin
                   if ( ( self._points[I].x <> ( Chart as TPdfPiecewise ).getPoint( I ).x ) or ( self._points[I].y <> ( Chart as TPdfPiecewise ).getPoint( I ).y ) ) then
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


  //-----------------------------------------------------------------------------
  //  Function for standardizing the area under the curve
  //-----------------------------------------------------------------------------
  class function TPdfPiecewise.standardize( var points: RPointArray; msg: PString = nil ): boolean;
    var
      i: integer;
      cumul: TARDoubleArray;
      width: double;
      total: double;

      // Checks that each successive X value is greater than the preceding one
      //----------------------------------------------------------------------
      function successiveXValuesOK( var points: RPointArray ): boolean;
        var
          i: integer;
          prevX: double;
        begin
          result := true;
          prevX := points[0].x;

          for i := 1 to high(points) do
            begin
              if( points[i].x <= prevX ) then
                begin
                  result := false;
                  break;
                end
              else
                prevX := points[i].x
              ;
            end
          ;
        end
      ;
      
      // Called recursively to remove multiple zeroes from the end of the array
      //-----------------------------------------------------------------------
      procedure trimEndPoints( var p: RPointArray );
        var
          pos: integer;
        begin
          pos := high( p ) - 1;

          if( sameValue( p[pos].y, 0.0, 0.00001 ) ) then
            begin
              setLength( p, length( p ) - 1 );
              p[high(p)].y := 0.0;
              
              trimEndPoints( p );
            end
          ;
        end
      ;


      // Called recursively to remove multiple zeroes from the start of the array
      //-------------------------------------------------------------------------
      procedure trimStartPoints( var p: RPointArray );
        var
          arr: RPointArray;
        begin
          if( sameValue( p[1].y, 0.0, 0.00001 ) ) then
            begin
              arr := copy( p, 1, length(p) - 1 );
              setLength( p, 0 );
              p := copy( arr );
              p[0].y := 0.0;

              trimStartPoints( p );
            end
          ;
        end
      ;

    begin // standardize
      result := true; // until proven otherwise
      
      // If there are fewer than 3 points, it's hopeless.  Exit now.
      //------------------------------------------------------------
      if( 3 > length( points ) ) then
        begin
          result := false;

          if ( nil <> msg ) then
            msg^ := tr( 'There are too few points to specify this function.' )
          ;
          exit;
        end
      ;

      // Make sure that 1st and last points are 0.
      //------------------------------------------
      if( 0.0 <> points[Low(points)].Y ) then
        begin
          result := false;

          if ( nil <> msg ) then
            msg^ :=tr( 'The first Y value of a piecewise PDF must be 0.' )
          ;
          exit;
        end
      ;
      if( 0.0 <> points[High(points)].Y ) then
        begin
          result := false;

          if ( nil <> msg ) then
            msg^ := tr( 'The last Y value of a piecewise PDF must be 0.' )
          ;
          exit;
        end
      ;

      // Trim excess 0's from the head and tail of the points list
      //----------------------------------------------------------
      trimStartPoints( points );
      trimEndPoints( points );

      // If there are fewer than 3 points after trimming, there is a problem.
      if( 3 > length( points ) ) then
        begin
          result := false;

          if ( nil <> msg ) then
            msg^ := tr( 'There are too few points to specify this function.' )
          ;
          exit;
        end
      ;

      // Make sure that all of the Y values are positive.
      //-------------------------------------------------
      for i := 0 to length( points ) - 1 do
        begin
          if( 0.0 > points[i].y ) then
            begin
              result := false;

              if ( nil <> msg ) then
                msg^ := tr( 'All Y values of a piecewise PDF must be greater than or equal to 0.' )
              ;
              exit;
            end
          ;
        end
      ;

      // Verify that each x is greater than the preceding one.
      //------------------------------------------------------
      if( not( successiveXValuesOK( points ) ) ) then
        begin
          result := false;

          if ( nil <> msg ) then
            msg^ := tr( 'Each X value must be greater than the preceding one.' )
          ;
          exit;
        end
      ;

      // Calculated the cumulative area
      //-------------------------------
      setLength( cumul, length( points ) );

      cumul[0] := 0.0;

      for i := 1 to length( points ) - 1 do
        begin
          width := points[i].x - points[i-1].x;
          cumul[i] := cumul[ i - 1 ] + ( ( ( points[i-1].y + points[i].y ) / 2 ) * width );
        end
      ;

      // Standardize y coordinates if necessary
      //---------------------------------------
      total := cumul[ length( points ) - 1 ];

      if( EqualsValue <> compareValue( total, 1.0, 0.00001 ) ) then
        begin
          for i := 0 to length( points ) - 1 do
            points[i].y := points[i].y / total
          ;
        end
      ;

      // Clean up and go home
      //---------------------
      setLength( cumul, 0 );
    end
  ;
  //-----------------------------------------------------------------------------


  procedure TPdfPiecewise.calculateSlopeAndCumul();
    var
      i: integer;
      width: double;
    begin
      if( not( validate() ) ) then
        exit
      ;

      // Free old values...
      setLength( _cumul, 0 );
      setLength( _slope, 0 );

      // ...then set new values.
      setLength( _slope, pointCount - 1 );
      setLength( _cumul, pointCount );

      _cumul[0] := 0.0;

      for i := 1 to pointCount - 1 do
        begin
          width := _points[i].x - _points[i-1].x;
          _slope[i-1] := ( _points[i].y - _points[i-1].y ) / width;
          _cumul[i] := _cumul[ i - 1 ] + ( ( ( _points[i-1].y + _points[i].y ) / 2 ) * width );
        end
      ;
    end
  ;

  (*
  // AR 10/30/09 This function isn't really safe.  Make it go away.
  procedure TPdfPiecewise.addPoint( x, y: double );
    begin
      setLength( _points, length(_points)+1 );
      _points[high(_points)].x := x;
      _points[high(_points)].y := y;
      freeDiscreteValues();
    end
  ;
  *)

  procedure TPdfPiecewise.clearPoints();
    begin
      setLength( _points, 0 );
      setLength( _slope, 0 );
      setLength( _cumul, 0 );
    end
  ;


  function TPdfPiecewise.getPoint( i: integer ): RPoint;
    begin
      result := _points[i];
    end
  ;


  function TPdfPiecewise.createPointArray(): T2DPointList;
    var
      arr: T2DPointList;
      i: integer;
    begin
      arr := T2DPointList.Create();

      for i := 0 to length( _points ) - 1 do
        arr.Append( T2DPoint.create( _points[i].x, _points[i].y ) )
      ;

      result := arr;
    end
  ;

  function TPdfPiecewise.createRecordPointArray(): RPointArray;
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


  function TPdfPiecewise.getPointCount(): integer;
    begin
      result := length( _points );
    end
  ;

  
  function TPdfPiecewise.validate( err: pstring = nil ): boolean;
    begin
      result := true;

      if( 3 > length( _points ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'There are too few points to specify this function.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( ( 0.0 > min ) or ( 1.0 < max ) ) then
            begin
              if( err <> nil ) then err^ := tr( 'Proportions and probabilities must be between 0 and 1, inclusive.' );
              result := false;
            end
          ;
        end
      else if( not( standardize( _points, err ) ) ) then
        begin
          // Standardize() will check the following conditions:
          //  - All y (probability) values must be greater than or equal to 0.
          //  - Each X must be greater than the preceding one.
          //  - There must be at least 3 points specified.
          //  - The X points on the ends must have y = 0.
          //  - The sum of the Y values must be greater than 0 (otherwise, it's a flat line).
          result := false;
        end
      ;
    end
  ;


  function TPdfPiecewise.getDescr(): string;
    begin
      result := tr( 'Piecewise' );
    end
  ;


  function TPdfPiecewise.getMean(): double;
    begin
      // If it were really necessary, the mean could be estimated by simulation.
      raise exception.Create( 'Mean cannot be calculated for a piecewise distribution.' );
      result := 0.0;
    end
  ;


  function TPdfPiecewise.getHasMean(): boolean;
    begin
      result := false;
    end
  ;


  function TPdfPiecewise.getMin(): double;
    begin
      result := _points[0].x;
    end
  ;

  function TPdfPiecewise.getMax(): double;
    begin
      result := _points[ length( _points ) - 1 ].x;
    end
  ;

  
  procedure TPdfPiecewise.debug();
    var
      str: string;
      i: integer;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
        str := str + 'Piecewise PDF (id ' + intToStr( self.id ) + ') :' + endl;
      for i := low( _points ) to high( _points ) do
        str := str +  '(' + usFloatToStr(_points[i].x) + ', ' + usFloatToStr(_points[i].y) + ')' + endl
      ;

      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfPiecewise.createCopy(): TChartFunction;
    begin
      result := TPdfPiecewise.create( self );
    end
  ;


  // Translated from N. Harvey's code for piecewise PDFs, included below.
  function TPdfPiecewise.probDensity( x: double ): double;
    var
      seg: integer; // the segment in which x lies
      low, high: integer; // for a binary search
    begin
      if( ( x <= _points[0].x ) or ( x >= _points[ pointCount-1 ].x ) ) then
        begin
          result := 0.0;
          exit;
        end
      ;

      // Find the segment in which x lies (binary search)
      low := 0;
      high := length( _points );

      while( high - low > 1 ) do
        begin
          seg := trunc( (low + high) / 2 );

          if( x > _points[seg].x ) then
            low := seg
          else
            high := seg
          ;
        end
      ;

      // Linear interpolation
      x := x - _points[ low ].x;
      result := _points[ low ].y + x * _slope[ low ];
    end
  ;


  // Translated from N. Harvey's code for piecewise PDFs.
  function TPdfPiecewise.cumulativeDistr( x: double ): double;
    var
      seg: integer;  // the segment in which x lies
      low, high: integer; // for a binary search
      y: double;
    begin
      if( x <= _points[0].x ) then
        result := 0.0
      else if( x >= _points[ pointCount - 1 ].x ) then
        result := 1.0
      else
        begin
          // Find the segment in which x lies (binary search)
          low := 0;
          high := pointCount;

          while( ( high - low ) > 1 ) do
            begin
              seg := trunc( ( low + high ) / 2 );

              if( x >= _points[ seg ].x ) then
                low := seg
              else
                high := seg
              ;
            end
          ;

          x := x - _points[low].x;
          y := _points[low].y + x * _slope[low];

          result := _cumul[low] + x * (( y + _points[low].y ) / 2 );
        end
      ;
    end
  ;


  function TPdfPiecewise.invCumulative( p: double ): double;
    var
      seg: integer;
      low, high: integer;
      x, y: double;
      slope: double;

      x1, x2: double;
      maxX: double;
      nSolutions: integer;
    begin
      if( 0.0 >= p ) then
        result := _points[0].x
      else  if( 1.0 <= p ) then
        result := _points[ pointCount - 1 ].x
      else
        begin
          // Find the segment in which area lies (binary search)
          low := 0;
          high := pointCount;

          while( high - low > 1 ) do
            begin
              seg := trunc( ( low + high ) / 2 );
              if( p > _cumul[ seg ] ) then
                low := seg
              else
                high := seg
              ;
            end
          ;

          p := p - _cumul[ low ];
          y := _points[ low ].y;
          slope := _slope[ low ];

          if( EqualsValue = compareValue( slope, 0.0, 0.000001 ) ) then
            x := p / y
          else
            begin
              // This case requires solving a quadratic equation
              nSolutions := gsl_poly_solve_quadratic( (slope / 2), y, (-1 * p), @x1, @x2 );
              if( 1 = nSolutions ) then
                x := x1
              else
                begin
                  // Check which one is the solution
                  maxX := _points[ low + 1 ].x - _points[ low ].x;
                  //x = ((x1 >= 0) && (x1 <= max_x) ? x1 : x2);
                  if( (x1 >= 0.0) and (x1 <= maxX) ) then
                    x := x1
                  else
                    x := x2
                  ;
                end
              ;
            end
          ;

          result := x + _points[low].x;
        end
      ;
    end
  ;


  // Translated from N. Harvey's code for piecewise PDFs.
  function TPdfPiecewise.rand(): double;
    var
      r: double;
      firstX, lastX: double;
    begin
      firstX := _points[0].x;
      lastX := _points[ pointCount - 1 ].x;

      r := 0.0;

      while( true ) do
        begin
          r := firstX + rngRand() * ( lastX - firstX );

          if( rngRand() <= probDensity( r ) ) then
            break
          ;
        end
      ;

      result := r;
    end
  ;


  function TPdfPiecewise.ssXml( const indent: integer ): string;
    var
      i, j: integer;
      valStr: string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;

      result := result + _indent + '  <piecewise>' + endl;

      for i := low( _points ) to high( _points ) do
        begin
          valStr := usFloatToStr( _points[i].x );
          result := result + _indent + '    <value><x>' + valStr + '</x>';

          for j := ( length( valStr ) + 8 + 9 ) to 28 do
            begin
              result := result + ' ';
            end
          ;
          result := result + '  ';

          result := result + '<p>' + usFloatToStr( _points[i].y ) + '</p></value>' + endl;
        end
      ;
      result := result + _indent + '  </piecewise>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfPiecewise.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
      i: integer;
      e: pointer;
      oldStyle, newStyle: boolean;

      procedure fillValuesFromOldXml();
        var
          hadError:bool;
          I:Integer;
          Count:Integer;
          value:String;
          P:String;
          Index:Integer;
          pairs:Integer;
        begin
          Count := Sdew.GetElementCount( fnElement );

          hadError := false;

          if (( Count mod 2 ) = 0 ) then
            begin
              Index := 0;
              pairs := Count div 2;
              For I := 1 to pairs do
                begin
                  if ( 'value' = String(Sdew.GetElementNameByIndex( fnElement, Index ) ) ) then
                    begin
                      value := Sdew.GetElementContents(Sdew.GetElementByIndex( fnElement, Index ));
                      if ( 'p' = String( Sdew.GetElementNameByIndex( fnElement, Index + 1) ) ) then
                        begin
                          P := Sdew.GetElementContents(Sdew.GetElementByIndex( fnElement, Index + 1 ));
                          Index := Index + 2;
                          setLength( _points, length(_points)+1 );
                          _points[high(_points)].x := usStrToFloat(value);
                          _points[high(_points)].y := usStrToFloat(P);
                        end
                      else
                        begin
                          hadError := true;
                          break;
                        end;
                    end
                  else
                    begin
                      hadError := true;
                      break;
                    end;
                end;
            end
          else
            hadError := true
          ;


          if ( not hadError ) then
            standardize( _points )
          else
            initialize()
          ;
        end
      ;

      procedure fillValuesFromNewXml();
        var
          i: integer;
          x, p: double;
          ve, xe, pe: pointer;
          hadError: boolean;
        begin
          hadError := false;

          for i := 0 to sdew.GetElementCount( fnElement ) - 1 do
            begin
              ve := sdew.GetElementByIndex( fnElement, i );
              if( 'value' = sdew.getElementName( ve ) ) then
                begin
                  xe := sdew.GetElementByName( ve, 'x' );
                  pe := sdew.GetElementByName( ve, 'p' );

                  if( nil <> xe ) then
                    x := usStrToFloat( sdew.GetElementContents( xe ), NaN )
                  else
                    x := NaN
                  ;
                  if( nil <> pe ) then
                    p := usStrToFloat( sdew.GetElementContents( pe ), NaN )
                  else
                    p := NaN
                  ;

                  if( ( not isNaN( x ) ) and ( not isNaN( p ) ) ) then
                    begin
                      setLength( _points, length(_points)+1 );
                      _points[high(_points)].x := x;
                      _points[high(_points)].y := p;
                    end
                  else
                    begin
                      hadError := true;
                      break;
                    end
                  ;
                end
              ;
            end
          ;

          if ( not hadError ) then
            standardize( _points )
          else
            initialize()
          ;
        end
      ;
    begin // importXml
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'piecewise' );

      if( nil = fnElement ) then
        begin
          result := false;
          exit;
        end
      ;

      // Old-style XML may have the name attached to the type tag
      if( strIsEmpty( self.name ) ) then
        begin
          nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
          if( nil <> nameAttribute ) then
            name := sdew.GetElementAttribute( fnElement, 'name' )
          ;
        end
      ;

      // Determine whether this is an old-style or a new-style XML block.
      // Then call the appropriate helper function.
      oldStyle := false;
      newStyle := false;

      for i := 0 to sdew.GetElementCount( fnElement ) - 1 do
        begin
          e := sdew.GetElementByIndex( fnElement, i );
          if( 'value' = sdew.getElementName( e ) ) then
            begin
              if( nil <> sdew.GetElementByName( e, 'x' ) ) then
                newStyle := true
              else
                oldStyle := true
              ;
            end
          ;
        end
      ;

      if( newStyle and oldStyle ) then
        raise exception.Create( 'Multiple XML schema versions seem to be mixed in TPdfPiecewise.importXml()' )
      else if( ( not newStyle ) and ( not oldStyle ) ) then
        raise exception.Create( 'XML schema version cannot be determined in TPdfPiecewise.importXml()' )
      else if( newStyle ) then
        fillValuesFromNewXml()
      else if( oldStyle ) then
        fillValuesFromOldXml()
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfPiecewise.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
    var
      q, q2: string;
      i: integer;
      dict: TQueryDictionary;
    begin
      dict := createDBFieldList( db );

      if( update ) then
        begin
          q := 'DELETE FROM `inChartDetail` WHERE `chartID` = ' + intToStr( id );
          db.execute( q );
        end
      ;

      dict['chartType']   := '"Piecewise"';

      if( not update ) then
        q := writeQuery( 'inChart', QInsert, dict )
      else
        q := writeQuery( 'inChart', QUpdate, dict, 'WHERE [chartID] = ' + intToStr( id ) );
      ;

      dict.Clear();
      freeAndNil( dict );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Point PDF
//-----------------------------------------------------------------------------
  constructor TPdfPoint.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;

  constructor TPdfPoint.create( c: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setPoint( c );
      xUnits := u;
    end
  ;


  constructor TPdfPoint.create( const src: TPdfPoint );
    begin
      dbcout( '*** Calling TPdfPoint copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setPoint( src._point );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfPoint.initialize();
    begin
      setPdfType( PdfPoint );

      setPoint( NaN );

      _infLeftTail := false;
      _infRightTail := false;
    end
  ;


  destructor TPdfPoint.destroy();
    begin
      //dbcout( 'Destroying TPdfPoint ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfPoint.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._point = (Chart as TPdfPoint)._point ) then
            result := true;
        end;
    end;

  function TPdfPoint.validate( err: pstring = nil ): boolean;
    begin
      if( isNaN( _point ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'Value must be specified.' );
          result := false;
        end
      (*
      // For some applications, negative values are legitimate and should be allowed.
      else if( _point < 0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Value must be equal to or greater than 0.' );
          result := false;
        end
      *)
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( ( 0.0 > min ) or ( 1.0 < max ) ) then
            begin
              if( err <> nil ) then err^ := tr( 'Proportions and probabilities must be between 0 and 1, inclusive.' );
              result := false;
            end
          else
            result := true
          ;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfPoint.probDensity( x: double ): double;
    begin
      raise exception.create( 'TPdfPoint.probDensity(): Point PDFs do not have a probability density' );
      result := NaN;
    end
  ;


  function TPdfPoint.cumulativeDistr( x: double ): double;
    begin
      raise exception.create( 'TPdfPoint.cumulativeDistr(): Point PDFs do not have a cumulative distribution at ' + usFloatToStr( x ) );
      result := NaN;
    end
  ;


  function TPdfPoint.invCumulative( p: double ): double;
    begin
      raise exception.create( 'TPdfPoint.invCumulative(): Point PDFs do not have a inverse cumulative function' );
      result := NaN;
    end
  ;


  function TPdfPoint.rand(): double;
    begin
      result := _point;
    end
  ;


  function TPdfPoint.getDescr(): string;
    begin
      result := tr( 'Point' ) + format( ' ( %f )', [_point] );
    end
  ;

  
  function TPdfPoint.getMean(): double;
    begin
      result := value;
    end
  ;


  function TPdfPoint.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfPoint.getMin(): double;
    begin
      result := value;
    end
  ;


  function TPdfPoint.getMax(): double;
    begin
      result := value;
    end
  ;


  function TPdfPoint.getMinD(): integer;
    begin
      result := trunc( RoundDblTo( value, 0 ) );
    end
  ;


  function TPdfPoint.getMaxD(): integer;
    begin
      result := trunc( RoundDblTo( value, 0 ) );
    end
  ;


  procedure TPdfPoint.initializeDiscreteVals();
    begin
      _discreteVals := TQIntegerDoubleMap.create();
      _discreteVals.Add( minD, 1.0 );

      standardizePosDiscreteVals();
    end
  ;


  function TPdfPoint.createPointArray(): T2DPointList;
    var
      a, b, c, x, y: double;
      arr: T2DPointList;
    begin
      arr := T2DPointList.Create();

      a := mode - 0.05;
      c := mode;
      b := mode + 0.05;

      arr.append( T2DPoint.create( a, 0 ) );

      x := c;
      y := 2*(x-a)/((b-a)*(c-a));
      arr.append( T2DPoint.create( c, y ) );

      arr.append( T2DPoint.create( b, 0 ) );

      result := arr;
    end
  ;


  procedure TPdfPoint.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Point PDF (id ' + intToStr( self.id ) + ') : '
        + 'value=' + usFloatToStr( _point ) + endl
        + ' units of ' + xUnitsString() + endl
      ;

      dbcout( str, true );
    end
  ;

  function TPdfPoint.getPoint(): double; begin result := _point; end;
  procedure TPdfPoint.setPoint( val: double ); begin _point := val; freeDiscreteValues(); end;


  function TPdfPoint.createCopy(): TChartFunction;
    begin
      result := TPdfPoint.create( self );
    end
  ;


  function TPdfPoint.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <point>' + endl;
      result := result + _indent + '    ' + usFloatToStr( _point ) + endl;
      result := result + _indent + '  </point>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfPoint.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'point' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setPoint( usStrToFloat( Sdew.GetElementContents( fnElement ), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfPoint.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Point"';
      dict['mode']        := usFloatToStr( _point );

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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Triangular PDF
//-----------------------------------------------------------------------------
  constructor TPdfTriangular.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfTriangular.create( a, c, b: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setA( a );
      setC( c );
      setB( b );
      xUnits := u;
    end
  ;


  constructor TPdfTriangular.create( const src: TPdfTriangular );
    begin
      dbcout( '*** Calling TPdfTriangular copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setA( src.a );
      setC( src.c );
      setB( src.b );
      xUnits := src.xUnits;
    end
  ;

  procedure TPdfTriangular.initialize();
    begin
      setPdfType( PdfTriangular );

      setA( NaN );
      setB( NaN );
      setC( NaN );

      _infLeftTail := false;
      _infRightTail := false;
    end
  ;


  destructor TPdfTriangular.destroy();
    begin
      //dbcout( 'Destroying TPdfTriangular ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfTriangular.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._a = (Chart as TPdfTriangular)._a ) and
             ( self._b = (Chart as TPdfTriangular)._b ) and
             ( self._c = (Chart as TPdfTriangular)._c ) then
            result := true;
        end;
    end;

  procedure TPdfTriangular.setA( val: double ); begin _a := val; freeDiscreteValues(); end;
  procedure TPdfTriangular.setC( val: double ); begin _c := val; freeDiscreteValues(); end;
  procedure TPdfTriangular.setB( val: double ); begin _b := val; freeDiscreteValues(); end;
  function TPdfTriangular.getA(): double; begin Result := _a; end;
  function TPdfTriangular.getC(): double; begin Result := _c; end;
  function TPdfTriangular.getB(): double; begin Result := _b; end;


  function TPdfTriangular.getDescr(): string;
    begin
      result := tr( 'Triangular' ) + format( ' ( %f, %f, %f )', [_a, _c, _b] );
    end
  ;


  function TPdfTriangular.getMean(): double;
    begin
      result := ( a + b + c ) / 3;
    end
  ;


  function TPdfTriangular.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfTriangular.getMin(): double;
    begin
      result := _a;
    end
  ;

  function TPdfTriangular.getMax(): double;
    begin
      result := _b;
    end
  ;


  function TPdfTriangular.validate( err: pstring = nil ): boolean;
    begin
      result := true;

      if( isNan(a) or isNan(b) or isNaN(c) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( a < 0 ) or ( c < 0 ) or ( b < 0 ) then
        begin
          if( err <> nil ) then err^ := tr( 'Values must be greater than or equal to 0.' );
          result := false;
        end
      else if not( (a < c) and (c < b) ) then
        begin
          if( err <> nil ) then err^ := tr( 'The minimum value must be less than the mode, which must be less than the maximum.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( ( 0.0 > min ) or ( 1.0 < max ) ) then
            begin
              if( err <> nil ) then err^ := tr( 'Proportions and probabilities must be between 0 and 1, inclusive.' );
              result := false;
            end
          ;
        end
      ;
    end
  ;


  function TPdfTriangular.createPointArray(): T2DPointList;
    var
      arr: T2DPointList;
    begin
      arr := T2DPointList.Create();
      arr.Append( T2DPoint.create( min, probDensity( min ) ) );
      arr.Append( T2DPoint.create( mode, probDensity( mode ) ) );
      arr.Append( T2DPoint.create( max, probDensity( max ) ) );

      result := arr;
    end
  ;


  {*
    Computes the probability density p(x) at x for this triangular distribution.

    @param x
    @return the probability density p(x)
  }
  function TPdfTriangular.probDensity( x: double ): double;
    var
      val: real;
    begin
      val := NaN;

      if( (x > b) or (x < a) ) then
        raise exception.Create( 'Value out of bounds in TPdfTriangular.probDensity' )
      else if( x <= c ) then
        val := 2*(x-a)/((b-a)*(c-a))
      else if( x > c ) then
        val := 2*(b-x)/((b-a)*(b-c))
      ;

      result := val;
    end
  ;


  {*
    Computes the cumulative area <= \a x for this triangular distribution.

    @param x
    @return the cumulative area <= \a x
  }
  function TPdfTriangular.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @aphi_triangular_cdf ) then
        result := aphi_triangular_cdf( x, self.min, self.mode, self.max )
      else
        result := NaN
      ;
    end
  ;


  {*
    The inverse cumulative distribution function for this triangular distribution.

    @param p 0 <= \a p <= 1.
    @return the value at which the cumulative distribution function = \a p.
  }
  function TPdfTriangular.invCumulative( p: double ): double;
    begin
      if( nil <> @aphi_triangular_inverse_cdf ) then
        result := aphi_triangular_inverse_cdf( p, self.min, self.mode, self.max )
      else
        result := NaN
      ;
    end
  ;


  function TPdfTriangular.rand(): double;
    begin
      if( nil <> @aphi_triangular_rand ) then
        result := aphi_triangular_rand( rngPtr(), self.min, self.mode, self.max )
      else
        result := NaN
      ;
    end
  ;


  procedure TPdfTriangular.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Triangular PDF (id ' + intToStr( self.id ) + ') :' + endl + 'a=' + usFloatToStr(a) + ', b=' + usFloatToStr(b) + ', c=' + usFloatToStr(c);
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfTriangular.createCopy(): TChartFunction;
    begin
      result := TPdfTriangular.create( self );
    end
  ;


  function TPdfTriangular.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <triangular>' + endl;
      result := result + _indent + '    <a>' + usFloatToStr( a ) + '</a>' + endl;
      result := result + _indent + '    <c>' + usFloatToStr( c ) + '</c>' + endl;
      result := result + _indent + '    <b>' + usFloatToStr( b ) + '</b>' + endl;
      result := result + _indent + '  </triangular>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfTriangular.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'triangular' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setA( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'a' )), NaN ) );
          setB( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'b' )), NaN ) );
          setC( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'c' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfTriangular.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Triangular"';
      dict['min']         := usFloatToStr( a );
      dict['mode']        := usFloatToStr( c );
      dict['max']         := usFloatToStr( b );

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

      // KEEP THIS CODE: it might be useful when I write the chart conversion function.
      (*
      mode: real;
      mode := valueAt( c );

      q := 'INSERT INTO [inChartDetail] ([chartID], [pointOrder], [x], [y])';
      q2 := q + ' VALUES (' + usFloatToStr( id ) + ', 0, ' + usFloatToStr( a ) + ', 0)';
      db.execute( q2 );

      q2 := q + ' VALUES (' + usFloatToStr( id ) + ', 1, ' + usFloatToStr( c ) + ', ' + usFloatToStr( mode ) + ')';
      db.execute( q2 );

      q2 := q + ' VALUES (' + usFloatToStr( id ) + ', 2, ' + usFloatToStr( b ) + ', 0)';
      db.execute( q2 );
      *)

    end
  ;
  {$ENDIF}
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Uniform PDF
//-----------------------------------------------------------------------------
  constructor TPdfUniform.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfUniform.create( a, b: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setA( a );
      setB( b );
      xUnits := u;
    end
  ;


  constructor TPdfUniform.create( const src: TPdfUniform );
    begin
      dbcout( '*** Calling TPdfUniform copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setA( src.a );
      setB( src.b );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfUniform.initialize();
    begin
      setPdfType( PdfUniform );

      setA( NaN );
      setB( NaN );

      _infLeftTail := false;
      _infRightTail := false;
    end
  ;


  destructor TPdfUniform.destroy();
    begin
      //dbcout( 'Destroying TPdfUniform ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfUniform.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._a = (Chart as TPdfUniform)._a ) and
             ( self._b = (Chart as TPdfUniform)._b ) then
            result := true;
        end;
    end;

  function TPdfUniform.getDescr(): string;
    begin
      result := tr( 'Uniform' ) + format( ' ( %f, %f )', [_a, _b] );
    end
  ;


  function TPdfUniform.getMean(): double;
    begin
      result := ( min + max ) / 2;
    end
  ;

  
  function TPdfUniform.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfUniform.getMin(): double;
    begin
      result := _a;
    end
  ;

  function TPdfUniform.getMax(): double;
    begin
      result := _b;
    end
  ;


  function TPdfUniform.createPointArray(): T2DPointList;
    var
      arr: T2DPointList;
    begin
      arr := T2DPointList.Create();
      arr.Append( T2DPoint.create( min, probDensity(min) ) );
      arr.Append( T2DPoint.create( max, probDensity(max) ) );

      result := arr;
    end
  ;

  procedure TPdfUniform.setA( val: double ); begin _a := val; freeDiscreteValues(); end;
  procedure TPdfUniform.setB( val: double ); begin _b := val; freeDiscreteValues(); end;
  function TPdfUniform.getA(): double; begin Result := _a; end;
  function TPdfUniform.getB(): double; begin Result := _b; end;

  function TPdfUniform.validate( err: pstring = nil ): boolean;
    begin
      result := true;

      if( isNaN( a ) or isNaN( b ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( ( a < 0 ) or ( b < 0 ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'Values must be greater than or equal to 0.' );
          result := false;
        end
      else if not( (a < b) ) then
        begin
          if( err <> nil ) then err^ := tr( 'The minimum value must be less than the maximum.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( ( 0.0 > min ) or ( 1.0 < max ) ) then
            begin
              if( err <> nil ) then err^ := tr( 'Proportions and probabilities must be between 0 and 1, inclusive.' );
              result := false;
            end
          ;
        end
      ;
    end
  ;


  function TPdfUniform.probDensity( x: double ): double;
    begin
      if( x < min ) or ( x > max ) then
        begin
          raise exception.Create( 'Value out of bounds in TPdfUniform.probDensity' );
          result := NaN;
        end
      else
        result := 1 / (max - min )
      ;
    end
  ;


  function TPdfUniform.cumulativeDistr( x: double ): double;
    begin
      if( nil <> @gslCdfUniformP ) then
        result := gslCdfUniformP( x, self.min, self.max )
      else
        result := NaN
      ;
    end
  ;


  function TPdfUniform.invCumulative( p: double ): double;
    begin
      if( nil <> @gslCdfUniformPinv ) then
        result := gslCdfUniformPInv( p, self.min, self.max )
      else
        result := NaN
      ;
    end
  ;


  function TPdfUniform.rand(): double;
    begin
      if( nil <> @gslRanUniform ) then
        result := gslRanUniform( gslRngPtr(), self.min, self.max )
      else
        result := NaN
      ;
    end
  ;


  class function TPdfUniform.rand( const min, max: double ): double;
    begin
      if( nil <> @gslRanUniform ) then
        result := gslRanUniform( gslRngPtr(), min, max )
      else
        result := NaN
      ;
    end
  ;


  procedure TPdfUniform.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Uniform PDF (id ' + intToStr( self.id ) + ') :' + endl + 'a=' + usFloatToStr(a) + ', b=' + usFloatToStr(b);
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfUniform.createCopy(): TChartFunction;
    begin
      result := TPdfUniform.create( self );
    end
  ;


  function TPdfUniform.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <uniform>' + endl;
      result := result + _indent + '    <a>' + usFloatToStr( a ) + '</a>' + endl;
      result := result + _indent + '    <b>' + usFloatToStr( b ) + '</b>' + endl;
      result := result + _indent + '  </uniform>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;
  

  function TPdfUniform.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'uniform' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setA( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'a' )), NaN ) );
          setB( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'b' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfUniform.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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
      dict['chartType']   := '"Uniform"';
      dict['min']         := usFloatToStr( a );
      dict['max']         := usFloatToStr( b );

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

      // KEEP THIS CODE: it might be useful when I write the chart conversion function.
      (*
      mode: real;
      mode := valueAt( c );

      q := 'INSERT INTO [inChartDetail] ([chartID], [pointOrder], [x], [y])';
      q2 := q + ' VALUES (' + usFloatToStr( id ) + ', 0, ' + usFloatToStr( a ) + ', 0)';
      db.execute( q2 );

      q2 := q + ' VALUES (' + usFloatToStr( id ) + ', 1, ' + usFloatToStr( c ) + ', ' + usFloatToStr( mode ) + ')';
      db.execute( q2 );

      q2 := q + ' VALUES (' + usFloatToStr( id ) + ', 2, ' + usFloatToStr( b ) + ', 0)';
      db.execute( q2 );
      *)

    end
  ;
  {$ENDIF}
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Weibull PDF
//-----------------------------------------------------------------------------
  constructor TPdfWeibull.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;

  constructor TPdfWeibull.create( alpha, beta: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setAlpha( alpha );
      setBeta( beta );
      xUnits := u;
    end
  ;


  constructor TPdfWeibull.create( const src: TPdfWeibull );
    begin
      dbcout( '*** Calling TPdfWeibull copy constructor.', DBSHOWMSG );
      inherited create( src );
      // initialize() is called by the inherited constructor

      setAlpha( src.alpha );
      setBeta( src.beta );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfWeibull.initialize();
    begin
      setPdfType( PdfWeibull );

      setAlpha( NaN );
      setBeta( NaN );

      _infLeftTail := false;
      _infRightTail := true;
    end
  ;


  destructor TPdfWeibull.destroy();
    begin
      //dbcout( 'Destroying TPdfWeibull ' + self.name, DBSHOWMSG );
      inherited destroy();
    end
  ;

  function TPdfWeibull.compare( Chart:TChartFunction ):boolean;
    begin
      result := false;
      if ( inherited compare( Chart ) ) then
        begin
          if ( self._alpha = (Chart as TPdfWeibull)._alpha ) and
             ( self._beta = (Chart as TPdfWeibull)._beta ) then
            result := true;
        end;
    end;

  procedure TPdfWeibull.setAlpha( val: double ); begin _alpha := val; freeDiscreteValues(); end;
  procedure TPdfWeibull.setBeta( val: double ); begin _beta := val; freeDiscreteValues(); end;
  function TPdfWeibull.getAlpha(): double; begin Result := _alpha; end;
  function TPdfWeibull.getBeta(): double; begin Result := _beta; end;


  function TPdfWeibull.validate( err: pstring = nil ): boolean;
    begin
      if( isNan( alpha ) or isNan( beta ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'All values must be specified.' );
          result := false;
        end
      else if( beta <= 0.1 ) then
        begin
          if( err <> nil ) then err^ := tr( 'The specified beta value is too small for practical calculation.' );
          result := false;
        end
      else if( alpha <= 0.1 ) then
        begin
          if( err <> nil ) then err^ := tr( 'The specified alpha value is too small for practical calculation.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Unbounded probability density functions cannot be used to represent proportions or probabilities.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;

  function TPdfWeibull.probDensity( x: double ): double;
    begin
      // NOTE: GSL functions accept parameters in
      // opposite order described by Vose and @Risk.
      if( nil <> @gslRanWeibullPdf ) then
        result := gslRanWeibullPdf( x, self.beta, self.alpha )
      else
        result := NaN
      ;
    end
  ;


  function TPdfWeibull.cumulativeDistr( x: double ): double;
    begin
      // NOTE: GSL functions accept parameters in
      // opposite order described by Vose and @Risk.
      if( nil <> @gslCdfWeibullP ) then
        result := gslCdfWeibullP( x, self.beta, self.alpha )
      else
        result := NaN
      ;
    end
  ;


  function TPdfWeibull.invCumulative( p: double ): double;
    begin
      // NOTE: GSL functions accept parameters in
      // opposite order described by Vose and @Risk.
      if( nil <> @gslCdfWeibullPinv ) then
        result := gslCdfWeibullPInv( p, self.beta, self.alpha )
      else
        result := NaN
      ;
    end
  ;


  function TPdfWeibull.rand(): double;
    begin
      // NOTE: GSL functions accept parameters in
      // opposite order described by Vose and @Risk.
      if( nil <> @gslRanWeibull ) then
        result := gslRanWeibull( gslRngPtr(), self.beta, self.alpha )
      else
        result := NaN
      ;
    end
  ;



  function TPdfWeibull.getDescr(): string;
    begin
      result := tr( 'Weibull' ) + format( ' ( %f, %f )', [_alpha, _beta] );
    end
  ;


  function TPdfWeibull.getMean(): double;
    begin
      result := ( beta / alpha ) * gamma( 1/ alpha );
    end
  ;


  function TPdfWeibull.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfWeibull.getMin(): double;
    begin
      result := 0.0;
    end
  ;

  function TPdfWeibull.getMax(): double;
    begin
      result := infinity;
    end
  ;


  procedure TPdfWeibull.debug();
    var
      str: string;
    begin
      // NOTE: GSL functions accept parameters in
      // opposite order described by Vose and @Risk.
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str + 'Weibull PDF (id ' + intToStr( self.id ) + ') :' + endl + 'alpha=' + usFloatToStr(alpha) + ', beta=' + usFloatToStr(beta);
      str := str + endl + '(Remember: alpha and beta are reverse in the interface relative to the XML.)';
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  function TPdfWeibull.createCopy(): TChartFunction;
    begin
      result := TPdfWeibull.create( self );
    end
  ;


  function TPdfWeibull.ssXml( const indent: integer ): string;
    begin
      // NOTE: GSL functions accept parameters in opposite order described by Vose
      // and @Risk, which are followed by NAADSM.  It may look like the GSL function
      // is called with the parameters in the wrong order, but it's actually OK.
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <weibull>' + endl;
      result := result + _indent + '    <alpha>' + usFloatToStr( alpha ) + '</alpha>' + endl;
      result := result + _indent + '    <beta>' + usFloatToStr( beta ) + '</beta>' + endl;
      result := result + _indent + '  </weibull>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfWeibull.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'weibull' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setAlpha( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'alpha' )), NaN ) );
          setBeta( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'beta' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfWeibull.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Weibull"';
      dict['alpha']       := usFloatToStr( alpha );
      dict['beta']        := usFloatToStr( beta );
      
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
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Base class for all discrete PDF types
//-----------------------------------------------------------------------------
  function TMassHistogramValues.at( const idx: integer ): TMassHistogramValue;
    begin
      result := inherited at( idx ) as TMassHistogramValue;
    end
  ;

  constructor TPdfDiscrete.create();
    begin
      inherited create();
      _histValues := nil;
      _cumulHistValues := nil;
    end
  ;

  destructor TPdfDiscrete.destroy();
    begin
      freeHistogramValues();
      inherited destroy();
    end
  ;


  procedure TPdfDiscrete.initializeDiscreteVals();
    var
      i: integer;
    begin
      _discreteVals := TQIntegerDoubleMap.create();

      for i := 0 to massHistogramValues.Count - 1 do
        _discreteVals.Add( trunc( massHistogramValues.at(i).x ), massHistogramValues.at(i).y )
      ;

      standardizePosDiscreteVals();
    end
  ;

  
  procedure TPdfDiscrete.freeHistogramValues();
    begin
      freeAndNil( _histValues );
      freeAndNil( _cumulHistValues );

      freeDiscreteValues();
    end
  ;

  function TPdfDiscrete.getIsDiscrete(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfDiscrete.getCumulHistValues(): TMassHistogramValues;
    begin
      if( nil = _cumulHistValues ) then
        getHistValues()
      ;

      result := _cumulHistValues;
    end
  ;


  procedure TPdfDiscrete.buildCumulHistValues();
    var
      i: integer;
      val: TMassHistogramValue;
    begin
      freeAndNil( _cumulHistValues );
      _cumulHistValues := TMassHistogramValues.create();

      val := TMassHistogramValue.create( _histValues.at( 0 ).x, _histValues.at( 0 ).y );
      _cumulHistValues.append( val );

      for i := 1 to _histValues.Count - 1 do
        begin
          val := TMassHistogramValue.create( _histValues.at(i).x, ( _cumulHistValues.at(i-1).y + _histValues.at(i).y ) );
          _cumulHistValues.append( val );
        end
      ;

    end
  ;

  function TPdfDiscrete.cumulativeDistr( k: integer; raiseExceptionOnFailure: boolean = true ): double;
    var
      i: integer;
    begin
      result := 0.0;

      if( k < cumulHistogramValues.first().x ) then
        result := 0.0
      else if( k > cumulHistogramValues.last().x ) then
        begin
          if( infiniteRightTail ) then
            result := 1.0
          else
            begin
              if( raiseExceptionOnFailure ) then
                raise exception.Create( 'Invalid k in TPdfDiscrete.cumulativeDistr()' )
              ;
              result := NaN;
            end
          ;
        end
      else
        begin
          for i := 0 to cumulHistogramValues.Count - 1 do
            begin
              if( k = cumulHistogramValues.at( i ).x ) then
                begin
                  result := cumulHistogramValues.at( i ).y;
                  break;
                end
              ;
            end
          ;
        end
      ;
    end
  ;


  (*
    Marian Talbert suggested that the invCumulative function should find matches based on
    "in-between" values, not based on "exact" (within epsilon) "greater-than" matches,
    as was done previously (see old function below).

    Here we start, not at the end of the cuml_hist list, but at the second to last item,
    and move toward the first.  If prob is greater than the y value
    of the item, then the prob falls between the current item and the next item in the list
    (this is why we start with the second to the last).  Thus, we assign the x value of the
    next item.

    If we continued to use the old implementation, we would almost never
    reach the last item, but when we did, we would assign a value larger than its x.
  *)
  function TPdfDiscrete.invCumulative( prob: double ): integer;
    var
      i: integer;
    begin
      if( ( 0.0 > prob ) or ( 1.0 < prob ) ) then
        begin
          raise exception.create( 'Invalid probability (' + usFloatToStr( prob ) + ') in TPdfDiscrete.invCumulative()' );
          result := 0;
        end
      else
        begin
          result := 0;

          if ( prob <= cumulHistogramValues.at( 0 ).y ) then
            begin
              result := trunc( RoundDblTo( cumulHistogramValues.at( 0 ).x, 0 ) );
            end
          else          
            for i := cumulHistogramValues.count - 2 downto 0 do
              begin
                if( ( prob > cumulHistogramValues.at( i ).y ) ) then
                  begin
                    result := trunc( RoundDblTo( cumulHistogramValues.at( i + 1 ).x, 0 ) );
                    //dbcout2( 'Found! ' + intToStr( result ) + ', ' + usFloatToStr( cumulHistogramValues.at( i + 1 ).x ) );
                    break;
                  end
                ;
              end
            ;
        end
      ;

      //dbcout2( '+++++ ' +  intToStr( result ) + ', ' + usFloatToStr( prob ) );
    end
  ;


  (*
  // Note: This function has been replaced by the more accurate one above.
  // Keep this code (for the time being) for reference.
  //
  function TPdfDiscrete.invCumulative( prob: double ): integer;
    var
      i: integer;
      test: TValueRelationship;
    begin
      if( ( 0.0 > prob ) or ( 1.0 < prob ) ) then
        begin
          raise exception.create( 'Invalid probability (' + usFloatToStr( prob ) + ') in TPdfDiscrete.invCumulative()' );
          result := 0;
        end
      else
        begin
          result := 0;
          
          for i := cumulHistogramValues.count - 1 downto 0 do
            begin
              test := compareValue( prob, cumulHistogramValues.at( i ).y, 0.0001 );

              if( ( EqualsValue = test ) or ( GreaterThanValue = test ) ) then
                begin
                  result := trunc( RoundDblTo( cumulHistogramValues.at( i ).x, 0 ) ) + 1;
                  //dbcout2( 'Found! ' + intToStr( result ) + ', ' + usFloatToStr( cumulHistogramValues.at( i ).x ) );
                  break;
                end
              ;
            end
          ;
        end
      ;
    end
  ;
  *)

  function TPdfDiscrete.randInvCumulative(): integer;
    begin
      result := invCumulative( rngRand() );
    end
  ;


  function TPdfDiscrete.asCsv(): string;
    var
      i: integer;
      points: T2DPointList;
    begin
      result := '';

      if( self.validate() ) then
        begin
          points := massHistogramValues;

          // Set the first and last Y values 0.
          // This will make conversion to piecewise easier.
          points[0].y := 0.0;
          points[points.Count - 1].y := 0.0;

          try
            result := 'x' + csvListSep() + ' y' + endl;

            for i := 0 to points.Count - 1 do
              begin
                result :=
                  result
                  + csvFloatToStr( points[i].x )
                  + csvListSep() + ' '
                  + csvFloatToStr( points[i].y )
                  + endl
                ;
              end
            ;
          except
            raise exception.Create( 'An exception occurred in TPdfDiscrete.asCsv()' );
            result := '';
          end;
        end
      ;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Bernoulli PDF
//-----------------------------------------------------------------------------
  constructor TPdfBernoulli.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfBernoulli.create( p: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setP( p );
      xUnits := u;
    end
  ;


  constructor TPdfBernoulli.create( const src: TPdfBernoulli );
    begin
      inherited create( src );
      // initialize() is called by the inherited constructor

      setP( src.p );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfBernoulli.initialize();
    begin
      setPdfType( PdfBernoulli );

      setP( NaN );

      _infLeftTail := false;
      _infRightTail := false;

      freeHistogramValues();
    end
  ;


  destructor TPdfBernoulli.destroy();
    begin
      inherited destroy();
    end
  ;


  procedure TPdfBernoulli.setParams( p: double; u: TChartUnitType = UUnknown );
    begin
      setP( p );
      xUnits := u;
      freeHistogramValues();
    end
  ;

  
  function TPdfBernoulli.validate( err: pstring = nil ): boolean;
    begin
      if( isNan( p ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'Probability of success must be specified.' );
          result := false;
        end
      else if( ( 0.0 > p ) or ( 1.0 < p ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Probability of success must be between 0 and 1, inclusive.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfBernoulli.compare( chart: TChartFunction ):boolean;
    begin
      result := false;

      if( inherited compare( chart ) ) then
        begin
          if( sameValue( self._p, ( chart as TPdfBernoulli ).p ) ) then
            result := true
          ;
        end
      ;
    end
  ;


  function TPdfBernoulli.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <bernoulli>' + endl;
      result := result + _indent + '    <p>' + usFloatToStr( p ) + '</p>' + endl;
      result := result + _indent + '  </bernoulli>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfBernoulli.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'bernoulli' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          _p := usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'p' )), NaN );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfBernoulli.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Bernoulli"';
      dict['p']           := usFloatToStr( p );

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


  function TPdfBernoulli.createCopy(): TChartFunction;
    begin
      result := TPdfBernoulli.create( self );
    end
  ;


  function TPdfBernoulli.probDensity( k: integer ): double;
    begin
      case k of
        0: result := 1 - p;
        1: result := p;
        else
          begin
            raise exception.create( 'k is out of bounds in TPdfBernoulli.probDensity()' );
            result := NaN;
          end
      end;
    end
  ;


  function TPdfBernoulli.rand(): double;
    begin
      if( p < TPdfUniform.rand( 0.0, 1.0 ) ) then
        result := 0
      else
        result := 1
      ;
    end
  ;


  class function TPdfBernoulli.rand( const p: double ): double;
    begin
      if( p < TPdfUniform.rand( 0.0, 1.0 ) ) then
        result := 0
      else
        result := 1
      ;
    end
  ;


  function TPdfBernoulli.getHistValues(): TMassHistogramValues;
    begin
      if( nil = _histValues ) then
        begin
          _histValues := TMassHistogramValues.create();
          _histValues.append( TMassHistogramValue.create( 0, 1.0 - p ) );
          _histValues.append( TMassHistogramValue.create( 1, p ) );
        end
      ;

      buildCumulHistValues();

      result := _histValues;
    end
  ;

  
  procedure TPdfBernoulli.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str
        + 'Bernoulli PDF (id ' + intToStr( self.id ) + ') :' + endl
        + ', p=' + usFloatToStr( p )
      ;
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  procedure TPdfBernoulli.setP( val: double ); begin _p := val; freeHistogramValues(); end;


  function TPdfBernoulli.getDescr(): string;
    begin
      result := tr( 'Bernoulli' ) + format( ' ( %f )', [ p ] );
    end
  ;


  function TPdfBernoulli.getMean(): double;
    begin
      result := p;
    end
  ;


  function TPdfBernoulli.getHasMean(): boolean;
    begin
      result := true;
    end
  ;

  function TPdfBernoulli.getMin(): double;
    begin
      if( sameValue( 1.0, p ) ) then
        result := 1
      else
        result := 0
      ;
    end
  ;

  function TPdfBernoulli.getMax(): double;
    begin
      if( sameValue( 0.0, p ) ) then
        result := 0
      else
        result := 1
      ;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Binomial PDF
//-----------------------------------------------------------------------------
  constructor TPdfBinomial.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfBinomial.create(n: integer; p: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setN( n );
      setP( p );
      xUnits := u;
    end
  ;


  constructor TPdfBinomial.create( const src: TPdfBinomial );
    begin
      inherited create( src );
      // initialize() is called by the inherited constructor

      setN( src.n );
      setP( src.p );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfBinomial.initialize();
    begin
      setPdfType( PdfBinomial );

      setN( -1 );
      setP( NaN );

      _infLeftTail := false;
      _infRightTail := false;

      freeHistogramValues();
    end
  ;


  destructor TPdfBinomial.destroy();
    begin
      inherited destroy();
    end
  ;


  procedure TPdfBinomial.setParams( n: integer; p: double; u: TChartUnitType = UUnknown );
    begin
      setN( n );
      setP( p );
      xUnits := u;
      freeHistogramValues();
    end
  ;

  
  function TPdfBinomial.validate( err: pstring = nil ): boolean;
    begin
      if( n < 1 ) then
        begin
          if( nil <> err ) then err^ := tr( 'Number of trials must be greater than 0.' );
          result := false;
        end
      else if( isNan( p ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'Probability of success must be specified.' );
          result := false;
        end
      else if( ( 0.0 > p ) or ( 1.0 < p ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Probability of success must be between 0 and 1, inclusive.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'This probability density function cannot be used to represent a proportion or a probability.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfBinomial.compare( chart: TChartFunction ):boolean;
    begin
      result := false;

      if( inherited compare( chart ) ) then
        begin
          if
            ( self._n = ( chart as TPdfBinomial ).n )
          and
            ( self._p = ( chart as TPdfBinomial ).p )
          then
            result := true
          ;
        end
      ;
    end
  ;


  function TPdfBinomial.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <binomial>' + endl;
      result := result + _indent + '    <n>' + intToStr( n ) + '</n>' + endl;
      result := result + _indent + '    <p>' + usFloatToStr( p ) + '</p>' + endl;
      result := result + _indent + '  </binomial>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfBinomial.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'binomial' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          _n := myStrToInt( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'n' )), -1 );
          _p := usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'p' )), NaN );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfBinomial.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Binomial"';
      dict['n']           := intToStr( n );
      dict['p']           := usFloatToStr( p );

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


  function TPdfBinomial.createCopy(): TChartFunction;
    begin
      result := TPdfBinomial.create( self );
    end
  ;


  function TPdfBinomial.probDensity( k: integer ): double;
    begin
      if( 1.0 = p ) then
        begin
          if( k = n ) then
            result := 1
          else
            result := 0
          ;
        end
      else if( nil <> @gslRanBinomialPdf ) then
        result := gslRanBinomialPdf( k, self.p, self.n )
      else
        result := NaN
      ;
    end
  ;


  function TPdfBinomial.rand(): double;
    var
      r: integer;
    begin
      if( 0 = self.n ) then
        result := 0
      else if( nil <> @gslRanBinomial ) then
        begin
          r := gslRanBinomial( gslRngPtr(), self.p, self.n );
          result := r;
        end
      else
        result := NaN
      ;
    end
  ;


  class function TPdfBinomial.rand( const n: integer; p: double ): double;
    var
      r: integer;
    begin
      if( 0 = n ) then
        result := 0
      else if( nil <> @gslRanBinomial ) then
        begin
          r := gslRanBinomial( gslRngPtr(), p, n );
          result := r;
        end
      else
        result := NaN
      ;
    end
  ;


  (*
    This function is based on C++ code written by Agner Fog and released under the GPL.
    See http://www.agner.org/random/
  *)
  class procedure TPdfBinomial.multinomialRand( n: integer; const p: TQDoubleVector; output: TQIntegerVector );
    var
      sum: double;
      s: double;
      i: integer;
      c: integer;
      x: integer;
    begin
      c := p.count;

      output.resize( c );

      if( 0 = c ) then
        exit
      ;

      // Compute sum of probabilities.  Do a little error checking along the way.
      sum := 0.0;
      for i := 0 to c - 1 do
        begin
          if( 0.0 > p[i] ) then
            begin
              raise exception.create( 'Probabilities cannot be negative in TPdfBinomial.multinomialRand()' );
              exit;
            end
          ;
          sum := sum + p[i];
        end
      ;

      if( 0.0 = sum ) then
        begin
          raise exception.create( 'Sum of probabilities cannot be zero in TPdfBinomial.multinomialRand()' );
          exit;
        end
      ;

      for i := 0 to ( c - 1 ) - 1 do
        begin
          // Generate output by calling binomial (c-1) times
          s := p[i];
          if( sum <= s ) then
            begin
              // This fixes two problems:
              // 1. prevent division by 0 when sum = 0
              // 2. prevent s/sum getting bigger than 1 in case of rounding errors
              x := n;
            end
          else
            x := trunc( TPdfBinomial.rand( n, s/sum ) )
          ;

          n := n - x;
          sum := sum - s;
          output[ i ] := x;
        end
      ;
      // Get the last one
      output[ c - 1 ] := n;

      (*
         // Original C code:
         for (i=0; i<colors-1; i++) {
            // generate output by calling binomial (colors-1) times
            s = source[i];
            if (sum <= s) {
               // this fixes two problems:
               // 1. prevent division by 0 when sum = 0
               // 2. prevent s/sum getting bigger than 1 in case of rounding errors
               x = n;
            }
            else {
               x = Binomial(n, s/sum);
            }
            n -= x; sum -= s;
            destination[i] = x;
         }
         // get the last one
         destination[i] = n;
      *)
    end
  ;


  function TPdfBinomial.getHistValues(): TMassHistogramValues;
    var
      i: integer;
      p: double;
    begin
      if( nil = _histValues ) then
        begin
          _histValues := TMassHistogramValues.create();

          // Start with the number of trials and work backwards, calculating
          // the probability of each outcome until a probability of 0 is reached.
          for i := n downto 0 do
            begin
              p := probDensity( i );
              if( 0 <= p ) then
                _histValues.prepend( TMassHistogramValue.create( i, p ) )
              else
                break
              ;
            end
          ;
        end
      ;

      buildCumulHistValues();

      result := _histValues;
    end
  ;

  
  procedure TPdfBinomial.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str
        + 'Binomial PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 'n=' + intToStr( n )
        + ', p=' + usFloatToStr( p )
      ;
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  procedure TPdfBinomial.setN( val: integer ); begin _n := val; freeHistogramValues(); end;
  procedure TPdfBinomial.setP( val: double ); begin _p := val; freeHistogramValues(); end;


  function TPdfBinomial.getDescr(): string;
    begin
      result := tr( 'Binomial' ) + format( ' ( %d, %f )', [n, p] );
    end
  ;


  function TPdfBinomial.getMean(): double;
    begin
      result := n * p;
    end
  ;


  function TPdfBinomial.getHasMean(): boolean;
    begin
      result := true;
    end
  ;

  function TPdfBinomial.getMin(): double;
    begin
      result := 0;
    end
  ;

  function TPdfBinomial.getMax(): double;
    begin
      result := n;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Discrete uniform PDF
//-----------------------------------------------------------------------------
  constructor TPdfDiscreteUniform.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfDiscreteUniform.create( min, max: integer; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setMin( min );
      setMax( max );
      xUnits := u;
    end
  ;


  constructor TPdfDiscreteUniform.create( const src: TPdfDiscreteUniform );
    begin
      inherited create( src );
      // initialize() is called by the inherited constructor

      setMin( trunc( src.min ) );
      setMax( trunc( src.max ) );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfDiscreteUniform.initialize();
    begin
      setPdfType( PdfDiscreteUniform );

      setMin( -1 );
      setMax( -1 );

      _infLeftTail := false;
      _infRightTail := false;
    end
  ;


  destructor TPdfDiscreteUniform.destroy();
    begin
      inherited destroy();
    end
  ;


  function TPdfDiscreteUniform.validate( err: pstring = nil ): boolean;
    begin
      result := true;

      if( min >= max ) then
        begin
          if( nil <> err ) then err^ := tr( 'The maximum must be greater than the minimum.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( ( 0.0 > min ) or ( 1.0 < max ) ) then
            begin
              if( err <> nil ) then err^ := tr( 'Proportions and probabililities must be between 0 and 1, inclusive.' );
              result := false;
            end
          ;
        end
      ;
    end
  ;


  function TPdfDiscreteUniform.compare( chart: TChartFunction ):boolean;
    begin
      result := false;

      if( inherited compare( chart ) ) then
        begin
          if
            ( self.min = ( chart as TPdfDiscreteUniform ).min )
          and
            ( self.max = ( chart as TPdfDiscreteUniform ).max )
          then
            result := true
          ;
        end
      ;
    end
  ;


  function TPdfDiscreteUniform.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <discrete-uniform>' + endl;
      result := result + _indent + '    <min>' + intToStr( trunc( min ) ) + '</min>' + endl;
      result := result + _indent + '    <max>' + intToStr( trunc( max ) ) + '</max>' + endl;
      result := result + _indent + '  </discrete-uniform>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfDiscreteUniform.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
      min, max: double;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'discrete-uniform' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          min := usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'min' ) ), NaN );
          max := usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'max' ) ), NaN );

          if( ( NaN <> min ) and ( NaN <> max ) ) then
            begin
              if( ( SameValue( roundDblTo( min, 0 ), min ) ) and ( SameValue( roundDblTo( max, 0 ), max ) ) ) then
                begin
                  _min := roundDbl( min );
                  _max := roundDbl( max );
                end
              else // One or both values aren't integers
                initialize()
              ;
            end
          else // One or both values aren't numbers
            initialize()
          ;
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfDiscreteUniform.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Discrete-uniform"';
      dict['dMin']          := intToStr( trunc( min ) );
      dict['dMax']         := intToStr( trunc( max ) );

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


  function TPdfDiscreteUniform.createCopy(): TChartFunction;
    begin
      result := TPdfDiscreteUniform.create( self );
    end
  ;


  function TPdfDiscreteUniform.probDensity( k: integer ): double;
    var
      range: integer;
    begin
      if( ( k > max ) or ( k < min ) ) then
        result := 0.0
      else
        begin
          range := abs( trunc( max ) - trunc( min ) ) + 1;
          result := 1/range;
        end
      ;
    end
  ;


  function TPdfDiscreteUniform.rand(): double;
    begin
      result := rngRandInt( trunc( min ), trunc( max ) + 1 );
    end
  ;


  function TPdfDiscreteUniform.getHistValues(): TMassHistogramValues;
    var
      range: integer;
      i: integer;
    begin
      if( nil = _histValues ) then
        begin
          _histValues := TMassHistogramValues.create();
          
          range := abs( trunc( max ) - trunc( min ) ) + 1;

          for i := trunc( min ) to trunc( max ) do
            _histValues.append( TMassHistogramValue.create( i, 1/range ) )
          ;
        end
      ;

      buildCumulHistValues();

      result := _histValues;
    end
  ;

  
  procedure TPdfDiscreteUniform.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str
        + 'Discrete uniform PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 'min=' + intToStr( trunc( min ) )
        + ', max=' + intToStr( trunc( max ) )
      ;
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;

  
  procedure TPdfDiscreteUniform.setMin( val: integer );
    begin
      _min := val;
      freeHistogramValues();
    end
  ;


  procedure TPdfDiscreteUniform.setMax( val: integer );
    begin
      _max := val;
      freeHistogramValues();
    end
  ;


  function TPdfDiscreteUniform.getDescr(): string;
    begin
      result := tr( 'Discrete uniform' ) + format( ' ( %d, %d )', [ trunc( min ), trunc( max ) ] );
    end
  ;


  function TPdfDiscreteUniform.getMean(): double;
    var
      i: integer;
      n, sum: integer;
    begin
      n := 0;
      sum := 0;
      for i := trunc( min ) to trunc( max ) do
        begin
          inc( n );
          inc( sum, i );
        end
      ;

      result := sum / n;
    end
  ;


  function TPdfDiscreteUniform.getHasMean(): boolean;
    begin
      result := true;
    end
  ;

  function TPdfDiscreteUniform.getMin(): double;
    begin
      result := _min;
    end
  ;

  function TPdfDiscreteUniform.getMax(): double;
    begin
      result := _max;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// TPdfHypergeometric PDF
//-----------------------------------------------------------------------------
  constructor TPdfHypergeometric.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfHypergeometric.create( n, d, m: integer; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setM( m );
      setD( d );
      setN( n );
      xUnits := u;
    end
  ;


  constructor TPdfHypergeometric.create( const src: TPdfHypergeometric );
    begin
      inherited create( src );
      // initialize() is called by the inherited constructor

      setM( src.m );
      setD( src.d );
      setN( src.n );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfHypergeometric.initialize();
    begin
      setPdfType( PdfHypergeometric );

      setM( -1 );
      setD( -1 );
      setN( -1 );

      _infLeftTail := false;
      _infRightTail := false;
    end
  ;


  destructor TPdfHypergeometric.destroy();
    begin
      inherited destroy();
    end
  ;


  function TPdfHypergeometric.validate( err: pstring = nil ): boolean;
    begin

      // All values must be greater than or equal to 0
      if( ( 0 > n ) or ( 0 > m ) or ( 0 > d ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'All parameters must be greater than or equal to 0.' );
          result := false;
        end

      // m must be greater than or equal to d
      else if( d > m ) then
        begin
          if( nil <> err ) then err^ := tr( 'D cannot be larger than M.' );
          result := false;
        end

      // m must be greater than or equal to n
      else if( n > m ) then
        begin
          if( nil <> err ) then err^ := tr( 'n cannot be larger than M.' );
          result := false;
        end

        // Don' do anything stupid
        else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
         begin
           if( nil <> err ) then err^ := tr( 'This probability density function cannot be used to represent a proportion or a probability.' );
           result := false;
        end

      else
        result := true
      ;

    end
  ;


  function TPdfHypergeometric.compare( chart: TChartFunction ):boolean;
    begin
      result := false;

      if( inherited compare( chart ) ) then
        begin
          if
            ( self.m = ( chart as TPdfHypergeometric ).m )
          and
            ( self.d = ( chart as TPdfHypergeometric ).d )
          and
            ( self.n = ( chart as TPdfHypergeometric ).n )
          then
            result := true
          ;
        end
      ;
    end
  ;


  function TPdfHypergeometric.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <hypergeometric>' + endl;
      result := result + _indent + '    <n>' + intToStr( n ) + '</n>' + endl;
      result := result + _indent + '    <d>' + intToStr( d ) + '</d>' + endl;
      result := result + _indent + '    <m>' + intToStr( m ) + '</m>' + endl;
      result := result + _indent + '  </hypergeometric>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfHypergeometric.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'hypergeometric' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setN( myStrToInt( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'n' )), -1, true ) );
          setD( myStrToInt( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'd' )), -1, true ) );
          setM( myStrToInt( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'm' )), -1, true ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfHypergeometric.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Hypergeometric"';
      dict['m']           := intToStr( m );
      dict['d']           := intToStr( d );
      dict['n']           := intToStr( n );

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


  function TPdfHypergeometric.createCopy(): TChartFunction;
    begin
      result := TPdfHypergeometric.create( self );
    end
  ;


  function TPdfHypergeometric.probDensity( k: integer ): double;
    begin
      if( nil <> @gslRanHypergeometricPdf ) then
        result := gslRanHypergeometricPdf( k, d, (m - d), n )
      else
        result := NaN
      ;
    end
  ;


  function TPdfHypergeometric.rand(): double;
    begin
      if( nil <> @gslRanHypergeometric ) then
        result := gslRanHypergeometric( gslRngPtr(), d, (m - d), n )
      else
        result := NaN
      ;
    end
  ;


  class function TPdfHypergeometric.rand( const n, d, m: integer ): double;
    begin
      if( nil <> @gslRanHypergeometric ) then
        result := gslRanHypergeometric( gslRngPtr(), d, (m - d), n )
      else
        result := NaN
      ;
    end
  ;


  (*
    This function is based on C++ code written by Agner Fog and released under the GPL.
    See http://www.agner.org/random/

    n is the total number of samples
    sum(d) is the total population
    d is the number in each bin
    output is the number of samples drawn from each bin
  *)
  class procedure TPdfHypergeometric.multiHypergeometricRand( n: integer; const d: TQIntegerVector; output: TQIntegerVector );
    var
      sum, x, y: integer;
      i: integer;
      c: integer;
    begin
      c := d.count;

      output.resize( c );

      if( 0 = c ) then
        exit
      ;

      if( 0 > n ) then
        begin
          raise exception.create( 'Parameter negative in TPdfHypergeometric.multiHypergeometricRand()' );
          exit;
        end
      ;

      // Compute total number in the population, doing some error checking along the way.
      sum := 0;
      for i := 0 to c - 1 do
        begin
          y := d.at( i );
          if( 0 > y ) then
            begin
              raise exception.create( 'Parameter negative in TPdfHypergeometric.multiHypergeometricRand()' );
              exit;
            end
          ;
          inc( sum, y );
        end
      ;

      if( n > sum ) then
        begin
          raise exception.create( 'Parameter negative in TPdfHypergeometric.multiHypergeometricRand()' );
          exit;
        end
      ;

      // Generate output by calling hypergeometric c - 1 times
      for i := 0 to ( c - 1 ) - 1 do
        begin
          y := d.at( i );
          x := trunc( RoundDblTo( TPdfHypergeometric.rand( n, y, sum ), 0 ) );
          dec( n, x );
          dec( sum, y );
          output[i] := x;
        end
      ;
     // Get the last one
     output[ c - 1 ] := n;

      (*
        // Original C code:
        void StochasticLib1::MultiHypergeometric (int32_t * destination, int32_t * source, int32_t n, int colors) {
           int32_t sum, x, y;
           int i;
           if (n < 0 || colors < 0) FatalError("Parameter negative in multihypergeo function");
           if (colors == 0) return;

           // compute total number of balls
           for (i=0, sum=0; i<colors; i++) {
              y = source[i];
              if (y < 0) FatalError("Parameter negative in multihypergeo function");
              sum += y;
           }
           if (n > sum) FatalError("n > sum in multihypergeo function");

           for (i=0; i<colors-1; i++) {
              // generate output by calling hypergeometric colors-1 times
              y = source[i];
              x = Hypergeometric(n, y, sum);
              n -= x;
              sum -= y;
              destination[i] = x;
           }
           // get the last one
           destination[i] = n;
        }
      *)
    end
  ;


  function TPdfHypergeometric.getHistValues(): TMassHistogramValues;
    var
      i: integer;
      nonZeroAlreadyFound: boolean;
      p: double;
    begin
      if( nil = _histValues ) then
        begin
          _histValues := TMassHistogramValues.create();

          // Start at 0, and try every integer until
          // the last non-zero probability is found.
          i := 0;
          nonZeroAlreadyFound := false;

          while( true ) do
            begin
              p := probDensity( i );

              if( 0.000001 < p ) then
                begin
                  nonZeroAlreadyFound := true;
                  _histValues.append( TMassHistogramValue.create( i, p ) );
                  inc( i );
                end
              else if( nonZeroAlreadyFound ) then
                break
              else
                inc( i )
              ;
            end
          ;

        end
      ;

      buildCumulHistValues();

      result := _histValues;
    end
  ;


  procedure TPdfHypergeometric.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str
        + 'Hypergeometric PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 'n=' + intToStr( n )
        + 'd=' + intToStr( d )
        + 'm=' + intToStr( m )
      ;
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;

  
  procedure TPdfHypergeometric.setM( val: integer );
    begin
      _m := val;
      freeHistogramValues();
    end
  ;
  
  
  procedure TPdfHypergeometric.setD( val: integer );
    begin
      _d := val;
      freeHistogramValues();
    end
  ;
  
  
   procedure TPdfHypergeometric.setN( val: integer );
     begin
       _n := val;
       freeHistogramValues();
     end
  ;


  function TPdfHypergeometric.getDescr(): string;
    begin
      result := tr( 'Hypergeometric' ) + format( ' ( %d, %d, %d )', [n, d, m] );
    end
  ;


  function TPdfHypergeometric.getMean(): double;
    begin
      result := (n * d ) / m;
    end
  ;

  
  function TPdfHypergeometric.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfHypergeometric.getMin(): double;
    begin
      result := Math.max( 0, n + d - m );
    end
  ;

  function TPdfHypergeometric.getMax(): double;
    begin
      result := Math.min( n, d );
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Negative Binomial PDF
//-----------------------------------------------------------------------------
  constructor TPdfNegativeBinomial.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfNegativeBinomial.create( s: integer; p: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setS( s );
      setP( p );
      xUnits := u;
    end
  ;


  constructor TPdfNegativeBinomial.create( const src: TPdfNegativeBinomial );
    begin
      inherited create( src );
      // initialize() is called by the inherited constructor

      setS( src.s );
      setP( src.p );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfNegativeBinomial.initialize();
    begin
      setPdfType( PdfNegativeBinomial );

      setS( -1 );
      setP( NaN );

      _infLeftTail := false;
      _infRightTail := true;

      freeHistogramValues();
    end
  ;


  destructor TPdfNegativeBinomial.destroy();
    begin
      inherited destroy();
    end
  ;


  procedure TPdfNegativeBinomial.setParams( s: integer; p: double; u: TChartUnitType = UUnknown );
    begin
      setS( s );
      setP( p );
      xUnits := u;
    end
  ;


  function TPdfNegativeBinomial.validate( err: pstring = nil ): boolean;
    begin
      if( s < 1 ) then
        begin
          if( nil <> err ) then err^ := tr( 'Number of successes must be greater than 0.' );
          result := false;
        end
      else if( isNan( p ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'Probability of success must be specified.' );
          result := false;
        end
      else if( not( isProbability( p ) ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'Probability of success must be between 0 and 1, inclusive.' );
          result := false;
        end
      else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
        begin
          if( nil <> err ) then err^ := tr( 'This probability density function cannot be used to represent a proportion or a probability.' );
          result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfNegativeBinomial.compare( chart: TChartFunction ):boolean;
    begin
      result := false;

      if( inherited compare( chart ) ) then
        begin
          if
            ( self._s = ( chart as TPdfNegativeBinomial ).s )
          and
            ( self._p = ( chart as TPdfNegativeBinomial ).p )
          then
            result := true
          ;
        end
      ;
    end
  ;


  function TPdfNegativeBinomial.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <negative-binomial>' + endl;
      result := result + _indent + '    <s>' + intToStr( s ) + '</s>' + endl;
      result := result + _indent + '    <p>' + usFloatToStr( p ) + '</p>' + endl;
      result := result + _indent + '  </negative-binomial>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfNegativeBinomial.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );


      fnElement := sdew.GetElementByName( element, 'negative-binomial' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          _s := myStrToInt( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 's' )), -1 );
          _p := usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'p' )), NaN );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfNegativeBinomial.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Negative-binomial"';
      dict['s']           := intToStr( s );
      dict['p']           := usFloatToStr( p );

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


  function TPdfNegativeBinomial.createCopy(): TChartFunction;
    begin
      result := TPdfNegativeBinomial.create( self );
    end
  ;


  function TPdfNegativeBinomial.probDensity( k: integer ): double;
    var
      ds: double;
    begin
      ds := self.s;

      if( nil <> @gslRanNegativeBinomialPdf ) then
        result := gslRanNegativeBinomialPdf( k, self.p, ds )
      else
        result := NaN
      ;
    end
  ;


  function TPdfNegativeBinomial.rand(): double;
    var
      ds: double;
    begin
      ds := self.s;

      if( nil <> @gslRanNegativeBinomial ) then
        result := gslRanNegativeBinomial( gslRngPtr(), self.p, ds )
      else
        result := NaN
      ;
    end
  ;


  function TPdfNegativeBinomial.getHistValues(): TMassHistogramValues;
    var
      i: integer;
      p: double;
      nonZeroAlreadyFound: boolean;
    begin
      if( nil = _histValues ) then
        begin
          _histValues := TMassHistogramValues.create();

          // Start at 0, and try every integer until
          // the last non-zero probability is found.
          i := 0;
          nonZeroAlreadyFound := false;

          while( true ) do
            begin
              p := probDensity( i );

              if( 0.000001 < p ) then
                begin
                  nonZeroAlreadyFound := true;
                  _histValues.append( TMassHistogramValue.create( i, p ) );
                  inc( i );
                end
              else if( nonZeroAlreadyFound ) then
                break
              else
                inc( i )
              ;
            end
          ;

        end
      ;

      buildCumulHistValues();

      result := _histValues;
    end
  ;


  procedure TPdfNegativeBinomial.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str
        + 'Negative binomial PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 's=' + intToStr( s )
        + ', p=' + usFloatToStr( p )
      ;
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;


  procedure TPdfNegativeBinomial.setS( val: integer ); begin _s := val; freeHistogramValues(); end;
  procedure TPdfNegativeBinomial.setP( val: double ); begin _p := val; freeHistogramValues(); end;


  function TPdfNegativeBinomial.getDescr(): string;
    begin
      result := tr( 'Negative binomial' ) + format( ' ( %d, %f )', [s, p] );
    end
  ;


  function TPdfNegativeBinomial.getMean(): double;
    begin
      result := s * (1 - p) / p;
    end
  ;


  function TPdfNegativeBinomial.getHasMean(): boolean;
    begin
      result := true;
    end
  ;

  function TPdfNegativeBinomial.getMin(): double;
    begin
      result := 0;
    end
  ;

  function TPdfNegativeBinomial.getMax(): double;
    begin
      result := infinity;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// Poisson PDF
//-----------------------------------------------------------------------------
  constructor TPdfPoisson.create();
    begin
      inherited create();
      // initialize() is called by the inherited constructor
    end
  ;


  constructor TPdfPoisson.create( mean: double; u: TChartUnitType = UUnknown );
    begin
      inherited create();
      // initialize() is called by the inherited constructor

      setMean( mean );
      xUnits := u;
    end
  ;


  constructor TPdfPoisson.create( const src: TPdfPoisson );
    begin
      inherited create( src );
      // initialize() is called by the inherited constructor

      setMean( src.mean );
      xUnits := src.xUnits;
    end
  ;


  procedure TPdfPoisson.initialize();
    begin
      setPdfType( PdfPoisson );

      setMean( NaN );

      _infLeftTail := false;
      _infRightTail := true;
    end
  ;


  destructor TPdfPoisson.destroy();
    begin
      inherited destroy();
    end
  ;


  function TPdfPoisson.validate( err: pstring = nil ): boolean;
    begin
      if( isNan( mean ) ) then
        begin
          if( err <> nil ) then err^ := tr( 'Mean must be specified.' );
          result := false;
        end
      else if( 0 >= mean ) then
        begin
          if( nil <> err ) then err^ := tr( 'The mean must be greater than 0.' );
          result := false;
        end
       else if( ( UProportion = xUnits ) or ( UProbability = XUnits ) ) then
         begin
           if( nil <> err ) then err^ := tr( 'This probability density function cannot be used to represent a proportion or a probability.' );
           result := false;
        end
      else
        result := true
      ;
    end
  ;


  function TPdfPoisson.compare( chart: TChartFunction ):boolean;
    begin
      result := false;

      if( inherited compare( chart ) ) then
        begin
          if
            ( self.mean = ( chart as TPdfPoisson ).mean )
          then
            result := true
          ;
        end
      ;
    end
  ;


  function TPdfPoisson.ssXml( const indent: integer ): string;
    begin
      setIndent( indent );
      result := '';

      if( 0 < length( self.name ) ) then
        result := result + _indent + '<probability-density-function name="' + encodeXml( self.name ) + '">' + endl
      else
        result := result + _indent + '<probability-density-function>' + endl
      ;
      result := result + _indent + '  <poisson>' + endl;
      result := result + _indent + '    <mean>' + usFloatToStr( mean ) + '</mean>' + endl;
      result := result + _indent + '  </poisson>' + endl;
      result := result + _indent + '  ' + chartUnitTypeAsSSXml( xUnits ) + endl;
      result := result + _indent + '</probability-density-function>' + endl;
    end
  ;


  function TPdfPoisson.importXml( element: Pointer; sdew: TSdew ): boolean;
    var
      fnElement, nameAttribute: Pointer;
    begin
      inherited importXml( element, sdew );

      fnElement := sdew.GetElementByName( element, 'poisson' );

      if( nil <> fnElement ) then
        begin
          // Old-style XML may have the name attached to the type tag
          if( strIsEmpty( self.name ) ) then
            begin
              nameAttribute := sdew.getAttributeByName( fnElement, 'name' );
              if( nil <> nameAttribute ) then
                name := sdew.GetElementAttribute( fnElement, 'name' )
              ;
            end
          ;

          setMean( usStrToFloat( Sdew.GetElementContents(Sdew.GetElementByName( fnElement, 'mean' )), NaN ) );
        end
      ;

      result := isValid;
    end
  ;


  {$IFDEF DATABASE_ENABLED}
  function TPdfPoisson.populateDatabase( db: TSqlDatabase; update: boolean = false ): integer;
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

      dict['chartType']   := '"Poisson"';
      dict['mean']        := usFloatToStr( mean );

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


  function TPdfPoisson.createCopy(): TChartFunction;
    begin
      result := TPdfPoisson.create( self );
    end
  ;


  function TPdfPoisson.probDensity( k: integer ): double;
    begin
      if( nil <> @gslRanPoissonPdf ) then
        result := gslRanPoissonPdf( k, self.mean )
      else
        result := NaN
      ;
    end
  ;


  function TPdfPoisson.rand(): double;
    begin
      if( nil <> @gslRanBinomial ) then
        result := gslRanPoisson( gslRngPtr(), self.mean )
      else
        result := NaN
      ;
    end
  ;


  function TPdfPoisson.getHistValues(): TMassHistogramValues;
    var
      i: integer;
      nonZeroAlreadyFound: boolean;
      p: double;
    begin
      if( nil = _histValues ) then
        begin
          _histValues := TMassHistogramValues.create();

          // Start at 0, and try every integer until
          // the last non-zero probability is found.
          i := 0;
          nonZeroAlreadyFound := false;

          while( true ) do
            begin
              p := probDensity( i );

              if( 0.000001 < p ) then
                begin
                  nonZeroAlreadyFound := true;
                  _histValues.append( TMassHistogramValue.create( i, p ) );
                  inc( i );
                end
              else if( nonZeroAlreadyFound ) then
                break
              else
                inc( i )
              ;
            end
          ;

        end
      ;

      buildCumulHistValues();

      result := _histValues;
    end
  ;


  procedure TPdfPoisson.debug();
    var
      str: string;
    begin
      if( strIsEmpty( name ) ) then
        str := '(unnamed)' + endl
      else
        str := name + endl
      ;
      str := str
        + 'Poisson PDF (id ' + intToStr( self.id ) + ') :' + endl
        + 'mean=' + usFloatToStr( mean )
      ;
      dbcout( str + endl + 'units of ' + xUnitsString(), true );
    end
  ;

  
  procedure TPdfPoisson.setMean( val: double );
    begin
      _mean := val;
      freeHistogramValues();
    end
  ;


  function TPdfPoisson.getDescr(): string;
    begin
      result := tr( 'Poisson' ) + format( ' ( %f )', [mean] );
    end
  ;


  function TPdfPoisson.getMean(): double;
    begin
      result := _mean;
    end
  ;

  
  function TPdfPoisson.getHasMean(): boolean;
    begin
      result := true;
    end
  ;


  function TPdfPoisson.getMin(): double;
    begin
      result := 0.0;
    end
  ;

  function TPdfPoisson.getMax(): double;
    begin
      result := infinity;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// XML import
//-----------------------------------------------------------------------------
  {* Record type used for callback functions (see below). }
  type RXmlPdfIndex = record
    name:  string;   /// name of the PDF type as spelled in XML import files
    index: TPdfType; /// corresponding TPdfType enumerated value
  end;

  {* Number of probability density functions that can be imported from XML. }
  const NumXmlImportedPdfTypes = 23;

  {* Probability density function names and their associated TPdfType equivalent. 
     These are used by createPdfFromXml as it parses the PDF portion of an XML file.
     Names must appear here exactly as they are indicated in the XML schema.
  }
  const XmlImportedPdfTypes:  array[ 0..NumXmlImportedPdfTypes - 1 ] of RXmlPdfIndex = (
                //  ADD  PdfSpGaniTalbert in NAADSM 5.0
    // Continuous types
    ( name: 'beta';               index: PdfBeta            ),
    ( name: 'beta-pert';          index: PdfBetaPert        ),
    ( name: 'exponential';        index: PdfExponential     ),
    ( name: 'gamma';              index: PdfGamma           ),
    ( name: 'gaussian';           index: PdfGaussian        ),
    ( name: 'histogram';          index: PdfHistogram       ),
    ( name: 'inverse-gaussian';   index: PdfInverseGaussian ),
    ( name: 'logistic';           index: PdfLogistic        ),
    ( name: 'loglogistic';        index: PdfLoglogistic     ),
    ( name: 'lognormal';          index: PdfLognormal       ),
    ( name: 'pareto';             index: PdfPareto          ),
    ( name: 'pearson5';           index: PdfPearson5        ),
    ( name: 'piecewise';          index: PdfPiecewise       ),
    ( name: 'point';              index: PdfPoint           ),
    ( name: 'triangular';         index: PdfTriangular      ),
    ( name: 'uniform';            index: PdfUniform         ),
    ( name: 'weibull';            index: PdfWeibull         ),

    // Discrete types
    ( name: 'bernoulli';          index: PdfBernoulli        ),
    ( name: 'binomial';           index: PdfBinomial         ),
    ( name: 'discrete-uniform';   index: PdfDiscreteUniform  ),
    ( name: 'hypergeometric';     index: PdfHypergeometric   ),
    ( name: 'negative-binomial';  index: PdfNegativeBinomial ),
    ( name: 'poisson';            index: PdfPoisson          )
  );


  // Exposed functions for PDF XML import.

  function createPdfFromXml( xml: string ): TPdf;
    var
      sdew: TSdew;
      element: Pointer;
    begin
      xml := stripXmlHeaders( xml );
      sdew := TSdew.createFromString( xml );
      element := sdew.GetRootTree();

      try
        // createPdfFromXml will raise an exception if no pdf is found.
        result := createPdfFromXml( element, sdew, nil ) as TPdf;
      except
        result := nil;
      end;
      
      sdew.Free();
    end
  ;


  function createPdfFromXml( element: pointer; sdew: TSdew ): TPdf;
    begin
      result := createPdfFromXml( element, sdew, nil ) as TPdf;
    end
  ;

  
  function createPdfFromXml( Element: pointer; Sdew:TSdew; unused: pointer ): TObject;
    var
      i,j: integer;

      ret_val: TPdf;
      pdfName: string;
      nameFound: boolean;
      index: TPdfType;
      Count: integer;
      pdfElement: pointer;
      tmp: pointer;
    begin
      ret_val := nil; // until set below

      if( 'probability-density-function' = sdew.GetElementName( element ) ) then
        begin
          // Do nothing: this is new-style XML, and we're ready to parse it.
        end
      else
        begin
          tmp := sdew.GetElementByName( element, 'probability-density-function' );

          if( nil <> tmp ) then
            begin
              // This is still new-style XML, but this function is being called from the
              // soon-to-be-obsolete TXMLReader, which sends the parent element rather than
              // the pdf element to this function to be parsed.
              element := tmp;
            end
          else
            begin
              // Assume that this is old-style XML.
              // Parsing will be pretty painless, except for histogram and piecewise distributions.
              // Functions that will parse these two distributions know what to expect, and nothing
              // more needs to happen here.
            end
          ;
        end
      ;

      
      index := PdfUndefined;
      Count := Sdew.GetElementCount( element );
      nameFound := false;
      
      for  j := 0 to Count - 1 do
        begin
          pdfElement := sdew.GetElementByIndex( element, j );
          pdfName := Trim( sdew.GetElementName( pdfElement ) );

          for i := 0 to length( XmlImportedPdfTypes ) - 1 do
            begin
              if ( 0 = AnsiCompareStr( XmlImportedPdfTypes[i].name, pdfName ) ) then
                begin
                  nameFound := true;
                  index := XmlImportedPdfTypes[i].index;
                  break;
                end
              ;
            end
          ;
          if ( nameFound ) then
            break;
        end
      ;

      if ( nameFound ) then
        begin
          case index of
            // Continuous types
            PDFBeta:              ret_val := TPdfBeta.create();
            PdfBetaPERT:          ret_val := TPdfBetaPERT.create();
            PdfExponential:       ret_val := TPdfExponential.create();
            PdfGamma:             ret_val := TPdfGamma.create();
            PdfGaussian:          ret_val := TPdfGaussian.create();
            PdfHistogram:         ret_val := TPdfHistogram.create();
            PdfInverseGaussian:   ret_val := TPdfInverseGaussian.create();
            PdfLogistic:          ret_val := TPdfLogistic.create();
            PdfLogLogistic:       ret_val := TPdfLogLogistic.create();
            PdfLognormal:         ret_val := TPdfLognormal.create();
            PdfPareto:            ret_val := TPdfPareto.create();
            PdfPearson5:          ret_val := TPdfPearson5.create();
            PdfPiecewise:         ret_val := TPdfPiecewise.create();
            PdfPoint:             ret_val := TPdfPoint.create();
            PdfTriangular:        ret_val := TPdfTriangular.create();
            PdfUniform:           ret_val := TPdfUniform.create();
            PdfWeibull:           ret_val := TPdfWeibull.create();
            
            // Discrete types
            PdfBernoulli:         ret_val := TPdfBernoulli.create();
            PdfBinomial:          ret_val := TPdfBinomial.create();
            PdfDiscreteUniform:   ret_val := TPdfDiscreteUniform.create();            
            PdfHypergeometric:    ret_val := TPdfHypergeometric.create();
            PdfNegativeBinomial:  ret_val := TPdfNegativeBinomial.create();
            PdfPoisson:           ret_val := TPdfPoisson.create();            
          else
            begin
              raise exception.create( 'There is a problem importing this PDF function, (' + pdfName +'), from the XML.  No similar function was found in this program.' );
            end
          ;
          end
        ;
        if( assigned( ret_val ) ) then
          begin
            if( not ret_val.importXml( element, Sdew ) )  then
              freeAndNil( ret_val )
            ;
          end
        end
      else
        begin
          raise exception.create( 'No PDF functions found in this section, ( ' + Sdew.GetElementName( Element ) +' ), of XML.' );
          ret_val := nil;
        end
      ;
      
      result := ret_val;
    end
  ;
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// DLL handling
//-----------------------------------------------------------------------------
  procedure freeLoadedPointers();
    begin
      // Pointers loaded from the GSL
      //-----------------------------
      // Continuous types
      //-----------------
      gslRanExponentialPDF := nil;
      gslCdfExponentialP := nil;
      gslCdfExponentialPinv := nil;
      gslRanExponential := nil;

      gslRanGammaPdf := nil;
      gslCdfGammaP := nil;
      gslCdfGammaPinv := nil;
      gslRanGamma := nil;

      gslRanGaussianPdf := nil;
      gslCdfGaussianP := nil;
      gslCdfGaussianPinv := nil;
      gslRanGaussian := nil;

      gslRanLogisticPDF := nil;
      gslCdfLogisticP := nil;
      gslCdfLogisticPinv := nil;
      gslRanLogistic := nil;

      gslRanLognormalPDF := nil;
      gslCdfLognormalP := nil;
      gslCdfLognormalPinv := nil;
      gslRanLognormal := nil;

      gslRanParetoPDF := nil;
      gslCdfParetoP := nil;
      gslCdfParetoPinv := nil;
      gslRanPareto := nil;

      gslCdfUniformP := nil;
      gslCdfUniformPinv := nil;
      gslRanUniform := nil;

      gslRanWeibullPDF := nil;
      gslCdfWeibullP := nil;
      gslCdfWeibullPinv := nil;
      gslRanWeibull := nil;

      // Discrete types
      //---------------
      gslRanBinomialPdf := nil;
      gslRanBinomial := nil;

      gslRanHypergeometricPdf := nil;
      gslRanHypergeometric := nil;

      gslRanNegativeBinomialPdf := nil;
      gslRanNegativeBinomial := nil;

      gslRanPoissonPdf := nil;
      gslRanPoisson := nil;
      

      // Pointers loaded from libAPHI
      //-----------------------------
      // Continuous types
      //-----------------
      aphi_beta_pdf := nil;
      aphi_beta_cdf := nil;
      aphi_beta_inverse_cdf := nil;
      aphi_beta_rand := nil;

      aphi_beta_pert_pdf := nil;
      aphi_beta_pert_cdf := nil;
      aphi_beta_pert_inverse_cdf := nil;
      aphi_beta_pert_rand := nil;

      aphi_create_histogram := nil;
      aphi_free_histogram := nil;
      aphi_histogram_mean := nil;
      aphi_histogram_pdf := nil;
      aphi_histogram_cdf := nil;
      aphi_histogram_inverse_cdf := nil;
      aphi_histogram_rand := nil;

      aphi_inverse_gaussian_pdf := nil;
      aphi_inverse_gaussian_cdf := nil;
      aphi_inverse_gaussian_inverse_cdf := nil;

      aphi_loglogistic_pdf := nil;
      aphi_loglogistic_cdf := nil;
      aphi_loglogistic_inverse_cdf := nil;
      aphi_loglogistic_rand := nil;

      aphi_pearson5_pdf := nil;
      aphi_pearson5_cdf := nil;
      aphi_pearson5_inverse_cdf := nil;
      aphi_pearson5_rand := nil;

      aphi_triangular_cdf := nil;
      aphi_triangular_inverse_cdf := nil;
      aphi_triangular_rand := nil;
    end
  ;


  function loadGslPointers(): boolean;
    var
      gslHandle: THandle; //Handle used to open the DLL.  Defined in unit Windows.
    begin
      try
        gslHandle := loadLibrary( 'libgsl.dll' );
      except
        result := false;
        freeLoadedPointers();
        exit;
      end;

      if( gslHandle >= 32 ) then // libraries were successfully loaded.  Assign function pointers now.
        begin
          result := true;

          //-------------------------------------------------------------------
          // CONTINUOUS TYPES
          //-------------------------------------------------------------------
          // Exponential PDF
          //----------------
          gslRanExponentialPdf := GetProcAddress( gslHandle, 'gsl_ran_exponential_pdf' );
          if( nil = @gslRanExponentialPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_exponential_pdf' + endl;
              result := false;
            end
          ;

          gslCdfExponentialP := GetProcAddress( gslHandle, 'gsl_cdf_exponential_P' );
          if( nil = @gslCdfExponentialP ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_exponential_P' + endl;
              result := false;
            end
          ;

          gslCdfExponentialPinv := GetProcAddress( gslHandle, 'gsl_cdf_exponential_Pinv' );
          if( nil = @gslCdfExponentialPinv ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_exponential_Pinv' + endl;
              result := false;
            end
          ;

          gslRanExponential := GetProcAddress( gslHandle, 'gsl_ran_exponential' );
          if( nil = @gslRanExponential ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_exponential' + endl;
              result := false;
            end
          ;

          // Gamma PDF
          //----------
          gslRanGammaPdf := GetProcAddress( gslHandle, 'gsl_ran_gamma_pdf' );
          if( nil = @gslRanGammaPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_gamma_pdf' + endl;
              result := false;
            end
          ;

          gslCdfGammaP := GetProcAddress( gslHandle, 'gsl_cdf_gamma_P' );
          if( nil = @gslCdfGammaP ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_gamma_P' + endl;
              result := false;
            end
          ;

          gslCdfGammaPinv := GetProcAddress( gslHandle, 'gsl_cdf_gamma_Pinv' );
          if( nil = @gslCdfGammaPinv ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_gamma_Pinv' + endl;
              result := false;
            end
          ;

          gslRanGamma := GetProcAddress( gslHandle, 'gsl_ran_gamma' );
          if( nil = @gslRanGamma ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_gamma' + endl;
              result := false;
            end
          ;

          // Gaussian PDF
          //-------------
          gslRanGaussianPdf := GetProcAddress( gslHandle, 'gsl_ran_gaussian_pdf' );
          if( nil = @gslRanGaussianPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_gaussian_pdf' + endl;
              result := false;
            end
          ;

          gslCdfGaussianP := GetProcAddress( gslHandle, 'gsl_cdf_gaussian_P' );
          if( nil = @gslCdfGaussianP ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_gaussian_P' + endl;
              result := false;
            end
          ;

          gslCdfGaussianPinv := GetProcAddress( gslHandle, 'gsl_cdf_gaussian_Pinv' );
          if( nil = @gslCdfGaussianPinv ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_gaussian_Pinv' + endl;
              result := false;
            end
          ;

          gslRanGaussian := GetProcAddress( gslHandle, 'gsl_ran_gaussian' );
          if( nil = @gslRanGaussian ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_gaussian' + endl;
              result := false;
            end
          ;

          // Logistic PDF
          //-------------
          gslRanLogisticPdf := GetProcAddress( gslHandle, 'gsl_ran_logistic_pdf' );
          if( nil = @gslRanLogisticPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_logistic_pdf' + endl;
              result := false;
            end
          ;

          gslCdfLogisticP := GetProcAddress( gslHandle, 'gsl_cdf_logistic_P' );
          if( nil = @gslCdfLogisticP ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_logistic_P' + endl;
              result := false;
            end
          ;

          gslCdfLogisticPinv := GetProcAddress( gslHandle, 'gsl_cdf_logistic_Pinv' );
          if( nil = @gslCdfLogisticPinv ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_logistic_Pinv' + endl;
              result := false;
            end
          ;

          gslRanLogistic := GetProcAddress( gslHandle, 'gsl_ran_logistic' );
          if( nil = @gslRanLogistic ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_logistic' + endl;
              result := false;
            end
          ;

          // Lognormal PDF
          //--------------
          gslRanLognormalPdf := GetProcAddress( gslHandle, 'gsl_ran_lognormal_pdf' );
          if( nil = @gslRanLognormalPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_lognormal_pdf' + endl;
              result := false;
            end
          ;

          gslCdfLognormalP := GetProcAddress( gslHandle, 'gsl_cdf_lognormal_P' );
          if( nil = @gslCdfLognormalP ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_lognormal_P' + endl;
              result := false;
            end
          ;

          gslCdfLognormalPinv := GetProcAddress( gslHandle, 'gsl_cdf_lognormal_Pinv' );
          if( nil = @gslCdfLognormalPinv ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_lognormal_Pinv' + endl;
              result := false;
            end
          ;

          gslRanLognormal := GetProcAddress( gslHandle, 'gsl_ran_lognormal' );
          if( nil = @gslRanLognormal ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_lognormal' + endl;
              result := false;
            end
          ;

          // Pareto PDF
          //-----------
          gslRanParetoPdf := GetProcAddress( gslHandle, 'gsl_ran_pareto_pdf' );
          if( nil = @gslRanParetoPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_pareto_pdf' + endl;
              result := false;
            end
          ;

          gslCdfParetoP := GetProcAddress( gslHandle, 'gsl_cdf_pareto_P' );
          if( nil = @gslCdfParetoP ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_pareto_P' + endl;
              result := false;
            end
          ;

          gslCdfParetoPinv := GetProcAddress( gslHandle, 'gsl_cdf_pareto_Pinv' );
          if( nil = @gslCdfParetoPinv ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_pareto_Pinv' + endl;
              result := false;
            end
          ;

          gslRanPareto := GetProcAddress( gslHandle, 'gsl_ran_pareto' );
          if( nil = @gslRanPareto ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_pareto' + endl;
              result := false;
            end
          ;

          // Piecewise PDF
          //--------------
          gsl_poly_solve_quadratic := GetProcAddress( gslHandle, 'gsl_poly_solve_quadratic' );
          if( nil = @gsl_poly_solve_quadratic ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_poly_solve_quadratic' + endl;
              result := false;
            end
          ;

          // Uniform PDF
          //------------
          gslCdfUniformP := getProcAddress( gslHandle, 'gsl_cdf_flat_P' );
          if( nil = @gslCdfUniformP ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_flat_P' + endl;
              result := false;
            end
          ;

          gslCdfUniformPinv := getProcAddress( gslHandle, 'gsl_cdf_flat_Pinv' );
          if( nil = @gslCdfUniformPinv ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_flat_Pinv' + endl;
              result := false;
            end
          ;

          gslRanUniform := GetProcAddress( gslHandle, 'gsl_ran_flat' );
          if( nil = @gslRanUniform ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_flat' + endl;
              result := false;
            end
          ;

          // Weibull PDF
          //-------------
          gslRanWeibullPdf := GetProcAddress( gslHandle, 'gsl_ran_weibull_pdf' );
          if( nil = @gslRanWeibullPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_weibull_pdf' + endl;
              result := false;
            end
          ;

          gslCdfWeibullP := GetProcAddress( gslHandle, 'gsl_cdf_weibull_P' );
          if( nil = @gslCdfWeibullP ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_weibull_P' + endl;
              result := false;
            end
          ;

          gslCdfWeibullPinv := GetProcAddress( gslHandle, 'gsl_cdf_weibull_Pinv' );
          if( nil = @gslCdfWeibullPinv ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_cdf_weibull_Pinv' + endl;
              result := false;
            end
          ;

          gslRanWeibull := GetProcAddress( gslHandle, 'gsl_ran_weibull' );
          if( nil = @gslRanWeibull ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_weibull' + endl;
              result := false;
            end
          ;
          //-------------------------------------------------------------------


          //-------------------------------------------------------------------
          // DISCRETE TYPES
          //-------------------------------------------------------------------
          // Binomial PDF
          //-------------
          gslRanBinomialPdf := GetProcAddress( gslHandle, 'gsl_ran_binomial_pdf' );
          if( nil = @gslRanBinomialPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_binomial_pdf' + endl;
              result := false;
            end
          ;

          gslRanBinomial := GetProcAddress( gslHandle, 'gsl_ran_binomial' );
          if( nil = @gslRanBinomial ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_binomial' + endl;
              result := false;
            end
          ;


          // Hypergeometric PDF
          //-------------------
          gslRanHypergeometricPdf := GetProcAddress( gslHandle, 'gsl_ran_hypergeometric_pdf' );
          if( nil = @gslRanHypergeometricPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_hypergeometric_pdf' + endl;
              result := false;
            end
          ;

          gslRanHypergeometric := GetProcAddress( gslHandle, 'gsl_ran_hypergeometric' );
          if( nil = @gslRanHypergeometric ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_hypergeometric' + endl;
              result := false;
            end
          ;


          // Negative Binomial PDF
          //----------------------
          gslRanNegativeBinomialPdf := GetProcAddress( gslHandle, 'gsl_ran_negative_binomial_pdf' );
          if( nil = @gslRanNegativeBinomialPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_negative_binomial_pdf' + endl;
              result := false;
            end
          ;

          gslRanNegativeBinomial := GetProcAddress( gslHandle, 'gsl_ran_negative_binomial' );
          if( nil = @gslRanNegativeBinomial ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_negative_binomial' + endl;
              result := false;
            end
          ;


          // Poisson PDF
          //------------
          gslRanPoissonPdf := getProcAddress( gslHandle, 'gsl_ran_poisson_pdf' );
          if( nil = @gslRanPoissonPdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_poisson_pdf' + endl;
              result := false;
            end
          ;

          gslRanPoisson := getProcAddress( gslHandle, 'gsl_ran_poisson' );
          if( nil = @gslRanPoisson ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION gsl_ran_poisson' + endl;
              result := false;
            end
          ;
          //-------------------------------------------------------------------
        end
      else
        begin
          dllLoadErrMsgs := dllLoadErrMsgs + 'GSL DLL could not be opened.' + endl;
          result := false;
        end
      ;

      dbcout( 'Gsl pdfs loaded', DBSHOWMSG );
      dbcout( result, DBSHOWMSG );
    end
  ;



  function loadAPHIPointers(): boolean;
    var
      pdfHandle: THandle; //Handle used to open the DLL.  Defined in unit Windows.

      lib_version: function(): pchar; cdecl;
      libVersion: string;
    begin
      try
        pdfHandle := loadLibrary( 'libaphi.dll' ); 
      except
        result := false;
        freeLoadedPointers();
        exit;
      end;

      if(  pdfHandle >= 32 ) then // libAPHI was successfully loaded.  Assign function pointers now.
        begin
          result := true;

          // Check the version of libAPHI
          //-----------------------------
          lib_version := getProcAddress( pdfHandle, 'lib_version' );
          if( nil = @lib_version ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION lib_version' + endl;
              result := false;
              exit;
            end
          else
            begin
              libVersion := lib_version();
              if( REQUIRED_APHI_DLL_VERSION <> libVersion ) then
                begin
                  dllLoadErrMsgs := dllLoadErrMsgs + 'INCOMPATIBLE LibAPHI DLL VERSION: ' + libVersion + endl;
                  result := false;
                end
              ;
            end
          ;

          //-------------------------------------------------------------------
          // CONTINUOUS TYPES
          //-------------------------------------------------------------------
          // Beta PDF
          //---------
          aphi_beta_pdf := getProcAddress( pdfHandle, 'aphi_beta_pdf' );
          if( nil = @aphi_beta_pdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_beta_pdf' + endl;
              result := false;
            end
          ;

          aphi_beta_cdf := getProcAddress( pdfHandle, 'aphi_beta_cdf' );
          if( nil = @aphi_beta_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_beta_cdf' + endl;
              result := false;
            end
          ;

          aphi_beta_inverse_cdf := getProcAddress( pdfHandle, 'aphi_beta_inverse_cdf' );
          if( nil = @aphi_beta_inverse_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_beta_inverse_cdf' + endl;
              result := false;
            end
          ;

          aphi_beta_rand := getProcAddress( pdfHandle, 'aphi_beta_rand' );
          if( nil = @aphi_beta_rand ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_beta_rand' + endl;
              result := false;
            end
          ;

          // BetaPERT PDF
          //-------------
          aphi_beta_pert_pdf := getProcAddress( pdfHandle, 'aphi_beta_pert_pdf' );
          if( nil = @aphi_beta_pert_pdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_beta_pert_pdf' + endl;
              result := false;
            end
          ;

          aphi_beta_pert_cdf := getProcAddress( pdfHandle, 'aphi_beta_pert_cdf' );
          if( nil = @aphi_beta_pert_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_beta_pert_cdf' + endl;
              result := false;
            end
          ;

          aphi_beta_pert_inverse_cdf := getProcAddress( pdfHandle, 'aphi_beta_pert_inverse_cdf' );
          if( nil = @aphi_beta_pert_inverse_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_beta_pert_inverse_cdf' + endl;
              result := false;
            end
          ;

          aphi_beta_pert_rand := getProcAddress( pdfHandle, 'aphi_beta_pert_rand' );
          if( nil = @aphi_beta_pert_rand ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_beta_pert_rand' + endl;
              result := false;
            end
          ;

          // Histogram PDF
          //--------------
          aphi_create_histogram := getProcAddress( pdfHandle, 'aphi_create_histogram' );
          if( nil = @aphi_create_histogram ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_create_histogram' + endl;
              result := false;
            end
          ;

          aphi_free_histogram := getProcAddress( pdfHandle, 'aphi_free_histogram' );
          if( nil = @aphi_free_histogram ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_free_histogram' + endl;
              result := false;
            end
          ;

          aphi_histogram_mean := getProcAddress( pdfHandle, 'aphi_histogram_mean' );
          if( nil = @aphi_histogram_mean ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_histogram_mean' + endl;
              result := false;
            end
          ;

          aphi_histogram_pdf := getProcAddress( pdfHandle, 'aphi_histogram_pdf' );
          if( nil = @aphi_histogram_pdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_histogram_pdf' + endl;
              result := false;
            end
          ;

          aphi_histogram_cdf := getProcAddress( pdfHandle, 'aphi_histogram_cdf' );
          if( nil = @aphi_histogram_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_histogram_cdf' + endl;
              result := false;
            end
          ;

          aphi_histogram_inverse_cdf := getProcAddress( pdfHandle, 'aphi_histogram_inverse_cdf' );
          if( nil = @aphi_histogram_inverse_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_histogram_inverse_cdf' + endl;
              result := false;
            end
          ;

          aphi_histogram_rand := getProcAddress( pdfHandle, 'aphi_histogram_rand' );
          if( nil = @aphi_histogram_rand ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_histogram_rand' + endl;
              result := false;
            end
          ;

          // Inverse Gaussian PDF
          //---------------------
          aphi_inverse_gaussian_pdf := getProcAddress( pdfHandle, 'aphi_inverse_gaussian_pdf' );
          if( nil = @aphi_inverse_gaussian_pdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_inverse_gaussian_pdf' + endl;
              result := false;
            end
          ;

          aphi_inverse_gaussian_cdf := getProcAddress( pdfHandle, 'aphi_inverse_gaussian_cdf' );
          if( nil = @aphi_inverse_gaussian_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_inverse_gaussian_cdf' + endl;
              result := false;
            end
          ;

          aphi_inverse_gaussian_inverse_cdf := getProcAddress( pdfHandle, 'aphi_inverse_gaussian_inverse_cdf' );
          if( nil = @aphi_inverse_gaussian_inverse_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_inverse_gaussian_inverse_cdf' + endl;
              result := false;
            end
          ;

          // Loglogistic PDF
          //----------------
          aphi_loglogistic_pdf := getProcAddress( pdfHandle, 'aphi_loglogistic_pdf' );
          if( nil = @aphi_loglogistic_pdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_loglogistic_pdf' + endl;
              result := false;
            end
          ;

          aphi_loglogistic_cdf := getProcAddress( pdfHandle, 'aphi_loglogistic_cdf' );
          if( nil = @aphi_loglogistic_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_loglogistic_cdf' + endl;
              result := false;
            end
          ;

          aphi_loglogistic_inverse_cdf := getProcAddress( pdfHandle, 'aphi_loglogistic_inverse_cdf' );
          if( nil = @aphi_loglogistic_inverse_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_loglogistic_inverse_cdf' + endl;
              result := false;
            end
          ;


          aphi_loglogistic_rand := getProcAddress( pdfHandle, 'aphi_loglogistic_rand' );
          if( nil = @aphi_loglogistic_rand ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_loglogistic_rand' + endl;
              result := false;
            end
          ;

          // Pearson5 PDF
          //-------------
          aphi_pearson5_pdf := getProcAddress( pdfHandle, 'aphi_pearson5_pdf' );
          if( nil = @aphi_pearson5_pdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_pearson5_pdf' + endl;
              result := false;
            end
          ;

          aphi_pearson5_cdf := getProcAddress( pdfHandle, 'aphi_pearson5_cdf' );
          if( nil = @aphi_pearson5_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_pearson5_cdf' + endl;
              result := false;
            end
          ;

          aphi_pearson5_inverse_cdf := getProcAddress( pdfHandle, 'aphi_pearson5_inverse_cdf' );
          if( nil = @aphi_pearson5_inverse_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_pearson5_inverse_cdf' + endl;
              result := false;
            end
          ;

          aphi_pearson5_rand := getProcAddress( pdfHandle, 'aphi_pearson5_rand' );
          if( nil = @aphi_pearson5_rand ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_pearson5_rand' + endl;
              result := false;
            end
          ;

          // Triangular PDF
          //---------------
          aphi_triangular_cdf := getProcAddress( pdfHandle, 'aphi_triangular_cdf' );
          if( nil = @aphi_triangular_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_triangular_cdf' + endl;
              result := false;
            end
          ;

          aphi_triangular_inverse_cdf := getProcAddress( pdfHandle, 'aphi_triangular_inverse_cdf' );
          if( nil = @aphi_triangular_inverse_cdf ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_triangular_inverse_cdf' + endl;
              result := false;
            end
          ;

          aphi_triangular_rand := getProcAddress( pdfHandle, 'aphi_triangular_rand' );
          if( nil = @aphi_triangular_rand ) then
            begin
              dllLoadErrMsgs := dllLoadErrMsgs + 'MISSING FUNCTION aphi_triangular_rand' + endl;
              result := false;
            end
          ;
          //-------------------------------------------------------------------
          
        end
      else
        begin
          dllLoadErrMsgs := dllLoadErrMsgs + 'LibAPHI DLL could not be loaded.';
          result := false;
        end
      ;

      dbcout( 'APHI pdfs loaded', DBSHOWMSG );
      dbcout( result, DBSHOWMSG );
    end
  ;
//-----------------------------------------------------------------------------


initialization

    // Load the needed functions from the GSL and libAPHI DLLs
    dllLoadErrMsgs := '';
    gslLoaded := loadGslPointers();
    libAPHILoaded := loadAPHIPointers();

    if( not( gslLoaded ) or not( libAPHILoaded ) ) then
      freeLoadedPointers()
    ;

end.

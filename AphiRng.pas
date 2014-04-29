{*
AphiRng.pas
-----------
Begin: 2008/11/10
Last revision: $Date: 2010-10-14 17:30:43 $ $Author: rhupalo $
Version: $Revision: 1.9 $
Project: APHI Delphi Library for Simulation Modeling
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2008, 2010 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

  This unit contains functions for random number generation.  The functions
  use a version of the SPRNG library compiled for Windows: see
  http://sprng.cs.fsu.edu and http://www.naadsm.org/opensource/sprng.
  
  The random number routines rngRand() and rngRandInt() can be used directly.
  The random number generator is also be used with a variety of routines for
  random number generation from probability density functions: see unit
  ProbDensityFunctions.pas for more information.
  
  This unit requires the following libraries:
    - APHI General Purpose Delphi Library, version 1.3.1.
      http://www.naadsm.org/opensource/delphi/ 
    
    - APHI Delphi Library for Stochastic Modeling, version 0.2.
      http://www.naadsm.org/opensource/libaphi/
    
    - SPRNG for Microsoft Windows, version 2.0a.
      http://www.naadsm.org/opensource/sprng/
      
    - The GNU Scientific Library for Windows version 1.6 or higher.
      http://gnuwin32.sourceforge.net/packages/gsl.htm
}

  (*
    Documentation generation tags begin with {* or ///
    Replacing these with (* or // foils the documentation generator
  *)

unit AphiRng;

interface

  // Unit/library management
  //------------------------

  // True if the required libraries were loaded and can be used.  Otherwise false.
  function libAphiRngLoaded(): boolean;

  // Empty string if the required libraries were loaded and can be used.  Otherwise, contains an error message.
  function libAphiRngLoadErrors(): string;

  // Routines for random number generation
  //--------------------------------------

  // Returns a random double in [0,1], i.e., 0 <= r < 1
  function rngRand(): double;

  // Returns a random integer in [min,maxPlusOne], i.e., min <= r <= max
  function rngRandInt( min, maxPlusOne: integer ): integer;


  // RNG seed
  //---------
  { Reinitializes the SPRNG random number generator with the indicated seed.
  Seed will be generated automatically if -1 is passed as the parameter value. }
  procedure setRngSeed( const seed: integer = -1 );

  /// Returns the current seed.
  function rngSeed(): integer;


  // Specialized functions for the RNG libraries
  //--------------------------------------------
  /// Returns a pointer to the GSL-compatible RNG.  Used by ProbDensityFunctions.
  function gslRngPtr(): cardinal;

  // Returns a pointer to the RNG. Used by ProbDensityFunctions.
  function rngPtr(): cardinal;

  // GMP functions
  //--------------
  // These two functions don't have anything to do with random number generation,
  // but this is a convenient place to test the multiple precision arithmetic functions
  // recently introduced in libAPHI.
  // FIX ME: if GMP functions are ever used for real work, these functions should be
  // moved to a more appropriate place.
  function gmpFactorial( val: cardinal ): cardinal;
  function gmpFactorialAsStr( val: cardinal ): string;


implementation

  uses
    // Standard Delphi units
    Windows, // Defines THandle
    Math,
    SysUtils,

    // APHI General Purpose Delphi Library
    DebugWindow,
    FunctionPointers,
    MyStrUtils
  ;

  const
    REQUIRED_APHI_DLL_VERSION: string = '0.5';  /// required version of libaphi.dll
    DBSHOWMSG: boolean = true; /// Set to true to display debugging messages for this unit.

  var
    _rngPtr: cardinal;              /// pointer to the RNG library
    _gslRngPtr: cardinal;           /// pointer to the GSL-compatible RNG library
    _libAphiRngLoaded: boolean;     /// indicates whether the required libraries were loaded
    _libAphiRngLoadErrors: string;  /// error messages generated if loading the required libraries fails

    aphi_rng_new_generator: TCFnInt_1_Int; //RAN_gen_t* RAN_new_generator( int seed );
    aphi_rng_num: TCFnDouble_1_Int; // double RAN_num( RAN_gen_t* );
    aphi_rng_generator_as_gsl: TCFnInt_1_Int; // gsl_rng* RAN_generator_as_gsl( RAN_gen_t* );
    aphi_rng_free_generator: TCFnVoid_1_Int; // void RAN_free_generator( RAN_gen_t* );
    aphi_rng_seed: TCFnInt_0; // int aphi_rng_seed( void );

    // FIX ME: if GMP functions are ever used for real work, these functions should be
    // moved to a more appropriate place.
    test_factorial: TCFnCardinal_1_Cardinal;
    test_factorial_as_str: TCFnCharP_1_Cardinal;

  {*
    Indicates whether the required libraries were loaded and can be used.
    @return true if libraries loaded sucessfully else false
  }
  function libAphiRngLoaded(): boolean;
    begin
      result := _libAphiRngLoaded;
    end
  ;

  {*
    Provides error messages generated when loading the required libraries is attempted.
    @return Empty string if the required libraries are available. Otherwise, contains an error message.
  }
  function libAphiRngLoadErrors(): string;
    begin
      result := _libAphiRngLoadErrors;
    end
  ;


  {*
     Returns a pointer to the RNG library.
     @return  pointer, note it is a cardinal data type, not PCardinal.
     @comment Used by ProbDensityFunctions.
  }
  function rngPtr(): cardinal;
    begin
      result := _rngPtr;
    end
  ;


  {*
     Returns a pointer to the GSL-compatible RNG library.
     @return  pointer, note it is a cardinal data type, not PCardinal.
     @comment Used by ProbDensityFunctions.
  }
  function gslRngPtr(): cardinal;
    begin
      result := _gslRngPtr;
    end
  ;
    

  {*
    Reinitializes the SPRNG random number generator with seed.
    @param seed whole number, if value is -1 then seed will be generated automatically
  }
  procedure setRngSeed( const seed: integer = -1 );
    begin
      if( _libAphiRngLoaded ) then
        begin
          aphi_rng_free_generator( _rngPtr );
          _rngPtr := aphi_rng_new_generator( seed );
          _gslRngPtr := aphi_rng_generator_as_gsl( _rngPtr );
        end
      ;
    end
  ;

  {*
    Returns the current seed value.
    @return the value of the current seed
  }
  function rngSeed(): integer;
    begin
      if( _libAphiRngLoaded ) then
        result := aphi_rng_seed()
      else
        result := -1
      ;
    end
  ;
  
  {*
    Generates a random number between 0 and 1
    @return a floating-point number (r) where 0 <= r < 1
  }
  function rngRand(): double;
    begin
      if( _libAphiRngLoaded ) then
        result := aphi_rng_num( _rngPtr )
      else
        begin
          raise exception.Create( 'Library not loaded in rngRand()' );
          result := NaN;
        end
      ;
    end
  ;

  {*
    Generates a random number between min and maxPlusOne - 1
    @param min lower bound of random number generation range
    @param maxPlusOne one more than the upper bound of random number range
    @return a whole number (r) where min <= r <= max
    @comment If min < 0 the input range is preserved but adjusted upwards making min 0.
  }
  function rngRandInt( min, maxPlusOne: integer ): integer;
    var
      range: integer;
      d: double;
      adj: integer;
    begin
      if( not _libAphiRngLoaded ) then
        begin
          raise exception.create( 'Library not loaded in rngRandInt()' );
          result := 0;
        end
      else if( min > maxPlusOne ) then
        begin
          raise exception.create( 'min > maxPlusOne in rngRandInt()' );
          result := 0;
        end
      else if( min = maxPlusOne ) then
        result := min
      else
        begin
          if( 0 > min ) then
            begin
              adj := -1 * min;
              min := min + adj;
              maxPlusOne := maxPlusOne + adj;

              //dbcout2( 'Adjusted range: ' + intToStr( min ) + ' to ' + intToStr( maxPlusOne ) );
            end
          else
            adj := 0
          ;

          range := maxPlusOne - min;
          d := rngRand();
          result := trunc( d * range + min ) - adj;
        end
      ;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// GMP functions
//-----------------------------------------------------------------------------
  // FIX ME: if GMP functions are ever used for real work, these functions should be
  // moved to a more appropriate place.

  {*
     Function has anything to do with random number generation, but this is a
     convenient place to test the multiple precision arithmetic functions introduced in libAPHI.
     @param val positive number for which the factorial is calculated
     @return  the product of all positive integers less than or equal to val
  }
  function gmpFactorial( val: cardinal ): cardinal;
    begin
      if( _libAphiRngLoaded ) then
        result := test_factorial( val )
      else
        begin
          raise exception.Create( 'Library not loaded in gmpFactorial()' );
          result := 0;
        end
      ;
    end
  ;


  {*
     Function has anything to do with random number generation, but this is a
     convenient place to test the multiple precision arithmetic functions introduced in libAPHI.
     @param val positive number for which the factorial is calculated
     @return  text string of the product of all positive integers less than or equal to val
  }
  function gmpFactorialAsStr( val: cardinal ): string;
    begin
      if( _libAphiRngLoaded ) then
        result := test_factorial_as_str( val )
      else
        begin
          raise exception.Create( 'Library not loaded in gmpFactorialAsStr()' );
          result := '';
        end
      ;
    end
  ;
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------
// DLL handling
//-----------------------------------------------------------------------------
  /// Sets the function pointers to nil in the event that the library does not sucessfully load.
  procedure freeLoadedPointers();
    begin
      aphi_rng_new_generator := nil;
      aphi_rng_num := nil;
      aphi_rng_generator_as_gsl := nil;
      aphi_rng_free_generator := nil;
      aphi_rng_seed := nil;

    // FIX ME: if GMP functions are ever used for real work, these functions should be
    // moved to a more appropriate place.
    test_factorial := nil;
    test_factorial_as_str := nil;
    end
  ;


  {*
    Does all the heavy lifting to load the dll and assign the function pointers.
    @result true if the library successfully loads, else false
  }
  function loadAPHIPointers(): boolean;
    var
      rngHandle: THandle; //Handle used to open the DLL.  Defined in unit Windows.

      lib_version: function(): pchar; cdecl;
      libVersion: string;
      lastError: integer;
    begin
      _libAphiRngLoadErrors := '';

      try
        setLastError( 0 );
        rngHandle := loadLibrary( 'libaphi.dll' );
      except
        result := false;
        lastError := getLastError();
        appendToString( _libAphiRngLoadErrors, 'Could not load library libaphi.dll: system error code ' + intToStr( lastError ) );
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
              appendToString( _libAphiRngLoadErrors, 'MISSING FUNCTION lib_version' );
              result := false;
              exit;
            end
          else
            begin
              libVersion := lib_version();
              if( REQUIRED_APHI_DLL_VERSION <> libVersion ) then
                begin
                  appendToString( _libAphiRngLoadErrors, 'INCOMPATIBLE DLL VERSION: ' + libVersion );
                  result := false;
                end
              ;
            end
          ;

          // Creating the generator
          //-----------------------
          aphi_rng_new_generator := getProcAddress( rngHandle, 'aphi_rng_new_generator' );
          if( nil = @aphi_rng_new_generator ) then
            begin
              appendToString( _libAphiRngLoadErrors, 'MISSING FUNCTION aphi_rng_new_generator' );
              result := false;
            end
          ;

          // Using the generator
          //--------------------
          aphi_rng_num := getProcAddress( rngHandle, 'aphi_rng_num' );
          if( nil = @aphi_rng_num ) then
            begin
              appendToString( _libAphiRngLoadErrors, 'MISSING FUNCTION aphi_rng_num' );
              result := false;
            end
          ;
          
          // Configuring the GSL-compatible generator
          //-----------------------------------------
          aphi_rng_generator_as_gsl := getProcAddress( rngHandle, 'aphi_rng_generator_as_gsl' );
          if( nil = @aphi_rng_generator_as_gsl ) then
            begin
              appendToString( _libAphiRngLoadErrors, 'MISSING FUNCTION aphi_rng_generator_as_gsl' );
              result := false;
            end
          ;          

          // Retrieving the RNG seed
          //------------------------
          aphi_rng_seed := getProcAddress( rngHandle, 'aphi_rng_seed' );
          if( nil = @aphi_rng_seed ) then
            begin
              appendToString( _libAphiRngLoadErrors, 'MISSING FUNCTION aphi_rng_seed' );
              result := false;
            end
          ;

          // Freeing the generator
          //-----------------------
          aphi_rng_free_generator := getProcAddress( rngHandle, 'aphi_rng_free_generator' );
          if( nil = @aphi_rng_free_generator ) then
            begin
              appendToString( _libAphiRngLoadErrors, 'MISSING FUNCTION aphi_rng_free_generator' );
              result := false;
            end
          ;

          // GMP functions
          //--------------
          // FIX ME: if GMP functions are ever used for real work, these functions should be
          // moved to a more appropriate place.
          test_factorial := getProcAddress( rngHandle, 'test_factorial' );
          if( nil = @test_factorial ) then
            begin
              appendToString( _libAphiRngLoadErrors, 'MISSING FUNCTION test_factorial' );
              result := false;
            end
          ;

          test_factorial_as_str := getProcAddress( rngHandle, 'test_factorial_as_str' );
          if( nil = @test_factorial_as_str ) then
            begin
              appendToString( _libAphiRngLoadErrors, 'MISSING FUNCTION test_factorial_as_str' );
              result := false;
            end
          ;
        end
      else
        begin
          result := false;
          lastError := getLastError();
          appendToString( _libAphiRngLoadErrors, 'Could not load library libaphi.dll: system error code ' + intToStr( lastError ) );
        end
      ;

      dbcout( 'APHI rng loaded', DBSHOWMSG );
      dbcout( result, DBSHOWMSG );
    end
  ;

initialization
  // Load the pointers
  _rngPtr := 0;
  _gslRngPtr := 0;

  _libAphiRngLoaded := loadAphiPointers();
    
  if( _libAphiRngLoaded ) then
    begin
      // Create the generator (it can be reinitialized later with a new seed)
      _rngPtr := aphi_rng_new_generator( -1 );
      
      // Set up the GSL-compatible generator
      _gslRngPtr := aphi_rng_generator_as_gsl( _rngPtr ); 
    end
  ;

finalization
  
  if( _libAphiRngLoaded ) then
    begin
      // Free the generator
      dbcout2( 'Freeing the generator...' );
      aphi_rng_free_generator( _rngPtr );
      dbcout2( 'Done.' );
    end
  else
    dbcout2( 'Something didn''t work.' )
  ;


end.




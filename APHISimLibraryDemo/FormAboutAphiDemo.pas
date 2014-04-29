unit FormAboutAphiDemo;

(*
FormAboutAphiDemo.pas/dfm
-------------------------
Begin: 2008/12/14
Last revision: $Date: 2010-05-20 17:32:53 $ $Author: areeves $
Version: $Revision: 1.4 $
Project: APHI Delphi Library for Simulation Modeling: Demo application
Website: http://www.naadsm.org/opensource/libaphi/
Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
--------------------------------------------------
Copyright (C) 2008 - 2010 Animal Population Health Institute, Colorado State University

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
*)

interface

  uses
    // Standard Delphi units
    Windows, 
    Messages, 
    SysUtils, 
    Variants, 
    Classes, 
    Graphics, 
    Controls, 
    Forms,
    Dialogs,
    StdCtrls,
    ExtCtrls
  ;

  type TFormAboutAphiDemo = class(TForm)
      lblDescription: TLabel;
      lblWebsite: TLabel;
      lblCopyright: TLabel;
      lblLicenseBlurb: TLabel;
      pnlButtons: TPanel;
      btnLicense: TButton;
      btnOK: TButton;

      procedure btnOKClick(Sender: TObject);
      procedure btnLicenseClick(Sender: TObject);
      procedure lblWebsiteClick(Sender: TObject);
      procedure FormCreate(Sender: TObject);

    protected
      // Construction/initialization/destruction
      procedure translateUI();

    public
      // Construction/initialization/destruction
      constructor create( AOwner: TComponent ); override;

    end
  ;

implementation

{$R *.dfm}

  uses
    // Standard Delphi units
    ShellAPI,

    // APHI General-Purpose Delphi Library
    MyStrUtils,
    DialogLongMessage,
    I88n
  ;

  constructor TFormAboutAphiDemo.create( AOwner: TComponent );
    begin
      inherited create( AOwner );
      translateUI();
    end
  ;


  procedure TFormAboutAphiDemo.FormCreate( Sender: TObject );
    begin
      Assert(not Scaled, 'You should set Scaled property of Form to False!');

      if( Screen.PixelsPerInch <> 96 ) then
        ScaleBy( Screen.PixelsPerInch, 96 )
      ;
    end
  ;


  procedure TFormAboutAphiDemo.translateUI();
    var
      str1, str2, str3: string;
    begin
      str1 := ''
        + 'This program is free software; you can redistribute it and/or modify'
        + ' it under the terms of the GNU General Public License as published by the'
        + ' Free Software Foundation; either version 2 of the License, or (at your option)'
        + ' any later version.  For complete license details, please see below.'
      ;

      self.Caption := tr( 'About the demo application' );
      btnOK.Caption := tr( '&OK' );
      btnLicense.Caption := tr( '&License...' );
      lblCopyright.Caption := tr( 'Copyright © 2008 - 2009 Animal Population Health Institute at Colorado State University' );
      lblLicenseBlurb.Caption := tr( str1 );

      str2 := ''
        + 'This application demonstrates how to use elements included in the APHI'
        + ' Delphi Library for Simulation Modeling.  This library, a collection of'
        + ' units and visual components written in the Delphi programming language,'
        + ' provides a set of probability density functions, random number generation'
        + ' capabilities, and a framework for the creation of graphical applications'
        + ' for stochastic simulation modeling.  The library is released under the'
        + ' terms of the GNU General Public license (see below).'
      ;

      str3 := ''
        + 'Source code for the APHI Delphi Library for Simulation Modeling is'
        + ' available from the following website:'
      ;

      lblDescription.Caption := tr( str2 ) + endl + endl + tr( str3 );
    end
  ;


  procedure TFormAboutAphiDemo.btnOKClick( Sender: TObject );
    begin
      self.Close();
    end
  ;


  procedure TFormAboutAphiDemo.btnLicenseClick( Sender: TObject );
    var
      frm: TDialogLongMessage;
    begin
      frm := TDialogLongMessage.create( self, tr( 'GNU General Public License' ) );
      frm.mmoLongMessage.Font.Name := 'Courier New';
      frm.setMessage( i88nLicense() );

      self.Hide();

      frm.ShowModal();
      frm.Free();

      self.Show();
    end
  ;


  procedure TFormAboutAphiDemo.lblWebsiteClick(Sender: TObject);
    begin
        ShellExecute(
          Application.Handle,
          PChar( 'open' ),
          PChar( lblWebsite.Caption ),
          PChar( 0 ),
          nil,
          SW_NORMAL
        );
    end
  ;



end.
 
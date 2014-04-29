object FormAboutAphiDemo: TFormAboutAphiDemo
  Left = 514
  Top = 102
  BorderStyle = bsDialog
  Caption = 'About the demo application'
  ClientHeight = 322
  ClientWidth = 448
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  Position = poOwnerFormCenter
  Scaled = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object lblDescription: TLabel
    Left = 16
    Top = 9
    Width = 417
    Height = 125
    Caption = 'lblDescription'
    Constraints.MaxHeight = 125
    Constraints.MaxWidth = 417
    Constraints.MinHeight = 125
    Constraints.MinWidth = 417
    WordWrap = True
  end
  object lblWebsite: TLabel
    Left = 16
    Top = 153
    Width = 214
    Height = 13
    Cursor = crHandPoint
    Caption = 'http://www.naadsm.org/opensource/libaphi/'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clBlue
    Font.Height = -11
    Font.Name = 'MS Sans Serif'
    Font.Style = [fsUnderline]
    ParentFont = False
    OnClick = lblWebsiteClick
  end
  object lblCopyright: TLabel
    Left = 16
    Top = 177
    Width = 425
    Height = 26
    Caption = 
      'Copyright '#169' 2008 - 2010 Animal Population Health Institute at Co' +
      'lorado State University'
    WordWrap = True
  end
  object lblLicenseBlurb: TLabel
    Left = 16
    Top = 208
    Width = 417
    Height = 61
    Caption = 
      'This program is free software; you can redistribute it and/or mo' +
      'dify it under the terms of the GNU General Public License as pub' +
      'lished by the Free Software Foundation; either version 2 of the ' +
      'License, or (at your option) any later version.  For complete li' +
      'cense details, please see below.'
    Constraints.MaxHeight = 61
    Constraints.MaxWidth = 417
    Constraints.MinHeight = 61
    Constraints.MinWidth = 417
    WordWrap = True
  end
  object pnlButtons: TPanel
    Left = 0
    Top = 284
    Width = 448
    Height = 38
    Align = alBottom
    BevelOuter = bvNone
    TabOrder = 0
    object btnLicense: TButton
      Left = 266
      Top = 5
      Width = 89
      Height = 25
      Caption = '&License...'
      TabOrder = 0
      OnClick = btnLicenseClick
    end
    object btnOK: TButton
      Left = 365
      Top = 5
      Width = 75
      Height = 25
      Caption = '&OK'
      Default = True
      TabOrder = 1
      OnClick = btnOKClick
    end
  end
end

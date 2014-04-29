object FormMain: TFormMain
  Left = 86
  Top = 21
  BorderStyle = bsDialog
  Caption = 'APHI Delphi Library for Simulation Modeling'
  ClientHeight = 234
  ClientWidth = 285
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  Position = poScreenCenter
  PixelsPerInch = 96
  TextHeight = 13
  object Panel1: TPanel
    Left = 0
    Top = 193
    Width = 285
    Height = 41
    Align = alBottom
    TabOrder = 0
    object btnTest: TButton
      Left = 112
      Top = 8
      Width = 75
      Height = 25
      Caption = 'Test RNG'
      TabOrder = 0
      OnClick = btnTestClick
    end
  end
  object Panel2: TPanel
    Left = 0
    Top = 0
    Width = 285
    Height = 168
    Align = alClient
    TabOrder = 1
    object lblUrl: TLabel
      Left = 1
      Top = 154
      Width = 283
      Height = 13
      Cursor = crHandPoint
      Align = alBottom
      Alignment = taCenter
      Caption = 'http://www.naadsm.org/opensource/libaphi/'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clBlue
      Font.Height = -11
      Font.Name = 'MS Sans Serif'
      Font.Style = [fsUnderline]
      ParentFont = False
      WordWrap = True
      OnClick = lblUrlClick
    end
    object lblDescr: TLabel
      Left = 11
      Top = 11
      Width = 263
      Height = 143
      Align = alClient
      Caption = 'lblDescr'
      WordWrap = True
    end
    object Panel3: TPanel
      Left = 1
      Top = 11
      Width = 10
      Height = 143
      Align = alLeft
      TabOrder = 0
    end
    object Panel4: TPanel
      Left = 274
      Top = 11
      Width = 10
      Height = 143
      Align = alRight
      TabOrder = 1
    end
    object Panel5: TPanel
      Left = 1
      Top = 1
      Width = 283
      Height = 10
      Align = alTop
      TabOrder = 2
    end
  end
  object Panel6: TPanel
    Left = 0
    Top = 168
    Width = 285
    Height = 25
    Align = alBottom
    TabOrder = 2
  end
end

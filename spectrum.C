// developed by Vitalii Petrov, vitalii17@bk.ru
//  based on examples from https://root.cern.ch/root/htmldoc/guides/users-guide/WritingGUI.html
#include <TApplication.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <fstream>
#include <vector>
#include <TGDoubleSlider.h>
#include "TGTextEntry.h"
#include <TGSlider.h>
#include <TGNumberEntry.h>
#include <TGraphErrors.h>
#include <TGLayout.h>
#include <TGWindow.h>
#include <TGLabel.h>
#include <TString.h>
#include <TGButtonGroup.h>
#include <TGFileBrowser.h>
#include <TSystem.h>
#include <utility>
#include <map>
#include <algorithm>

EColor color[] = {kGreen, kCyan, kOrange, kMagenta, kBlack, kBlue, kGreen, kYellow, kMagenta, kOrange};
struct PeakParams
{
  vector<int> MaxBins;
  vector<int> nBins;
  int UsedBins;
};

void FindPeak(TH1D *hist, struct PeakParams *PPar)
{
  int FirstBin = PPar->UsedBins;
  int PeakBinsCounter = 0;
  int nBins = 0;
  double LocalMax = 0;
  int LocalMaxBin;
  while (PeakBinsCounter < 2)
  {
    if (FirstBin + nBins >= hist->GetNbinsX())
      return;
    if (hist->GetBinContent(FirstBin + nBins) > LocalMax)
    {
      PeakBinsCounter = 0;
      LocalMax = hist->GetBinContent(FirstBin + nBins);
      LocalMaxBin = FirstBin + nBins;
    }
    else
      PeakBinsCounter++;
    nBins++;
  }
  int nIncreasingBins = nBins - 2;
  PeakBinsCounter = 0;
  double LocalMin = hist->GetBinContent(FirstBin + nBins - 1);
  while (PeakBinsCounter < 1)
  {
    if (hist->GetBinContent(FirstBin + nBins) < LocalMin)
    {
      PeakBinsCounter = 0;
      LocalMin = hist->GetBinContent(FirstBin + nBins);
    }
    else
      PeakBinsCounter++;
    nBins++;
  }
  int nDecreasingBins = nBins - nIncreasingBins;
  int nBinsForHalfPeak = min(nDecreasingBins, nIncreasingBins);
  if (nBins == 0)
    nBins = 1;
  PPar->UsedBins = FirstBin + nBins;
  if (nBinsForHalfPeak > 2 && LocalMax - hist->GetBinContent(FirstBin) > 30)
  {
    PPar->MaxBins.push_back(LocalMaxBin);
    PPar->nBins.push_back((nDecreasingBins + nIncreasingBins) / 2);
  }
  return;
}

double FindClosestPeak(TH1D *hist, double ExpectedBin)
{
  int nBins = 0;
  Double_t X[10], Y[10];
  for (int iBin = 0; iBin < 10; iBin++)
  {
    if (ExpectedBin - 10 / 2 + iBin < 0)
      continue;
    X[nBins] = hist->GetBinCenter(1 + ExpectedBin - 10 / 2 + iBin);
    Y[nBins] = hist->GetBinContent(1 + ExpectedBin - 10 / 2 + iBin);
    nBins++;
  }
  int max = distance(Y, max_element(Y, Y + nBins));
  return X[max];
}

enum ETestCommandIdentifiers
{
  None,
  LoadSpect,
  LoadBckg,
  LoadCalib,
  LoadEff,
  ShowSource,
  ShowBckg,
  ShowSourceBckg,
  ShowCalib,
  ShowEff,
  SaveImage,
  SaveCalib,
  SaveEffic,
};
class MyMainFrame : public TGMainFrame
{
  RQ_OBJECT("MyMainFrame")
private:
  TGDoubleHSlider *fSlider0, *fSlider;
  TGNumberEntry *fNumber, *fNumberEff, *fnGamma, *fBorderLength, *fHistBinCombined, *fNumberOfPeacks, *fMeasureTime;
  TGHorizontalFrame *Frame1, *Frame2, *Frame3;
  TGLabel *fLabelPosition, *fLabelWidth, *fLabelIntegral, *fLabel2, *fLabelB;
  TGraphErrors *fGraph, *fGraphEff;
  TF1 *FitEff;
  TGTextButton *fEnergyChannel, *fCalibration, *fEff, *fLogy;
  TGMenuBar *fMenuLoad, *fMenuShow, *fMenuSave, *fMenuXaxisDim;
  TGLayoutHints *fMyBarLayout;
  TRootEmbeddedCanvas *fEcanvas;
  TCanvas *cAuto;
  Bool_t EnergyNotChannel, EnergyNotChannelOneTime, CalibrationOrNot, Eff, Logy, GraphEffChanged, EffSliderFirstMoved, Rebinned;

public:
  TGTransientFrame *fMain;
  MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~MyMainFrame();
  void ApplyChanges();
  void LoadSpectrFile(int ForBckg);
  void LoadCalibrationFile();
  void LoadEffFile();
  void DoDrawSource();
  void DoDrawBckg();
  void DoDrawSourceWithoutBckg();
  void ChangeSlider();
  void SetXaxisDimAs(Int_t id);
  void SetLoadAs(Int_t id);
  void SetShowAs(Int_t id);
  void SetSaveAs(Int_t id);
  void SetSpectrumAs(Int_t id);
  void DoDraw();
  void Rebin();
  void DoDrawCalibrationGraph();
  void DoDrawEff();
  void Save();
  void SaveCal();
  void SaveEff();
  void AddPoint();
  void AddPointEff();
  void ChangeEnergyChannelLabel();
  void ChangeLogyLabel();
  void ChangeEffLabel();
  TH1D *hLoaded, *hSignal, *hBackground, *TempHist;
  vector<double> PeakPosition, PeakPositionUncertancy, PeakIntegral, PeakIntegralUncertancy;
  Double_t a, b, FitEffError;                                                                                     // y = ax + b
  Double_t SliderPositions[6], Slider0Positions[6];                                                               // y = ax + b
  int iPeak, iPeak2, NumberOfPeacks, istring, HistBinCombined0, HistBinCombined, AutoCanvasCreated, WasOpened[4]; // 0 hLoaded, 1 hSignal, 2 hBackground, 3 Calibration
  int nPeaksToAddCal, nPeaksToAddEff;
  Int_t NowOpen;
};

void MyMainFrame::ApplyChanges()
{
  if (GraphEffChanged)
  {
    int nPoints = fGraphEff->GetN();
    double x, y, xerr, yerr;
    vector<double> xArr, yArr, xerrArr, yerrArr;
    for (int i = 0; i < nPoints; i++)
    {
      fGraphEff->GetPoint(i, x, y);
      yerr = fGraphEff->GetErrorY(i);
      xArr.push_back(x);
      yArr.push_back(y);
      yerrArr.push_back(yerr);
    }
    fGraphEff->Set(0);
    for (int i = 0; i < nPoints; i++)
    {
      int ind = min_element(xArr.begin(), xArr.end()) - xArr.begin();
      fGraphEff->SetPoint(i, (xArr[ind] - b) / a, yArr[ind]);
      fGraphEff->SetPointError(i, 0, yerrArr[ind]);
      xArr.erase(xArr.begin() + ind);
      yArr.erase(yArr.begin() + ind);
      yerrArr.erase(yerrArr.begin() + ind);
    }
    GraphEffChanged = kFALSE;
  }
}

void MyMainFrame::ChangeEnergyChannelLabel()
{
  if (!WasOpened[3])
  {
    cout << "Warning: No calibration applied!\n";
    return;
  }
  fEnergyChannel->SetState(kButtonDown);
  if (!EnergyNotChannel)
  {
    fEnergyChannel->SetText("&Disable Calibation");
    EnergyNotChannel = kTRUE;
    EnergyNotChannelOneTime = kTRUE;
  }
  else
  {
    fEnergyChannel->SetText("&Enable Calibation");
    EnergyNotChannel = kFALSE;
    EnergyNotChannelOneTime = kFALSE;

    TempHist = new TH1D("TempHist", "", hLoaded->GetNbinsX(), (hLoaded->GetXaxis()->GetXmin() - b) / a, (hLoaded->GetXaxis()->GetXmax() - b) / a);
    for (int ibin = 0; ibin < hLoaded->GetNbinsX(); ibin++)
      TempHist->SetBinContent(ibin, hLoaded->GetBinContent(ibin));
    delete hLoaded;
    hLoaded = (TH1D *)TempHist->Clone();
    delete TempHist;
    EnergyNotChannelOneTime = kFALSE;
  }
  fEnergyChannel->SetState(kButtonUp);
  SetShowAs(NowOpen);
}

void MyMainFrame::ChangeEffLabel()
{
  fEff->SetState(kButtonDown);
  if (!Eff)
  {
    fEff->SetText("&No Efficiency");
    fnGamma->SetNumber(1.00);
    Eff = kTRUE;
  }
  else
  {
    fEff->SetText("&Consider Efficiency");
    Eff = kFALSE;
  }
  fEff->SetState(kButtonUp);
  SetShowAs(NowOpen);
}
void MyMainFrame::ChangeLogyLabel()
{
  fLogy->SetState(kButtonDown);
  if (!Logy)
  {
    fLogy->SetText("&NoLog");
    Logy = kTRUE;
  }
  else
  {
    fLogy->SetText("&LogY");
    Logy = kFALSE;
  }
  fLogy->SetState(kButtonUp);
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h)
    : TGMainFrame(p, w, h)
{
  for (int i = 0; i < 4; i++)
  {
    WasOpened[i] = 0;
  }
  NowOpen = None;
  AutoCanvasCreated = 0;
  // Create a main frame
  EnergyNotChannelOneTime = kFALSE;
  EnergyNotChannel = kFALSE;
  CalibrationOrNot = kTRUE;
  Eff = kFALSE;
  Logy = kFALSE;
  GraphEffChanged = kFALSE;
  EffSliderFirstMoved = kFALSE;

  Frame1 = new TGHorizontalFrame(this, 1500, 300);
  Frame2 = new TGHorizontalFrame(this, 1500, 300);
  Frame3 = new TGHorizontalFrame(this, 1500, 300);
  fMyBarLayout = new TGLayoutHints(kLHintsLeft, 5, 4, 3, 0);

  fMenuLoad = new TGMenuBar(Frame1, 1, 1, kHorizontalFrame);
  fMenuShow = new TGMenuBar(Frame1, 1, 1, kHorizontalFrame);
  fMenuSave = new TGMenuBar(Frame1, 1, 1, kHorizontalFrame);
  fMenuXaxisDim = new TGMenuBar(Frame1, 1, 1, kHorizontalFrame);

  TGPopupMenu *fSetLoad = new TGPopupMenu(gClient->GetRoot());
  fSetLoad->AddEntry("&Source Spectrum", LoadSpect);
  fSetLoad->AddSeparator();
  fSetLoad->AddEntry("&Background Spectrum", LoadBckg);
  fSetLoad->AddSeparator();
  fSetLoad->AddEntry("&Calibation", LoadCalib);
  fSetLoad->AddSeparator();
  fSetLoad->AddEntry("&Efficiency", LoadEff);
  fSetLoad->Connect("Activated(Int_t)", "MyMainFrame", this,
                    "SetLoadAs(Int_t)");
  fMenuLoad->AddPopup("&Load |", fSetLoad, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));

  TGPopupMenu *fSetShow = new TGPopupMenu(gClient->GetRoot());
  fSetShow->AddEntry("&Source", ShowSource);
  fSetShow->AddSeparator();
  fSetShow->AddEntry("&Background", ShowBckg);
  fSetShow->AddSeparator();
  fSetShow->AddEntry("&Source - Bckg", ShowSourceBckg);
  fSetShow->AddSeparator();
  fSetShow->AddEntry("&Calibation", ShowCalib);
  fSetShow->AddSeparator();
  fSetShow->AddEntry("&Efficiency", ShowEff);
  fSetShow->Connect("Activated(Int_t)", "MyMainFrame", this,
                    "SetShowAs(Int_t)");
  fMenuShow->AddPopup("&Show |", fSetShow, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));

  TGPopupMenu *fSetSave = new TGPopupMenu(gClient->GetRoot());
  fSetSave->AddEntry("&Current image", SaveImage);
  fSetSave->AddSeparator();
  fSetSave->AddEntry("&Calibation data", SaveCalib);
  fSetSave->AddSeparator();
  fSetSave->AddEntry("&Efficiency data", SaveEffic);
  fSetSave->Connect("Activated(Int_t)", "MyMainFrame", this,
                    "SetSaveAs(Int_t)");
  fMenuSave->AddPopup("&Save |", fSetSave, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));

  TGLabel *fTextLabelBinning = new TGLabel(Frame1, "  Combine channels:");
  HistBinCombined0 = 1;
  HistBinCombined = 1;
  fHistBinCombined = new TGNumberEntry(Frame1, 1, 10, 100, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMinMax, 1, 10000);
  fHistBinCombined->Connect("ValueSet(int)", "MyMainFrame", this, "Rebin()");
  (fHistBinCombined->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "Rebin()");
  TGLabel *fTextLabelBackgroundFit = new TGLabel(Frame1, "  N Bins for background fit:");
  fBorderLength = new TGNumberEntry(Frame1, 0, 10, 100, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 1, 10000);
  fBorderLength->Connect("ValueSet(Int_t)", "MyMainFrame", this, "DoDraw()");
  (fBorderLength->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "DoDraw()");
  fBorderLength->SetIntNumber(10);
  TGLabel *fTextLabelNumberOfPeacks = new TGLabel(Frame1, "  Number of Peaks:");
  fNumberOfPeacks = new TGNumberEntry(Frame1, 0, 10, 100, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 1, 10);
  fNumberOfPeacks->Connect("ValueSet(Int_t)", "MyMainFrame", this, "DoDraw()");
  (fNumberOfPeacks->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "DoDraw()");
  fNumberOfPeacks->SetIntNumber(1);
  TGLabel *fTextLabelMeasureTime = new TGLabel(Frame1, "   Time of Measurement (s):");
  fMeasureTime = new TGNumberEntry(Frame1, 0, 10, 100, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 1, 1e7);
  fMeasureTime->Connect("ValueSet(Int_t)", "MyMainFrame", this, "DoDraw()");
  (fMeasureTime->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "DoDraw()");
  fMeasureTime->SetIntNumber(1);

  TGTextButton *exit = new TGTextButton(Frame1, "&Exit", "gApplication->Terminate(0)");

  Frame1->AddFrame(fMenuLoad, fMyBarLayout);
  Frame1->AddFrame(fMenuShow, fMyBarLayout);
  Frame1->AddFrame(fMenuSave, fMyBarLayout);
  Frame1->AddFrame(fTextLabelBinning, new TGLayoutHints(kLHintsLeft, 50, 5, 8, 0));
  Frame1->AddFrame(fHistBinCombined, fMyBarLayout);
  Frame1->AddFrame(fTextLabelBackgroundFit, new TGLayoutHints(kLHintsLeft, 0, 5, 8, 0));
  Frame1->AddFrame(fBorderLength, fMyBarLayout);
  Frame1->AddFrame(fTextLabelNumberOfPeacks, new TGLayoutHints(kLHintsLeft, 0, 5, 8, 0));
  Frame1->AddFrame(fNumberOfPeacks, fMyBarLayout);
  Frame1->AddFrame(fTextLabelMeasureTime, new TGLayoutHints(kLHintsLeft, 0, 5, 8, 0));
  Frame1->AddFrame(fMeasureTime, fMyBarLayout);
  Frame1->AddFrame(exit, fMyBarLayout);

  TGLabel *fTextLabelEnergy = new TGLabel(Frame2, "Energy Calibration.       Insert Energy:");
  fNumber = new TGNumberEntry(Frame2, 0, 10, 999, TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 99999);
  fNumber->Connect("ValueSet(Float_t)", "MyMainFrame", this, "AddPoint()");
  (fNumber->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "AddPoint()");
  TGTextButton *AddPointButton = new TGTextButton(Frame2, "&Add Calibration point");
  AddPointButton->Connect("Clicked()", "MyMainFrame", this, "AddPoint()");
  fEnergyChannel = new TGTextButton(Frame2, "&Enable Calibation");
  fEnergyChannel->Connect("Clicked()", "MyMainFrame", this, "ChangeEnergyChannelLabel()");

  Frame2->AddFrame(fTextLabelEnergy, new TGLayoutHints(kLHintsLeft, 10, 5, 3, 0));
  Frame2->AddFrame(fNumber, new TGLayoutHints(kLHintsLeft, 0, 5, 3, 0));
  Frame2->AddFrame(AddPointButton, new TGLayoutHints(kLHintsLeft, 102, 5, 3, 0));
  Frame2->AddFrame(fEnergyChannel, fMyBarLayout);

  fLabelPosition = new TGLabel(Frame2, "              Position:                                 ");
  fLabelWidth = new TGLabel(Frame2, "           FWHM:                                         ");
  Frame2->AddFrame(fLabelPosition, new TGLayoutHints(kLHintsCenterX, 5, 200, 3, 0));
  Frame2->AddFrame(fLabelWidth, new TGLayoutHints(kLHintsCenterX, 5, 200, 3, 0));

  TGLabel *fTextLabelEff = new TGLabel(Frame3, "Efficiency Calibration.   Insert N_dec:");
  fNumberEff = new TGNumberEntry(Frame3, 0, 10, 999, TGNumberFormat::kNESReal, TGNumberFormat::kNELLimitMinMax, 0, 1e10);
  fNumberEff->Connect("ValueSet(Float_t)", "MyMainFrame", this, "AddPointEff()");
  (fNumberEff->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "AddPointEff()");
  TGLabel *fTextLabelNgamma = new TGLabel(Frame3, "N_g:");
  fnGamma = new TGNumberEntry(Frame3, 0, 5, 999, TGNumberFormat::kNESRealFour, TGNumberFormat::kNELLimitMinMax, 0, 9999);
  fnGamma->Connect("ValueSet(Float_t)", "MyMainFrame", this, "AddPointEff()");
  (fnGamma->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "AddPointEff()");
  fnGamma->SetNumber(0.9999);
  TGTextButton *AddEffPointButton = new TGTextButton(Frame3, "&Add Efficiency point ");
  AddEffPointButton->Connect("Clicked()", "MyMainFrame", this, "AddPointEff()");
  fEff = new TGTextButton(Frame3, "&Consider Efficiency");
  fEff->Connect("Clicked()", "MyMainFrame", this, "ChangeEffLabel()");
  fLabelIntegral = new TGLabel(Frame3, "                       Integral:                                                                ");

  Frame3->AddFrame(fTextLabelEff, new TGLayoutHints(kLHintsLeft, 10, 5, 3, 0));
  Frame3->AddFrame(fNumberEff, new TGLayoutHints(kLHintsLeft, 0, 5, 3, 0));
  Frame3->AddFrame(fTextLabelNgamma, fMyBarLayout);
  Frame3->AddFrame(fnGamma, new TGLayoutHints(kLHintsLeft, 0, 5, 3, 0));
  Frame3->AddFrame(AddEffPointButton, fMyBarLayout);
  Frame3->AddFrame(fEff, fMyBarLayout);
  Frame3->AddFrame(fLabelIntegral, new TGLayoutHints(kLHintsExpandX, 5, 0, 3, 0));

  AddFrame(Frame1, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
  AddFrame(Frame2, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
  AddFrame(Frame3, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
  // Create canvas widget
  fEcanvas = new TRootEmbeddedCanvas("Ecanvas", this, 1500, 1000);
  AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX |
                                           kLHintsExpandY,
                                       10, 10, 10, 1));

  // Create a horizontal frame widget with slider
  TGGroupFrame *hframeRange = new TGGroupFrame(this, "Displayed range");
  TGGroupFrame *hframeFit = new TGGroupFrame(this, "Fitting function");
  fSlider0 = new TGDoubleHSlider(hframeRange);
  fSlider0->Connect("PositionChanged()", "MyMainFrame", this, "ChangeSlider()");
  fSlider0->SetRange(0, 16000);
  fSlider0->SetPosition(0, 16000);
  for (int i = 0; i < 3; i++)
  {
    Slider0Positions[i * 2] = 0;
    Slider0Positions[i * 2 + 1] = 16000;
  }
  fSlider = new TGDoubleHSlider(hframeFit);
  fSlider->Connect("PositionChanged()", "MyMainFrame", this, "ChangeSlider()");
  fSlider->SetRange(0, 1000);
  fSlider->SetPosition(0, 100);
  for (int i = 0; i < 3; i++)
  {
    SliderPositions[i * 2] = 0;
    SliderPositions[i * 2 + 1] = 100;
  }
  hframeRange->AddFrame(fSlider0, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5));
  hframeFit->AddFrame(fSlider, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5));
  AddFrame(hframeRange, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5));
  AddFrame(hframeFit, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5));

  fGraph = new TGraphErrors();
  fGraphEff = new TGraphErrors();
  // Set a name to the main frame
  SetWindowName("Spectrum");

  // Map all subwindows of main frame
  MapSubwindows();

  // Initialize the layout algorithm
  Resize(GetDefaultSize());

  // Map main frame
  MapWindow();
}

void MyMainFrame::ChangeSlider()
{
  ApplyChanges();
  if (NowOpen != ShowEff)
    if (fSlider->GetMinPosition() < fSlider0->GetMinPosition() || fSlider->GetMaxPosition() > fSlider0->GetMaxPosition())
    {
      fSlider->SetPosition(fSlider0->GetMinPosition() + (fSlider0->GetMaxPosition() - fSlider0->GetMinPosition()) / 3, fSlider0->GetMaxPosition() - (fSlider0->GetMaxPosition() - fSlider0->GetMinPosition()) / 3);
    }
  double MinRangeSlider, MaxRangeSlider;

  if (NowOpen == None)
  {
    cout << "Warning: Load the data file!\n";
    return;
  }
  if (NowOpen == ShowSource)
  {
    MinRangeSlider = hLoaded->GetXaxis()->GetXmin() - hLoaded->GetBinCenter(2) + hLoaded->GetBinCenter(1);
    MaxRangeSlider = hLoaded->GetXaxis()->GetXmax() + hLoaded->GetBinCenter(2) - hLoaded->GetBinCenter(1);
    fSlider0->SetRange(MinRangeSlider, MaxRangeSlider);
    fSlider->SetRange(Slider0Positions[0], Slider0Positions[1]);
    Slider0Positions[0] = fSlider0->GetMinPosition();
    Slider0Positions[1] = fSlider0->GetMaxPosition();
    SliderPositions[0] = fSlider->GetMinPosition();
    SliderPositions[1] = fSlider->GetMaxPosition();
    SetShowAs(NowOpen);
  }
  if (NowOpen == ShowBckg)
  {
    MinRangeSlider = hLoaded->GetXaxis()->GetXmin() - hLoaded->GetBinCenter(2) + hLoaded->GetBinCenter(1);
    MaxRangeSlider = hLoaded->GetXaxis()->GetXmax() + hLoaded->GetBinCenter(2) - hLoaded->GetBinCenter(1);
    fSlider0->SetRange(MinRangeSlider, MaxRangeSlider);
    fSlider->SetRange(Slider0Positions[0], Slider0Positions[1]);
    Slider0Positions[0] = fSlider0->GetMinPosition();
    Slider0Positions[1] = fSlider0->GetMaxPosition();
    SliderPositions[0] = fSlider->GetMinPosition();
    SliderPositions[1] = fSlider->GetMaxPosition();
    SetShowAs(NowOpen);
  }
  if (NowOpen == ShowSourceBckg)
  {
    MinRangeSlider = hLoaded->GetXaxis()->GetXmin() - hLoaded->GetBinCenter(2) + hLoaded->GetBinCenter(1);
    MaxRangeSlider = hLoaded->GetXaxis()->GetXmax() + hLoaded->GetBinCenter(2) - hLoaded->GetBinCenter(1);
    fSlider0->SetRange(MinRangeSlider, MaxRangeSlider);
    fSlider->SetRange(Slider0Positions[0], Slider0Positions[1]);
    Slider0Positions[0] = fSlider0->GetMinPosition();
    Slider0Positions[1] = fSlider0->GetMaxPosition();
    SliderPositions[0] = fSlider->GetMinPosition();
    SliderPositions[1] = fSlider->GetMaxPosition();
    SetShowAs(NowOpen);
  }
  if (NowOpen == ShowCalib)
  {
    MinRangeSlider = TMath::MinElement(fGraph->GetN(), fGraph->GetX());
    MaxRangeSlider = TMath::MaxElement(fGraph->GetN(), fGraph->GetX());
    fSlider0->SetRange(MinRangeSlider - (MaxRangeSlider - MinRangeSlider) * 0.01, MaxRangeSlider + (MaxRangeSlider - MinRangeSlider) * 0.01);
    Slider0Positions[2] = fSlider0->GetMinPosition();
    Slider0Positions[3] = fSlider0->GetMaxPosition();
    SliderPositions[2] = fSlider->GetMinPosition();
    SliderPositions[3] = fSlider->GetMaxPosition();
    SetShowAs(NowOpen);
  }
  if (NowOpen == ShowEff)
  {
    EffSliderFirstMoved = kTRUE;
    Slider0Positions[4] = fSlider0->GetMinPosition();
    Slider0Positions[5] = fSlider0->GetMaxPosition();
    SliderPositions[4] = fSlider->GetMinPosition();
    SliderPositions[5] = fSlider->GetMaxPosition();
    NowOpen = None;
    SetShowAs(ShowEff);
  }
  return;
}

void MyMainFrame::LoadSpectrFile(int ForBckg)
{
  ApplyChanges();
  // insert data from file to hist
  double x, istring = 0;
  vector<double> data;
  TGFileInfo fi;
  fMain = new TGTransientFrame();
  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
  if (!fi.fFilename)
    return;
  // printf("Open file: %s (dir: %s)\n", fi.fFilename, fi.fIniDir);
  ifstream inFile;
  inFile.open(fi.fFilename);
  if (!inFile)
  {
    cout << "Unable to open file\n";
    // exit(1); // terminate with error
  }

  while (inFile >> x)
  {
    data.push_back(x);
    istring++;
  }
  inFile.close();
  if (!WasOpened[0])
    WasOpened[0] = 1;
  else
    delete hLoaded;
  hLoaded = new TH1D("hLoaded", ";Channel;", istring / 2, data[0], data[istring - 2]);
  for (int i = 0; i < istring; i += 2)
  {
    hLoaded->Fill(data[i], abs(data[i + 1]));
  }
  if (EnergyNotChannel)
  {
    fEnergyChannel->SetText("&Enable Calibration");
    EnergyNotChannel = kFALSE;
    EnergyNotChannelOneTime = kFALSE;
  }
  NowOpen = ForBckg ? ShowBckg : ShowSource;
  DoDraw();
  if (!ForBckg)
  {
    if (WasOpened[1])
      delete hSignal;
    else
      WasOpened[1] = 1;
    hSignal = (TH1D *)hLoaded->Clone();
  }
  else
  {
    if (WasOpened[2])
      delete hBackground;
    else
      WasOpened[2] = 1;
    hBackground = (TH1D *)hLoaded->Clone();
  }
  return;
}

void MyMainFrame::SetLoadAs(Int_t id)
{
  ApplyChanges();
  switch (id)
  {
  case LoadSpect:
    LoadSpectrFile(0);
    break;
  case LoadBckg:
    LoadSpectrFile(1);
    break;
  case LoadCalib:
    LoadCalibrationFile();
    break;
  case LoadEff:
    LoadEffFile();
    break;
  }
}

void MyMainFrame::SetShowAs(Int_t id)
{
  ApplyChanges();

  if (id == None)
  {
    cout << "Warning: Load file with data!";
    return;
  }
  if (id == ShowSource)
  {
    NowOpen = ShowSource;
    fSlider0->SetPosition(Slider0Positions[0], Slider0Positions[1]);
    fSlider->SetPosition(SliderPositions[0], SliderPositions[1]);
    DoDrawSource();
  }
  if (id == ShowBckg)
  {
    NowOpen = ShowBckg;
    fSlider0->SetPosition(Slider0Positions[0], Slider0Positions[1]);
    fSlider->SetPosition(SliderPositions[0], SliderPositions[1]);
    DoDrawBckg();
  }
  if (id == ShowSourceBckg)
  {
    NowOpen = ShowSourceBckg;
    fSlider0->SetPosition(Slider0Positions[0], Slider0Positions[1]);
    fSlider->SetPosition(SliderPositions[0], SliderPositions[1]);
    DoDrawSourceWithoutBckg();
  }
  if (id == ShowCalib)
  {
    NowOpen = ShowCalib;
    fSlider0->SetPosition(Slider0Positions[2], Slider0Positions[3]);
    fSlider->SetPosition(-99999, 99999);
    DoDrawCalibrationGraph();
  }
  if (id == ShowEff)
  {
    if (NowOpen == ShowEff)
      EffSliderFirstMoved = kFALSE;
    NowOpen = ShowEff;
    fSlider0->SetPosition(Slider0Positions[4], Slider0Positions[5]);
    fSlider->SetPosition(SliderPositions[4], SliderPositions[5]);
    DoDrawEff();
  }
  return;
}

void MyMainFrame::SetSaveAs(Int_t id)
{
  ApplyChanges();

  if (id == SaveImage)
    Save();
  if (id == SaveCalib)
    SaveCal();
  if (id == SaveEffic)
    SaveEff();
  return;
}

void MyMainFrame::DoDrawSource()
{
  if (!WasOpened[1])
  {
    cout << "WARNING: No source file selected!\n";
    LoadSpectrFile(0);
    hSignal = (TH1D *)hLoaded->Clone();
    WasOpened[1] = 1;
  }

  hLoaded = (TH1D *)hSignal->Clone();
  if (EnergyNotChannel)
    EnergyNotChannelOneTime = kTRUE;
  if (HistBinCombined > 1)
    Rebinned = kTRUE;
  DoDraw();
}
void MyMainFrame::DoDrawBckg()
{
  if (!WasOpened[2])
  {
    cout << "WARNING: No background file selected!\n";
    LoadSpectrFile(1);
    hBackground = (TH1D *)hLoaded->Clone();
    WasOpened[2] = 1;
  }

  hLoaded = (TH1D *)hBackground->Clone();
  if (EnergyNotChannel)
    EnergyNotChannelOneTime = kTRUE;
  if (HistBinCombined > 1)
    Rebinned = kTRUE;
  DoDraw();
}
void MyMainFrame::DoDrawSourceWithoutBckg()
{
  if (!WasOpened[1])
  {
    cout << "WARNING: No source file selected!\n";
    LoadSpectrFile(0);
    return;
  }
  if (!WasOpened[2])
  {
    cout << "WARNING: No background file selected!\n";
    LoadSpectrFile(1);
    return;
  }
  hLoaded = (TH1D *)hBackground->Clone();
  hLoaded->Scale(-1.0);
  hLoaded->Add(hSignal);
  for (int ibin = 0; ibin < hLoaded->GetNbinsX(); ibin++)
    if (hLoaded->GetBinContent(ibin + 1) < 0)
      hLoaded->SetBinContent(ibin + 1, 0);
  if (EnergyNotChannel)
    EnergyNotChannelOneTime = kTRUE;
  if (HistBinCombined > 1)
    Rebinned = kTRUE;
  DoDraw();
}
void MyMainFrame::Rebin()
{
  HistBinCombined = fHistBinCombined->GetNumberEntry()->GetNumber();
  Rebinned = kTRUE;
  DoDraw();
}
void MyMainFrame::DoDraw()
{
  ApplyChanges();
  if (!WasOpened[0])
    return;
  if (EnergyNotChannelOneTime)
  {
    TempHist = new TH1D("TempHist", "", hLoaded->GetNbinsX(), hLoaded->GetXaxis()->GetXmin() * a + b, hLoaded->GetXaxis()->GetXmax() * a + b);
    for (int ibin = 0; ibin < hLoaded->GetNbinsX(); ibin++)
      TempHist->SetBinContent(ibin, hLoaded->GetBinContent(ibin));
    delete hLoaded;
    hLoaded = (TH1D *)TempHist->Clone();
    delete TempHist;
    EnergyNotChannelOneTime = kFALSE;
  }

  if (Rebinned)
  {
    int NbinsNew = hLoaded->GetNbinsX() / HistBinCombined;
    TempHist = new TH1D("TempHist", "", NbinsNew, hLoaded->GetXaxis()->GetXmin(), hLoaded->GetXaxis()->GetXmax());
    for (int ibin = 0; ibin < hLoaded->GetNbinsX(); ibin++)
      TempHist->Fill(hLoaded->GetBinCenter(ibin), hLoaded->GetBinContent(ibin));
    TempHist->Scale(1.0 / HistBinCombined);
    delete hLoaded;
    hLoaded = (TH1D *)TempHist->Clone();
    delete TempHist;
    Rebinned = kFALSE;
  }
  fSlider0->SetRange(hLoaded->GetXaxis()->GetXmin() - hLoaded->GetBinCenter(2) + hLoaded->GetBinCenter(1), hLoaded->GetXaxis()->GetXmax() + hLoaded->GetBinCenter(2) - hLoaded->GetBinCenter(1));
  fSlider->SetRange(Slider0Positions[0], Slider0Positions[1]);
  Double_t MinPos = SliderPositions[0], MaxPos = SliderPositions[1];
  Double_t borderLength = fBorderLength->GetNumberEntry()->GetNumber();
  if (EnergyNotChannel)
    borderLength *= a;
  TF1 *fitLineUnderPeak = new TF1("fitLineUnderPeak", "pol1", MinPos - borderLength, MaxPos + borderLength);
  fitLineUnderPeak->SetLineColor(kGreen);
  hLoaded->SetStats(0);
  hLoaded->SetTitle(EnergyNotChannel ? ";Energy, keV;" : ";Channel;");
  hLoaded->GetXaxis()->SetRangeUser(Slider0Positions[0], Slider0Positions[1]);
  TGraph *helpGraph = new TGraph();
  int Binwidth = hLoaded->GetBinCenter(2) - hLoaded->GetBinCenter(1);
  for (int i = 0; i < fBorderLength->GetNumberEntry()->GetNumber(); ++i)
    helpGraph->SetPoint(helpGraph->GetN(), hLoaded->GetBinCenter(hLoaded->GetXaxis()->FindBin(MinPos) - i), hLoaded->GetBinContent(hLoaded->GetXaxis()->FindBin(MinPos) - i)); //-i
  for (int i = 0; i < fBorderLength->GetNumberEntry()->GetNumber(); ++i)
    helpGraph->SetPoint(helpGraph->GetN(), hLoaded->GetBinCenter(hLoaded->GetXaxis()->FindBin(MaxPos) + i), hLoaded->GetBinContent(hLoaded->GetXaxis()->FindBin(MaxPos) + i)); //+i
  helpGraph->Fit("fitLineUnderPeak", "RQ");
  hLoaded->Add(fitLineUnderPeak, -1);
  NumberOfPeacks = fNumberOfPeacks->GetNumberEntry()->GetNumber();
  TF1 *fitGaus = new TF1("fitGaus", "gaus", MinPos, MaxPos);
  TF1 *fit[10];
  TString FuncName = "";
  for (int i = 0; i < NumberOfPeacks; i++)
  {
    fit[i] = new TF1(Form("fit%d", i), "[0]*exp(-0.5*((x-[1])/[2])**2)", MinPos, MaxPos);
    if (i > 0)
      FuncName += "+";
    FuncName += Form("fit%d", i);
  }
  TF1 *SumOfGaus = new TF1("SumOfGaus", FuncName, MinPos, MaxPos);
  for (int i = 0; i < NumberOfPeacks; i++)
  {
    SumOfGaus->SetParameter(i * 3, 100);
    SumOfGaus->SetParameter(1 + i * 3, (MinPos + MaxPos) / 2);
    SumOfGaus->SetParameter(2 + i * 3, 3);
  }

  SumOfGaus->SetLineColor(kYellow);
  if (NumberOfPeacks == 1)
    hLoaded->Fit("fitGaus", "RQ");
  else
    hLoaded->Fit("SumOfGaus", "RQ");
  hLoaded->Add(fitLineUnderPeak, 1);
  // hLoaded->SetMinimum(1);
  if (Logy)
    gPad->SetLogy();
  else
    gPad->SetLogy(kFALSE);
  hLoaded->Draw("hist");
  // helpGraph->Draw("same");
  // SumOfGaus->Draw("same");
  fitLineUnderPeak->Draw("same");
  // visual-only
  // Double_t a1 = fitLineUnderPeak->GetParameter(1), a2 = fitLineUnderPeak->GetParameter(0);
  // TF1 *visualfit = new TF1("visualfit", Form("fit + %f*x+%f", a1, a2), MinPos, MaxPos);
  TF1 *visualfit;

  PeakPosition.clear();
  PeakPositionUncertancy.clear();
  PeakIntegral.clear();
  PeakIntegralUncertancy.clear();
  nPeaksToAddCal = NumberOfPeacks;
  nPeaksToAddEff = NumberOfPeacks;
  if (NumberOfPeacks == 1)
  {
    visualfit = new TF1("visualfit", "fitLineUnderPeak+fitGaus", MinPos, MaxPos);

    PeakPosition.push_back(fitGaus->GetParameter(1));
    PeakPositionUncertancy.push_back(fitGaus->GetParameter(2) * 2.35);
    PeakIntegral.push_back(fitGaus->Integral(MinPos, MaxPos) / (EnergyNotChannel ? a : 1.));
    PeakIntegralUncertancy.push_back(TMath::Sqrt(PeakIntegral[0]));
    PeakIntegral[0] /= fMeasureTime->GetNumberEntry()->GetNumber();
    PeakIntegralUncertancy[0] /= fMeasureTime->GetNumberEntry()->GetNumber();
    if (Eff)
    {
      if (EffSliderFirstMoved)
      {
        PeakIntegral[0] /= FitEff->Eval(PeakPosition[0]);
        PeakIntegralUncertancy[0] = sqrt(pow(PeakIntegralUncertancy[0] / FitEff->Eval(PeakPosition[0]), 2) + pow(PeakIntegral[0] / FitEff->Eval(PeakPosition[0]) * FitEffError, 2));
      }
      else
      {
        if (EnergyNotChannel)
          PeakPosition[0] = (PeakPosition[0] - b) / a;
        cout << "fGraphEff: " << fGraphEff->Eval(PeakPosition[0]) << endl;
        PeakIntegralUncertancy[0] /= fGraphEff->Eval(PeakPosition[0]);
        PeakIntegral[0] /= fGraphEff->Eval(PeakPosition[0]);
        if (EnergyNotChannel)
          PeakPosition[0] = a * PeakPosition[0] + b;
      }
      PeakIntegral[0] /= fnGamma->GetNumberEntry()->GetNumber();
      PeakIntegralUncertancy[0] /= fnGamma->GetNumberEntry()->GetNumber();
    }
    fLabelPosition->SetText(Form("          Position: %.2f                    ", PeakPosition[0]));
    fLabelWidth->SetText(Form("     FWHM: %.2f          ", PeakPositionUncertancy[0]));
    fLabelIntegral->SetText(Form("          Integral: %.2f +- %.2f      ", PeakIntegral[0], PeakIntegralUncertancy[0]));
  }
  else
  {
    // FuncName = "+fitLineUnderPeak";
    visualfit = new TF1("visualfit", "fitLineUnderPeak+SumOfGaus", MinPos, MaxPos);
    visualfit->SetParameter(0, fitLineUnderPeak->GetParameter(0));
    visualfit->SetParameter(1, fitLineUnderPeak->GetParameter(1));
    for (int i = 0; i < 3 * NumberOfPeacks; i++)
    {
      visualfit->SetParameter(2 + i, SumOfGaus->GetParameter(i));
    }
    TF1 *SinglePeak[10], *visualSinglePeak[10];
    TString PositionsText = "Positions:";
    TString IntegralsText = "Integrals:";
    TString FWHMsText = "FWHMs:";
    vector<double> MultiplePeaksParams[3];
    for (int i = 0; i < NumberOfPeacks; i++)
    {
      SinglePeak[i] = new TF1(Form("SinglePeak%d", i), "fitGaus", MinPos, MaxPos);
      // SinglePeak[i]->SetParameter(0, fitLineUnderPeak->GetParameter(0));
      // SinglePeak[i]->SetParameter(1, fitLineUnderPeak->GetParameter(1));
      SinglePeak[i]->SetParameter(0, SumOfGaus->GetParameter(i * 3));
      SinglePeak[i]->SetParameter(1, SumOfGaus->GetParameter(1 + i * 3));
      SinglePeak[i]->SetParameter(2, SumOfGaus->GetParameter(2 + i * 3));
      visualSinglePeak[i] = new TF1("visualSinglePeak", "fitLineUnderPeak+fitGaus", MinPos, MaxPos);
      visualSinglePeak[i]->SetParameter(0, fitLineUnderPeak->GetParameter(0));
      visualSinglePeak[i]->SetParameter(1, fitLineUnderPeak->GetParameter(1));
      visualSinglePeak[i]->SetParameter(2, SumOfGaus->GetParameter(i * 3));
      visualSinglePeak[i]->SetParameter(3, SumOfGaus->GetParameter(1 + i * 3));
      visualSinglePeak[i]->SetParameter(4, SumOfGaus->GetParameter(2 + i * 3));
      visualSinglePeak[i]->SetFillColor(color[i]);
      visualSinglePeak[i]->SetFillStyle(3345);
      visualSinglePeak[i]->SetLineColor(color[i]);
      visualSinglePeak[i]->DrawCopy("same");
      MultiplePeaksParams[0].push_back(SinglePeak[i]->GetParameter(1));
      MultiplePeaksParams[1].push_back(SinglePeak[i]->Integral(MinPos, MaxPos) / (EnergyNotChannel ? a : 1.));
      MultiplePeaksParams[2].push_back(SinglePeak[i]->GetParameter(2) * 2.35);
    }

    for (int i = 0; i < NumberOfPeacks; i++)
    {
      int ind = min_element(MultiplePeaksParams[0].begin(), MultiplePeaksParams[0].end()) - MultiplePeaksParams[0].begin();
      PeakPosition.push_back(MultiplePeaksParams[0][ind]);
      PeakIntegral.push_back(MultiplePeaksParams[1][ind]);
      PeakPositionUncertancy.push_back(MultiplePeaksParams[2][ind]);
      PeakIntegralUncertancy.push_back(TMath::Sqrt(PeakIntegral.back()));
      for (int j = 0; j < 3; j++)
        MultiplePeaksParams[j].erase(MultiplePeaksParams[j].begin() + ind);

      PeakIntegral.back() /= fMeasureTime->GetNumberEntry()->GetNumber(); 
      PeakIntegralUncertancy.back() /= fMeasureTime->GetNumberEntry()->GetNumber(); 
      if (Eff)
      {
        if (EffSliderFirstMoved)
        {
          PeakIntegral.back() /= FitEff->Eval(PeakPosition.back());
          PeakIntegralUncertancy.back() = sqrt(pow(PeakIntegralUncertancy.back() / FitEff->Eval(PeakPosition.back()), 2) + pow(PeakIntegral.back() / FitEff->Eval(PeakPosition.back()) * FitEffError, 2));
        }
        else
        {
          if (EnergyNotChannel)
            PeakPosition.back() = (PeakPosition.back() - b) / a;
          PeakIntegralUncertancy.back() /= fGraphEff->Eval(PeakPosition.back());
          PeakIntegral.back() /= fGraphEff->Eval(PeakPosition.back());
          if (EnergyNotChannel)
            PeakPosition.back() = a * PeakPosition.back() + b;
        }
        PeakIntegralUncertancy.back() /= fnGamma->GetNumberEntry()->GetNumber();
        PeakIntegral.back() /= fnGamma->GetNumberEntry()->GetNumber();
      }
      PositionsText += Form("  %.2f;", PeakPosition.back());
      IntegralsText += Form("  %.2f +- %.2f;", PeakIntegral.back(), PeakIntegralUncertancy.back());
      FWHMsText += Form("  %.2f;", PeakPositionUncertancy.back());
    }

    fLabelPosition->SetText(PositionsText);
    fLabelWidth->SetText(FWHMsText);
    fLabelIntegral->SetText(IntegralsText);
  }
  visualfit->SetLineColor(kRed);
  visualfit->Draw("same");

  //  Parent frame Layout() method will redraw the label showing the new value.
  Frame1->Layout();

  fEcanvas->GetCanvas()->Modified();
  fEcanvas->GetCanvas()->Update();
}

void MyMainFrame::Save()
{
  const char *filetypes[] = {"*.png", "*.png", "All files", "*", 0, 0};
  TGFileInfo fi;
  fi.fFileTypes = filetypes;
  new TGFileDialog(gClient->GetRoot(), fMain, kFDSave, &fi);
  TCanvas *c2 = fEcanvas->GetCanvas();
  c2->SaveAs(fi.fFilename, "RECREATE");
}

void MyMainFrame::SaveCal()
{
  // multimap<double, double> sortedECal;
  double x, y, yerr;
  vector<double> xArr, yArr, yerrArr;
  for (int i = 0; i < fGraph->GetN(); i++)
  {
    fGraph->GetPoint(i, x, y);
    yerr = fGraph->GetErrorY(i);
    // sortedECal.insert(make_pair(x, i));
    xArr.push_back(x);
    yArr.push_back(y);
    yerrArr.push_back(yerr);
  }
  const char *filetypes[] = {"*.txt", "*.txt", "All files", "*", 0, 0};
  TGFileInfo fi;
  fi.fFileTypes = filetypes;
  new TGFileDialog(gClient->GetRoot(), fMain, kFDSave, &fi);
  ofstream outfile(fi.fFilename);
  // ofstream outfile2("CalibrationWithFWHM.txt");
  for (int i = 0; i < xArr.size(); i++)
  {
    int ind = min_element(xArr.begin(), xArr.end()) - xArr.begin();
    outfile << Form("%.2f", xArr[ind]) << " " << Form("%.2f", yArr[ind]) << " " << Form("%.2f", yerrArr[ind]) << endl;
    xArr.erase(xArr.begin() + ind);
    yArr.erase(yArr.begin() + ind);
    yerrArr.erase(yerrArr.begin() + ind);
  }
  outfile.close();
}

void MyMainFrame::SaveEff()
{
  const char *filetypes[] = {"*.txt", "*.txt", "All files", "*", 0, 0};
  TGFileInfo fi;
  fi.fFileTypes = filetypes;
  new TGFileDialog(gClient->GetRoot(), fMain, kFDSave, &fi);
  ofstream outfile(fi.fFilename);

  double x, y, xerr, yerr;
  int nPoints = fGraphEff->GetN();
  for (int i = 0; i < nPoints; i++)
  {
    fGraphEff->GetPoint(i, x, y);
    yerr = fGraphEff->GetErrorY(i);
    if (EnergyNotChannel)
      outfile << x * a + b << " " << y << " " << yerr << endl;
    else
      outfile << x << " " << y << " " << yerr << endl;
  }
  outfile.close();
}
//__________________________________________________________________________________________________________________________________________

void MyMainFrame::LoadCalibrationFile()
{
  ApplyChanges();
  // insert data from file to hist
  double nPoints = 0, x, y, yerr;
  vector<double> dataX;
  vector<double> dataY;
  vector<double> dataYerr;
  TGFileInfo fi;
  fMain = new TGTransientFrame();
  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
  if (!fi.fFilename)
    return;
  // printf("Open file: %s (dir: %s)\n", fi.fFilename, fi.fIniDir);
  ifstream inFile;
  inFile.open(fi.fFilename);
  if (!inFile)
  {
    cout << "Unable to open file\n";
    // exit(1); // terminate with error
  }

  while (inFile >> x >> y >> yerr)
  {
    dataX.push_back(x);
    dataY.push_back(y);
    dataYerr.push_back(yerr);
    nPoints++;
  }
  inFile.close();
  fGraph->Set(0);
  for (int i = 0; i < nPoints; i++)
    fGraph->SetPoint(i, dataX[i], dataY[i]);

  TF1 *fitline = new TF1("fitline", "pol1", 0, TMath::MaxElement(fGraph->GetN(), fGraph->GetX()) * 1.1);
  fGraph->Fit(fitline, "QW");
  a = fitline->GetParameter(1);
  for (int i = 0; i < nPoints; i++)
    fGraph->SetPointError(i, dataYerr[i] / a, dataYerr[i]);

  SetShowAs(ShowCalib);
}

void MyMainFrame::LoadEffFile()
{
  ApplyChanges();
  // insert data from file to hist
  double nPoints = 0, x, y, yerr;
  vector<double> dataX, dataY, dataYerr;
  TGFileInfo fi;
  fMain = new TGTransientFrame();
  new TGFileDialog(gClient->GetRoot(), fMain, kFDOpen, &fi);
  if (!fi.fFilename)
    return;
  // printf("Open file: %s (dir: %s)\n", fi.fFilename, fi.fIniDir);
  ifstream inFile;
  inFile.open(fi.fFilename);
  if (!inFile)
  {
    cout << "Unable to open file\n";
    // exit(1); // terminate with error
  }

  while (inFile >> x >> y >> yerr)
  {
    if (EnergyNotChannel)
      x = (x - b) / a;
    dataX.push_back(x);
    dataY.push_back(y);
    dataYerr.push_back(yerr);
    nPoints++;
  }
  inFile.close();
  fGraphEff->Set(0);
  for (int i = 0; i < nPoints; i++)
  {
    fGraphEff->SetPoint(i, dataX[i], dataY[i]);
    fGraphEff->SetPointError(i, 0, dataYerr[i]);
  }
  EffSliderFirstMoved = kFALSE;
  SetShowAs(ShowEff);
}

void MyMainFrame::AddPoint()
{
  ApplyChanges();
  if (nPeaksToAddCal <= 0)
  {
    cout << "No peaks to add!\n";
    return;
  }
  fGraph->SetPoint(fGraph->GetN(), PeakPosition[PeakPosition.size() - nPeaksToAddCal], fNumber->GetNumberEntry()->GetNumber());
  fGraph->SetPointError(fGraph->GetN() - 1, PeakPositionUncertancy[PeakPosition.size() - nPeaksToAddCal], 0);
  nPeaksToAddCal--;
}

void MyMainFrame::AddPointEff()
{
  ApplyChanges();
  if (Eff) // To divide integral by Ng
  {
    DoDraw();
    return;
  }
  if (nPeaksToAddEff <= 0)
  {
    cout << "No peaks to add!\n";
    return;
  }
  if (EnergyNotChannel)
    fGraphEff->SetPoint(fGraphEff->GetN(), (PeakPosition[PeakPosition.size() - nPeaksToAddEff] - b) / a, PeakIntegral[PeakIntegral.size() - nPeaksToAddEff] / (fNumberEff->GetNumberEntry()->GetNumber() * fnGamma->GetNumberEntry()->GetNumber()));
  else
    fGraphEff->SetPoint(fGraphEff->GetN(), PeakPosition[PeakPosition.size() - nPeaksToAddEff], PeakIntegral[PeakIntegral.size() - nPeaksToAddEff] / (fNumberEff->GetNumberEntry()->GetNumber() * fnGamma->GetNumberEntry()->GetNumber()));
  fGraphEff->SetPointError(fGraphEff->GetN() - 1, 0, sqrt(PeakIntegral[PeakIntegral.size() - nPeaksToAddEff]) / (fNumberEff->GetNumberEntry()->GetNumber() * fnGamma->GetNumberEntry()->GetNumber()));
  nPeaksToAddEff--;
}

void MyMainFrame::DoDrawCalibrationGraph()
{
  ApplyChanges();
  if (fGraph->GetN() < 2)
  {
    cout << "Warning: Calibration Graph have not enough points!";
    return;
  }
  TF1 *fitline = new TF1("fitline", "pol1", 0, TMath::MaxElement(fGraph->GetN(), fGraph->GetX()) * 1.1);
  fGraph->Fit(fitline, "QW");
  fGraph->SetMarkerStyle(21);
  fitline->SetTitle(";Channel;Energy ,keV");
  a = fitline->GetParameter(1);
  b = fitline->GetParameter(0);

  TGraphErrors *CalibrationGraphDev = new TGraphErrors();
  double x, y, xerr;
  for (int i = 0; i < fGraph->GetN(); i++)
  {
    fGraph->GetPoint(i, x, y);
    xerr = fGraph->GetErrorX(i);
    fGraph->SetPointError(i, xerr, xerr * a);
    CalibrationGraphDev->SetPoint(i, x, y - a * x - b);
    CalibrationGraphDev->SetPointError(i, xerr, xerr * a);
  }
  fEcanvas->GetCanvas()->Clear();
  TPad *pad[2];
  pad[0] = new TPad("pad0", "", 0, 0.4, 1, 1.0);
  pad[0]->SetBottomMargin(.06); // Upper and lower plot are joined
  pad[0]->Draw();
  pad[1] = new TPad("pad1", "", 0, 0.0, 1, 0.31);
  pad[1]->SetTopMargin(0);      // Upper and lower plot are joined
  pad[1]->SetBottomMargin(0.3); // Upper and lower plot are joined
  pad[1]->Draw();

  double MinRangeSlider = TMath::MinElement(fGraph->GetN(), fGraph->GetX());
  double MaxRangeSlider = TMath::MaxElement(fGraph->GetN(), fGraph->GetX());
  fSlider0->SetRange(MinRangeSlider - (MaxRangeSlider - MinRangeSlider) * 0.01, MaxRangeSlider + (MaxRangeSlider - MinRangeSlider) * 0.01);
  fitline->GetXaxis()->SetRangeUser(Slider0Positions[2], Slider0Positions[3]);
  pad[0]->cd();
  fitline->GetYaxis()->SetTitleSize(0.06);
  fitline->GetYaxis()->SetTitleOffset(0.5);
  fitline->GetYaxis()->SetLabelSize(0.06);
  fitline->SetTitle(";;E, keV");
  fitline->GetXaxis()->SetTitleSize(0.04);
  fitline->GetXaxis()->SetLabelSize(0.06);
  fitline->Draw();
  fGraph->Draw("p,same");
  pad[1]->cd();
  TF1 *ZeroLine = new TF1("ZeroLine", "0", 0, TMath::MaxElement(fGraph->GetN(), fGraph->GetX()) * 1.1);
  ZeroLine->SetLineColor(kBlack);
  ZeroLine->SetMinimum(TMath::MinElement(CalibrationGraphDev->GetN(), CalibrationGraphDev->GetY()) * 1.15);
  ZeroLine->SetMaximum(TMath::MaxElement(CalibrationGraphDev->GetN(), CalibrationGraphDev->GetY()) * 1.15);
  ZeroLine->SetLineStyle(2);
  // ZeroLine->GetXaxis()->SetTitle("Channel");
  ZeroLine->SetTitle(";Channel;#Delta E");
  ZeroLine->GetXaxis()->SetTitleSize(0.13);
  ZeroLine->GetXaxis()->SetLabelSize(0.12);
  ZeroLine->GetYaxis()->SetTitleSize(0.1);
  ZeroLine->GetYaxis()->SetLabelSize(0.1);
  ZeroLine->GetYaxis()->SetTitleOffset(0.3);
  ZeroLine->GetXaxis()->SetRangeUser(Slider0Positions[2], Slider0Positions[3]);
  ZeroLine->Draw("");
  CalibrationGraphDev->SetMarkerStyle(20);
  CalibrationGraphDev->Draw("p,same");

  fEcanvas->GetCanvas()->cd();
  fEcanvas->GetCanvas()->Modified();
  fEcanvas->GetCanvas()->Update();
  WasOpened[3] = 1;
}

void MyMainFrame::DoDrawEff()
{
  ApplyChanges();
  int nPoints = fGraphEff->GetN();
  double x, y, xerr, yerr;
  vector<double> xArr, yArr, xerrArr, yerrArr;
  if (EnergyNotChannel)
  {
    for (int i = 0; i < nPoints; i++)
    {
      fGraphEff->GetPoint(i, x, y);
      yerr = fGraphEff->GetErrorY(i);
      xArr.push_back(x);
      yArr.push_back(y);
      yerrArr.push_back(yerr);
    }
    fGraphEff->Set(0);
    for (int i = 0; i < nPoints; i++)
    {
      int ind = min_element(xArr.begin(), xArr.end()) - xArr.begin();
      fGraphEff->SetPoint(i, xArr[ind] * a + b, yArr[ind]);
      fGraphEff->SetPointError(i, 0, yerrArr[ind]);
      xArr.erase(xArr.begin() + ind);
      yArr.erase(yArr.begin() + ind);
      yerrArr.erase(yerrArr.begin() + ind);
    }
    GraphEffChanged = kTRUE;
  }
  fGraphEff->SetMarkerStyle(23);
  double GraphXmin = TMath::MinElement(fGraphEff->GetN(), fGraphEff->GetX());
  double GraphXmax = TMath::MaxElement(fGraphEff->GetN(), fGraphEff->GetX());
  double GraphYmax = TMath::MaxElement(fGraphEff->GetN(), fGraphEff->GetY());
  // double GraphYmaxInd = TMath::LocMax(fGraphEff->GetN(), fGraphEff->GetY());
  // double GraphYmax = fGraphEff->GetPointY(GraphYmaxInd);

  fSlider0->SetRange(GraphXmin - (GraphXmax - GraphXmin) * 0.01, GraphXmax + (GraphXmax - GraphXmin) * 0.01);
  fSlider->SetRange(GraphXmin - (GraphXmax - GraphXmin) * 0.01, GraphXmax + (GraphXmax - GraphXmin) * 0.01);
  // fSlider->SetRange(Slider0Positions[4], Slider0Positions[5]);
  double BreakingFitPoint1 = Slider0Positions[5];
  double BreakingFitPoint2 = SliderPositions[4];

  TF1 *FitEff1 = new TF1("FitEff1", "[0]*pow(x-[1],2)+[2]", GraphXmin - 0.01, BreakingFitPoint1 * 2);
  TF1 *FitEff2 = new TF1("FitEff2", "[0]/(x+[1])+[2]", BreakingFitPoint2, GraphXmax + 0.01);
  // FitEff1->SetParLimits(1, 0 , BreakingFitPoint1);
  // FitEff2->SetParLimits(1, BreakingFitPoint2 , GraphXmax);
  // FitEff2->SetParLimits(2, 0, fGraphEff->GetPointY(TMath::LocMax(fGraphEff->GetN(), fGraphEff->GetX())));
  TGraphErrors *helpGraphEff = new TGraphErrors();
  for (int i = 0; i < nPoints; i++)
  {
    fGraphEff->GetPoint(i, x, y);
    yerr = fGraphEff->GetErrorY(i);
    if (x < BreakingFitPoint1)
    {
      helpGraphEff->SetPoint(helpGraphEff->GetN(), x, y);
      helpGraphEff->SetPointError(helpGraphEff->GetN() - 1, 0, yerr);
    }
    else
      break;
  }
  for (int i = 0; i < nPoints; i++)
  {
    fGraphEff->GetPoint(i, x, y);
    yerr = fGraphEff->GetErrorY(i);
    if (x < BreakingFitPoint1)
    {
      helpGraphEff->SetPoint(helpGraphEff->GetN(), 2 * BreakingFitPoint1 - x, y);
      helpGraphEff->SetPointError(helpGraphEff->GetN() - 1, 0, yerr);
    }
    else
      break;
  }
  helpGraphEff->Fit("FitEff1", "0QR");
  fGraphEff->Fit("FitEff2", "0QR");

  if (FitEff1->GetParameter(2) != FitEff2->GetParameter(2))
    BreakingFitPoint2 = FitEff2->GetParameter(0) / (FitEff1->GetParameter(2) - FitEff2->GetParameter(2)) - FitEff2->GetParameter(1);
  if (BreakingFitPoint2 - BreakingFitPoint1 == 0)
    BreakingFitPoint2 -= 0.01;

  BreakingFitPoint1 = Slider0Positions[4];
  BreakingFitPoint2 = SliderPositions[5];
  double middleA = (FitEff2->Eval(BreakingFitPoint2) - FitEff1->Eval(BreakingFitPoint1)) / (BreakingFitPoint2 - BreakingFitPoint1);
  double middleB = FitEff2->Eval(BreakingFitPoint2) - middleA * BreakingFitPoint2;

  FitEff = new TF1("FitEff", Form("(x<%f)*(FitEff1)+(x>=%f&&x<=%f)*(%f*x+%f)+(x>%f)*(FitEff2)", BreakingFitPoint1, BreakingFitPoint1, BreakingFitPoint2, middleA, middleB, BreakingFitPoint2),
                   GraphXmin - (GraphXmax - GraphXmin) * 0.01, GraphXmax + (GraphXmax - GraphXmin) * 0.01);
  FitEff->SetParameter(0, FitEff1->GetParameter(0));
  FitEff->SetParameter(1, FitEff1->GetParameter(1));
  FitEff->SetParameter(2, FitEff1->GetParameter(2));
  FitEff->SetParameter(3, FitEff2->GetParameter(0));
  FitEff->SetParameter(4, FitEff2->GetParameter(1));
  FitEff->SetParameter(5, FitEff2->GetParameter(2));
  FitEff->SetLineColor(kRed);
  FitEff->SetNpx(1000);
  FitEff->SetTitle(EnergyNotChannel ? ";Energy ,keV;Efficiency" : ";Channel;Efficiency");
  fGraphEff->SetTitle(EnergyNotChannel ? ";Energy ,keV;Efficiency" : ";Channel;Efficiency");
  fGraphEff->SetLineColor(kRed);
  // fGraphEff->Fit("FitEff", "0Q");
  // FitEff->GetXaxis()->SetRangeUser(Slider0Positions[4], Slider0Positions[5]);
  FitEff->GetYaxis()->SetRangeUser(0, GraphYmax * 1.1);
  FitEff->SetMinimum(.00001);
  if (!EffSliderFirstMoved)
    FitEff->SetLineColor(kWhite);
  FitEff->Draw("");
  fGraphEff->Draw(EffSliderFirstMoved ? "pe, same" : "lpe");
  // TF1 *VerticalLine1 = new TF1("VerticalLine1", Form("(x-%f)*5", BreakingFitPoint1), GraphXmin - (GraphXmax - GraphXmin) * 0.01, GraphXmax + (GraphXmax - GraphXmin) * 0.01);
  // TF1 *VerticalLine2 = new TF1("VerticalLine2", Form("(x-%f)*5", BreakingFitPoint2), GraphXmin - (GraphXmax - GraphXmin) * 0.01, GraphXmax + (GraphXmax - GraphXmin) * 0.01);
  // VerticalLine1->SetLineColor(kBlack);
  // VerticalLine2->SetLineColor(kBlack);
  // VerticalLine1->SetLineStyle(8);
  // VerticalLine2->SetLineStyle(8);
  // VerticalLine1->SetNpx(100000);
  // VerticalLine2->SetNpx(100000);
  // VerticalLine1->Draw("same");
  // VerticalLine2->Draw("same");

  // helpGraphEff->SetMarkerStyle(22);
  // helpGraphEff->SetMarkerColor(kGreen);
  // helpGraphEff->Draw("p, same");
  // FitEff3->SetLineColor(kGreen);
  // FitEff3->Draw("same");
  FitEffError = 0;
  for (int i = 0; i < fGraphEff->GetN(); i++)
  {
    fGraphEff->GetPoint(i, x, y);
    FitEffError += pow(y - FitEff->Eval(x), 2);
  }
  FitEffError = sqrt(FitEffError);
  FitEffError /= fGraphEff->GetN();
  TGraph *GraphShade = new TGraph(100 * 2);
  for (int i = 0; i < 100; i++)
  {
    x = GraphXmin - (GraphXmax - GraphXmin) * 0.01 + (GraphXmax - GraphXmin) * 1.01 / 100 * i;
    y = FitEff->Eval(x) + FitEffError;
    GraphShade->SetPoint(i, x, y);
    x = GraphXmin - (GraphXmax - GraphXmin) * 0.01 + (GraphXmax - GraphXmin) * 1.01 / 100 * (100 - i);
    y = FitEff->Eval(x) - FitEffError;
    GraphShade->SetPoint(i + 100, x, y);
  }
  GraphShade->SetFillStyle(3244);
  GraphShade->SetFillColor(kBlue);
  if (EffSliderFirstMoved)
    GraphShade->Draw("f, same");
  fEcanvas->GetCanvas()->Modified();
  fEcanvas->GetCanvas()->Update();
}

MyMainFrame::~MyMainFrame()
{
  // Clean up used widgets: frames, buttons, layout hints
  Cleanup();
}
void spectrum()
{
  // Popup the GUI...
  new MyMainFrame(gClient->GetRoot(), 20, 20);
}
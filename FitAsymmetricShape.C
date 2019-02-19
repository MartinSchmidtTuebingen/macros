#include <iostream>
#include "TMath.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TTree.h"

#include "TRandom3.h"
#include "THnSparse.h"
#include "TStopwatch.h"

#include "THnSparseDefinitions.h"

const Double_t OneOverSqrt2 = 1./TMath::Sqrt2();


/*
//_____________________________________________________________________________
// This function is more or less a copy of AliAnalysisTaskPID::SetConvolutedGaussLambdaParameter
Bool_t SetConvolutedGaussLambdaParameter(Double_t lambda)
{
  // Set the lambda parameter of the convolution to the desired value. Automatically
  // calculates the parameters for the transition (restricted) gauss -> convoluted gauss.
  fConvolutedGaussTransitionPars[2] = lambda;
  
  // Save old parameters and settings of function to restore them later:
  Double_t* oldFuncParams = new Double_t[fkConvolutedGausNPar];
  for (Int_t i = 0; i < fkConvolutedGausNPar; i++)
    oldFuncParams[i] = fConvolutedGausDeltaPrime->GetParameter(i);
  Int_t oldFuncNpx = fConvolutedGausDeltaPrime->GetNpx();
  Double_t oldFuncRangeLow = 0, oldFuncRangeUp = 100;
  fConvolutedGausDeltaPrime->GetRange(oldFuncRangeLow, oldFuncRangeUp);
    
  // Choose some sufficiently large range
  const Double_t rangeStart = 0.5;
  const Double_t rangeEnd = 2.0;
  
  // To get the parameters for the transition, just choose arbitrarily input parameters for mu and sigma
  // (it makes sense to choose typical values). The ratio sigma_gauss / sigma_convolution is independent
  // of mu and as well the difference mu_gauss - mu_convolution.
  Double_t muInput = 1.0;
  Double_t sigmaInput = fgkSigmaReferenceForTransitionPars;
  
  
  // Step 1: Generate distribution with input parameters
  const Int_t nPar = 3;
  Double_t inputPar[nPar] = { muInput, sigmaInput, lambda };

  TH1D* hInput = new TH1D("hInput", "Input distribution", 2000, rangeStart, rangeEnd);

  fConvolutedGausDeltaPrime->SetParameters(inputPar);
  fConvolutedGausDeltaPrime->SetRange(rangeStart, rangeEnd);
  fConvolutedGausDeltaPrime->SetNpx(2000);
  
  TH1* hProbDensity = fConvolutedGausDeltaPrime->GetHistogram();
  
  for (Int_t i = 0; i < 50000000; i++)
    hInput->Fill(hProbDensity->GetRandom());
  
  // Step 2: Fit generated distribution with restricted gaussian
  Int_t maxBin = hInput->GetMaximumBin();
  Double_t max = hInput->GetBinContent(maxBin);
  
  UChar_t usedBins = 1;
  if (maxBin > 1) {
    max += hInput->GetBinContent(maxBin - 1);
    usedBins++;
  }
  if (maxBin < hInput->GetNbinsX()) {
    max += hInput->GetBinContent(maxBin + 1);
    usedBins++;
  }
  max /= usedBins;
  
  // NOTE: The range (<-> fraction of maximum) should be the same
  // as for the spline and eta maps creation
  const Double_t lowThreshold = hInput->GetXaxis()->GetBinLowEdge(hInput->FindFirstBinAbove(0.1 * max));
  const Double_t highThreshold = hInput->GetXaxis()->GetBinUpEdge(hInput->FindLastBinAbove(0.1 * max));
  
  TFitResultPtr fitResGaussFirstStep = hInput->Fit("gaus", "QNRS", "", lowThreshold, highThreshold);
  
  TFitResultPtr fitResGauss;
  
  if ((Int_t)fitResGaussFirstStep == 0) {
    TF1 fGauss("fGauss", "[0]*TMath::Gaus(x, [1], [2], 1)", rangeStart, rangeEnd);
    fGauss.SetParameter(0, fitResGaussFirstStep->GetParams()[0]);
    fGauss.SetParError(0, fitResGaussFirstStep->GetErrors()[0]);
    fGauss.SetParameter(2, fitResGaussFirstStep->GetParams()[2]);
    fGauss.SetParError(2, fitResGaussFirstStep->GetErrors()[2]);

    fGauss.FixParameter(1, fitResGaussFirstStep->GetParams()[1]);
    fitResGauss = hInput->Fit(&fGauss, "QNS", "s", rangeStart, rangeEnd);
  }
  else {
    fitResGauss = hInput->Fit("gaus", "QNRS", "same", rangeStart, rangeEnd);
  }
  
  // Step 3: Use parameters from gaussian fit to obtain parameters for the transition "restricted gauss" -> "convoluted gauss"
  
  // 3.1 The ratio sigmaInput / sigma_gaussFit ONLY depends on lambda (which is fixed per period) -> Calculate this first
  // for an arbitrary (here: typical) sigma. The ratio is then ~the same for ALL sigma for given lambda!
  if ((Int_t)fitResGauss != 0) {
    AliError("Not able to calculate parameters for the transition \"restricted gauss\" -> \"convoluted gauss\": Gauss Fit failed!\n");
    fConvolutedGausDeltaPrime->SetParameters(oldFuncParams);
    fConvolutedGausDeltaPrime->SetNpx(oldFuncNpx);
    fConvolutedGausDeltaPrime->SetRange(oldFuncRangeLow, oldFuncRangeUp);
    
    delete hInput;
    delete [] oldFuncParams;
    
    return kFALSE; 
  }
  
  if (fitResGauss->GetParams()[2] <= 0) {
    AliError("Not able to calculate parameters for the transition \"restricted gauss\" -> \"convoluted gauss\": Sigma of gauss fit <= 0!\n");
    fConvolutedGausDeltaPrime->SetParameters(oldFuncParams);
    fConvolutedGausDeltaPrime->SetNpx(oldFuncNpx);
    fConvolutedGausDeltaPrime->SetRange(oldFuncRangeLow, oldFuncRangeUp);
    
    delete hInput;
    delete [] oldFuncParams;
    
    return kFALSE;
  }
  
  // sigma correction factor
  fConvolutedGaussTransitionPars[1] = sigmaInput / fitResGauss->GetParams()[2];
  
  // 3.2 Now that sigma und lambda are determined, one can calculate mu by shifting the maximum to the desired position,
  // i.e. the maximum of the input distribution should coincide with that of the re-calculated distribution,
  // which is achieved by getting the same mu for the same sigma.
  // NOTE: For fixed lambda, the shift is proportional to sigma and independent of mu!
  // So, one can calculate the shift for an arbitrary fixed (here: typical)
  // sigma and then simply use this shift for any other sigma by scaling it correspondingly!!!
  
  // Mu shift correction:
  // Shift in mu (difference between mean of gauss and mean of convolution) is proportional to sigma!
  // Thus, choose a reference sigma (typical -> 0.05), such that for arbitrary sigma one can simple scale 
  // this shift correction with sigma / referenceSigma.
  fConvolutedGaussTransitionPars[0] = (fitResGauss->GetParams()[1] - muInput);
  
  fConvolutedGausDeltaPrime->SetParameters(oldFuncParams);
  fConvolutedGausDeltaPrime->SetNpx(oldFuncNpx);
  fConvolutedGausDeltaPrime->SetRange(oldFuncRangeLow, oldFuncRangeUp);
  
  delete hInput;
  delete [] oldFuncParams;

  return kTRUE;
}



//_____________________________________________________________________________
// This function is more or less a copy of AliAnalysisTaskPID::SetParamsForConvolutedGaus
Bool_t SetParamsForConvolutedGaus(Double_t gausMean, Double_t gausSigma) 
{
  // Set parameters for convoluted gauss using parameters for a pure gaussian.
  // If SetConvolutedGaussLambdaParameter has not been called before to initialise the translation parameters,
  // some default parameters will be used and an error will show up.
  
  if (fConvolutedGaussTransitionPars[1] < -998) {
    AliError("Transition parameters not initialised! Default parameters will be used. Please call SetConvolutedGaussLambdaParameter(...) before any calculations!");
    SetConvolutedGaussLambdaParameter(2.0);
    AliError(Form("Parameters set to:\n[0]: %f\n[1]: %f\n[2]: %f\n", fConvolutedGaussTransitionPars[0],           
                  fConvolutedGaussTransitionPars[1], fConvolutedGaussTransitionPars[2]));
  }
  
  Double_t par[fkConvolutedGausNPar];
  par[2] = fConvolutedGaussTransitionPars[2];
  par[1] = fConvolutedGaussTransitionPars[1] * gausSigma;
  // maxX perfectly shifts as par[0] (scaled by sigma) -> Can shift maxX so that it sits at the right place.
  par[0] = gausMean - fConvolutedGaussTransitionPars[0] * par[1] / fgkSigmaReferenceForTransitionPars;
  
  fConvolutedGausDeltaPrime->SetParameters(par);
  
  fConvolutedGausDeltaPrime->SetNpx(20); // Small value speeds up following algorithm (valid, since extrema far apart)

  // Accuracy of 10^-5 is enough to get 0.1% precise peak for MIPS w.r.t. to dEdx = 2000 of protons
  // (should boost up the algorithm, because 10^-10 is the default value!)
  Double_t maxX= fConvolutedGausDeltaPrime->GetMaximumX(TMath::Max(0.001, gausMean - 2. * gausSigma), 
                                                        gausMean + 6. * gausSigma, 1.0E-5);
  
  const Double_t maximum = fConvolutedGausDeltaPrime->Eval(maxX);
  const Double_t maximumFraction = maximum * fAccuracyNonGaussianTail;
  
  // Estimate lower boundary for subsequent search:
  Double_t lowBoundSearchBoundLow = TMath::Max(1e-4, maxX - 5. * gausSigma);
  Double_t lowBoundSearchBoundUp = maxX;
  
  Bool_t lowerBoundaryFixedAtZero = kFALSE;
  
  while (fConvolutedGausDeltaPrime->Eval(lowBoundSearchBoundLow) >= maximumFraction) {
    if (lowBoundSearchBoundLow <= 0) {
      // This should only happen to low dEdx particles with very few clusters and therefore large sigma, such that the gauss goes below zero deltaPrime
      if (maximum <= 0) { // Something is weired
        printf("Error generating signal: maximum is <= 0!\n");
        return kFALSE;
      }
      else {
        const Double_t valueAtZero = fConvolutedGausDeltaPrime->Eval(0);
        if (valueAtZero / maximum > 0.05) {
          // Too large fraction below zero deltaPrime. Signal generation cannot be reliable in this case
          printf("Error generating signal: Too large fraction below zero deltaPrime: convGauss(0) / convGauss(max) =  %e / %e = %e!\n",
                 valueAtZero, maximum, valueAtZero / maximum);
          return kFALSE;
        }
      }
      
      lowerBoundaryFixedAtZero = kTRUE;
      break;
    }
    
    lowBoundSearchBoundUp -= gausSigma;
    lowBoundSearchBoundLow -= gausSigma;
    
    if (lowBoundSearchBoundLow < 0) {
      lowBoundSearchBoundLow = 0;
      lowBoundSearchBoundUp += gausSigma;
    }
  }
  
  // Determine lower boundary inside estimated range. For small values of the maximum: Need more precision, since finer binning!
  Double_t rangeStart = lowerBoundaryFixedAtZero ? 0 :
                        fConvolutedGausDeltaPrime->GetX(maximumFraction, lowBoundSearchBoundLow, lowBoundSearchBoundUp, (maxX < 0.4) ? 1e-5 : 0.001);
  
  // .. and the same for the upper boundary
  Double_t rangeEnd = 0;
  // If distribution starts beyond upper boundary, everything ends up in the overflow bin. So, just reduce range and Npx to minimum
  if (rangeStart > fkDeltaPrimeUpLimit) {
    rangeEnd = rangeStart + 0.00001;
    fConvolutedGausDeltaPrime->SetRange(rangeStart,rangeEnd);
    fConvolutedGausDeltaPrime->SetNpx(4);
  }
  else {
    // Estimate upper boundary for subsequent search:
    Double_t upBoundSearchBoundUp = maxX + 5 * gausSigma;
    Double_t upBoundSearchBoundLow = maxX;
    while (fConvolutedGausDeltaPrime->Eval(upBoundSearchBoundUp) >= maximumFraction) {
      upBoundSearchBoundUp += gausSigma;
      upBoundSearchBoundLow += gausSigma;
    }
  
    //  For small values of the maximum: Need more precision, since finer binning!
    rangeEnd = fConvolutedGausDeltaPrime->GetX(maximumFraction, upBoundSearchBoundLow, upBoundSearchBoundUp, (maxX < 0.4) ? 1e-5 : 0.001);
    
    fConvolutedGausDeltaPrime->SetRange(rangeStart,rangeEnd);
    fConvolutedGausDeltaPrime->SetNpx(fDeltaPrimeAxis->FindFixBin(rangeEnd) - fDeltaPrimeAxis->FindFixBin(rangeStart) + 1);
    
    fConvolutedGausDeltaPrime->SetNpx(fDeltaPrimeAxis->FindFixBin(rangeEnd) - fDeltaPrimeAxis->FindFixBin(rangeStart) + 1);
  }
  
  
  return kTRUE;
}


//_____________________________________________________________________________
Int_t FitAsymmetricShapeHisto(TH1* h, Int_t fixPars, Double_t shiftPar = 0.002, Double_t sigmaPar = 0.9, Double_t lambdaPar = 0.5)
{
  if (!h) {
    std::cout << "Failed to load data histo!" << std::endl;
    return -1;
  }
  h->SetName("h");
  
  Bool_t isDelta = kFALSE; // kFALSE = isDeltaPrime
  
  if (h->GetXaxis()->GetBinLowEdge(1) < 0)
    isDelta = kTRUE;
  
  TCanvas* c = new TCanvas("c", "c",  100,10,1380,800);
  c->SetGridx(kTRUE);
  c->SetGridy(kTRUE);
  c->SetLogy(kTRUE);
  c->cd();
  
  h->SetDrawOption("e");
  h->DrawCopy("e");
  
  // Step 1: Fit distribution with restricted gaussian
  
  // Find maximum value of distribution (average over 3 bins for better stability)
  Int_t maxBin = h->GetMaximumBin();
  Double_t max = h->GetBinContent(maxBin);
  UChar_t usedBins = 1;
  if (maxBin > 1) {
    max += h->GetBinContent(maxBin - 1);
    usedBins++;
  }
  if (maxBin < h->GetNbinsX()) {
    max += h->GetBinContent(maxBin + 1);
    usedBins++;
  }
  max /= usedBins;
  
  
  // Set fit range to first/last bin below/above a certain fraction of the maximum
  // NOTE: The range (<-> fraction of maximum) should be the same
  // as for the spline and eta maps creation
  Double_t rangeStartSecondStep = 0.5;
  Double_t rangeEndSecondStep = 2.0;
  Double_t lowThreshold = h->GetXaxis()->GetBinLowEdge(h->FindFirstBinAbove(0.1 * max));
  Double_t highThreshold = h->GetXaxis()->GetBinUpEdge(h->FindLastBinAbove(0.1 * max));
  
  // Fit the restricted gaussian
  TFitResultPtr fitResGaussFirstStep = h->Fit("gaus", "NRS", "", lowThreshold, highThreshold);
  
  TFitResultPtr fitResGauss;
  
  // Do a second gaussian fit in the full range (or rather: sufficiently broad interval that contains more or less the whole distribution).
  // If the first fit succeeded, fix the mean value of the second fit to that of the first one and set for the amplitude and the width as initial values the results from the first fit.
  // If the first fit failed, just fit the gaussian in the full range without any restrictions for the parameters (fall back solution, in case the restricted range has caused some troubles)
  if ((Int_t)fitResGaussFirstStep == 0) {
    TF1 fGauss("fGauss", "[0]*TMath::Gaus(x, [1], [2], 1)", rangeStartSecondStep, rangeEndSecondStep);
    fGauss.SetParameter(0, fitResGaussFirstStep->GetParams()[0]);
    fGauss.SetParError(0, fitResGaussFirstStep->GetErrors()[0]);
    fGauss.SetParameter(2, fitResGaussFirstStep->GetParams()[2]);
    fGauss.SetParError(2, fitResGaussFirstStep->GetErrors()[2]);

    fGauss.FixParameter(1, fitResGaussFirstStep->GetParams()[1]);
    fitResGauss = h->Fit(&fGauss, "S+", "same", rangeStartSecondStep, rangeEndSecondStep);
  }
  else {
    fitResGauss = h->Fit("gaus", "RS+", "same", rangeStartSecondStep, rangeEndSecondStep);
  }
  
  const Int_t nPar = 4;
  TF1* fitFunc = new TF1("fitFunc", conv0, h->GetXaxis()->GetBinLowEdge(h->FindFirstBinAbove(20.)),
                         h->GetXaxis()->GetBinUpEdge(h->FindLastBinAbove(20.)), nPar);
  
  printf("Using params for pp @ 7 TeV!!!!!!!!!!!!\n");
  Double_t parAnalytical[nPar] = { h->Integral("width"),
                         fitResGauss->GetParams()[1] - 0.025553 * 0.046838 / 0.05, 
                         0.8945 * fitResGauss->GetParams()[2],
                         2.0 };
  
  // Maximum of fitFunc should sit at Mean(GaussFit) - 0.002 (empirically determined) -> Shift mean correspondingly
  fitFunc->SetParameters(par);
  
  fitFunc->SetRange(TMath::Max(0.001, fitResGauss->GetParams()[1] * 0.7), fitResGauss->GetParams()[1] * 1.3);
  Double_t maxX = fitFunc->GetMaximumX();
  if (maxX <= 0) {
    printf("Error: Maximum of fitFunc <= 0!\n");
    return -1;
  }
  
  // Maybe really SHIFT the mean of the convoluted gaus instead of multiply? Because the mean really only shifts the whole shape! -> Results ~equal
  //par[1] -= maxX - (fitResGauss->GetParams()[1] - 0.002);
  par[1] *= (fitResGauss->GetParams()[1] - 0.002) / maxX;
  
  if (fixPars >= 1) {
    if (fixPars == 2) {
      fitFunc->FixParameter(1, par[1]);
      fitFunc->FixParameter(2, par[2]);
    }
    fitFunc->FixParameter(3, par[3]);
  }
  
  
  fitFunc->SetLineColor(kBlue);
  
  TFitResultPtr fitRes = h->Fit(fitFunc, "RS+");
  
  if ((Int_t)fitRes == 0) {
    printf("Fit successful!\n\n");
    for (Int_t i = 0; i < nPar; i++)  {
      printf("par[%d] = %f +- %f\n", i, fitRes->GetParams()[i], fitRes->GetErrors()[i]);
    }
    printf("\nMaximum of distribution: %f\n", fitFunc->GetMaximumX());
    printf("\nchi^2/ndf = %f / %d = %f\n\n", fitRes->Chi2(), fitRes->Ndf(), fitRes->Ndf() > 0 ? fitRes->Chi2() / fitRes->Ndf() : -1);
  }
  else
    printf("Fit failed!\n");
  
  fitFunc->SetRange(h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
  
  
  Double_t maximumFraction = fitFunc->Eval(maxX) * 1e-8;

  fitFunc->SetNpx(20); // Small value speeds up following algorithm (valid, since extrema far apart)
  Double_t rangeStart = fitFunc->GetX(maximumFraction, TMath::Max(0.001, maxX - 30 * par[2]), maxX, 0.001);
  Double_t rangeEnd = fitFunc->GetX(maximumFraction, maxX, maxX + 30 * par[2], 0.001);
  printf("old: %e, %e, %e\n", maximumFraction, rangeStart, rangeEnd);
  fitFunc->SetRange(rangeStart,rangeEnd);
  fitFunc->SetNpx((rangeEnd - rangeStart) / h->GetXaxis()->GetBinWidth(1) + 1);
  
  fitFunc->Draw("same");
  
  
  
  
  
  
  TF1* fitFuncAnalytical = new TF1("fitFuncAnalytical", convolutedGaussAnalytical, rangeStart, rangeEnd, nPar);
  
  // Maximum of fitFunc should sit at Mean(GaussFit) - 0.002 (empirically determined) -> Shift mean correspondingly
  fitFuncAnalytical->SetParameters(parAnalytical);
  
  // Maximum shifts proportional to sigma and is linear in mu (~mean of gauss)
  // Maximum should be typically located within [gaussMean, gaussMean + 3 gaussSigma]. 
  // -> Larger search range for safety reasons (also: sigma and/or mean might be not 100% accurate).
  //fitFuncAnalytical->SetRange(TMath::Max(0.001, fitResGauss->GetParams()[1] * (1. - 2 * fitResGauss->GetParams()[2])),
  //                            fitResGauss->GetParams()[1] * (1. + 6 * fitResGauss->GetParams()[2]));
  Double_t maxXAnalytical = fitFuncAnalytical->GetMaximumX(TMath::Max(0.001, fitResGauss->GetParams()[1] * (1. - 2 * fitResGauss->GetParams()[2])),
                                                           fitResGauss->GetParams()[1] * (1. + 6 * fitResGauss->GetParams()[2]),
                                                           0.0001);
  if (maxXAnalytical <= 0) {
    printf("Error: Maximum of fitFuncAnalytical <= 0!\n");
    return -1;
  }
  
  // maxX perfectly shifts as par[1] -> Shift maxX to heuristically found value
  parAnalytical[1] -= maxXAnalytical - (fitResGauss->GetParams()[1] - 0.002);
  
  if (fixPars >= 1) {
    fitFuncAnalytical->FixParameter(3, lambdaPar);
    if (fixPars == 2) {
      fitFuncAnalytical->FixParameter(2, sigmaPar * fitResGauss->GetParams()[2]);
      fitFuncAnalytical->FixParameter(1, fitResGauss->GetParams()[1] - shiftPar * fitFuncAnalytical->GetParameter(2) / 0.05);
    }
  }
  
                                  
  Double_t maximumFractionAnalytical = fitFuncAnalytical->Eval(maxXAnalytical) * 1e-8;

  fitFuncAnalytical->SetNpx(20); // Small value speeds up following algorithm (valid, since extrema far apart)
  Double_t rangeStartAnalytical = fitFuncAnalytical->GetX(maximumFractionAnalytical,
                                                          TMath::Max(0.001, maxXAnalytical - 30 * parAnalytical[2]),
                                                          maxXAnalytical, 0.001);
  Double_t rangeEndAnalytical   = fitFuncAnalytical->GetX(maximumFractionAnalytical, maxXAnalytical, 
                                                          maxXAnalytical + 30 * parAnalytical[2], 0.001);
  printf("\nanalytical: %e, %e, %e\n", maximumFractionAnalytical, rangeStartAnalytical, rangeEndAnalytical);
  fitFuncAnalytical->SetRange(rangeStartAnalytical,rangeEndAnalytical);
  fitFuncAnalytical->SetNpx((rangeEndAnalytical - rangeStartAnalytical) / h->GetXaxis()->GetBinWidth(1) + 1);
  
  
  fitFuncAnalytical->SetLineColor(kMagenta);
  
  printf("\n\n");
  TFitResultPtr fitRes2 = h->Fit(fitFuncAnalytical, "RS+");
  
  if ((Int_t)fitRes2 == 0) {
    printf("Fit successful!\n\n");
    for (Int_t i = 0; i < nPar; i++)  {
      printf("par[%d] = %f +- %f\n", i, fitRes2->GetParams()[i], fitRes2->GetErrors()[i]);
    }
    printf("\nMaximum of distribution: %f\n", fitFuncAnalytical->GetMaximumX());
    printf("\nchi^2/ndf = %f / %d = %f\n\n", fitRes2->Chi2(), fitRes2->Ndf(), fitRes2->Ndf() > 0 ? fitRes2->Chi2() / fitRes2->Ndf() : -1);
    
    Double_t shift = (fitResGauss->GetParams()[1] - fitRes2->GetParams()[1]) * 0.05 / fitRes2->GetParams()[2];
    Double_t sigmaScale = fitRes2->GetParams()[2] / fitResGauss->GetParams()[2];
    printf("\n\nsigma_fit / sigma_gaus = %f\nshift (reference is 0.05) = %f\nreal shift = %f\ntail = sigma / lambda = %f\n\n", 
           sigmaScale,
           shift,
           shift * sigmaScale * fitResGauss->GetParams()[2] / 0.05,
           fitRes2->GetParams()[2] / fitRes2->GetParams()[3]);
  }
  else
    printf("Fit failed!\n");
  
  gRandom = new TRandom3(0);

  TH1D* hGenerated = new TH1D((TH1D&)*h);
  hGenerated->SetName("hGenerated");
  hGenerated->SetLineColor(kRed);
  hGenerated->SetMarkerColor(kRed);
  hGenerated->Reset();
  

  
  TStopwatch timer;
  timer.Start();
  {
    TH1* hProbDensity = fitFunc->GetHistogram();
  
    Double_t random = -1;
    for (Int_t i = 0; i < 20000000; i++)  {
      random = hProbDensity->GetRandom();
    }
  }
  timer.Stop();
  printf("Convoluted, non-analytical:\n");
  timer.Print();
  
  timer.Reset();
  timer.Start();
  {
    TH1* hProbDensity = fitFuncAnalytical->GetHistogram();
  
    Double_t random = -1;
    for (Int_t i = 0; i < 20000000; i++)  {
      random = hProbDensity->GetRandom();
    }
  }
  timer.Stop();
  printf("Analytical:\n");
  timer.Print();
  
  TH1* hProbDensity = fitFuncAnalytical->GetHistogram();
  
  Double_t random = -1;
  for (Int_t i = 0; i < 20000000; i++)  {
    random = hProbDensity->GetRandom();
    hGenerated->Fill(random);
    //hGenerated->Fill(fitFunc->GetRandom());// MUCH slower than using the binned version via the histogram
  }
  
  //timer.Stop();
  
  printf("Integral generated in range %e\n\n", hGenerated->Integral(hGenerated->GetXaxis()->FindBin(0.6), hGenerated->GetXaxis()->FindBin(2.0)));
  
  
  
  if (hGenerated->Integral() > 0) {
    hGenerated->Scale(h->Integral() / hGenerated->Integral());
    hGenerated->Draw("e same");
  }
  
  TH1D* hDiff = new TH1D((TH1D&)*hGenerated);
  hDiff->Add(h, -1);
  hDiff->Scale(-1);
  hDiff->Divide(h);
  
  new TCanvas();
  hDiff->Draw("e");
  
  return 0;
}
*/


//_____________________________________________________________________________
Double_t convolutedGaussAnalytical(const Double_t* xx, const Double_t* par)
{
  const Double_t amp = par[0];
  const Double_t mean = par[1];
  const Double_t sigma = par[2];
  const Double_t lambda = par[3];
  
  return amp*lambda/sigma*TMath::Exp(-lambda/sigma*(xx[0]-mean)+lambda*lambda*0.5)*0.5*TMath::Erfc((-xx[0]+mean+sigma*lambda)/sigma*OneOverSqrt2);
}


//_____________________________________________________________________________
TFitResultPtr FitAsymmetricShapeHisto(TH1* h, Double_t lambdaInitialValue = 2.0, Double_t lambdaLowerbound = 0.2, Double_t lambdaUpperBound = 10.0, Bool_t verbose = kTRUE)
{
  if (!h) {
    std::cout << "Failed to load data histo!" << std::endl;
    return -1;
  }
  
  TCanvas* c = 0x0;
  
  if (verbose) {
    h->SetName("h");

    c = new TCanvas("c", "c",  100,10,1380,800);
    c->Divide(1,2);
    c->GetPad(1)->SetGridx(kTRUE);
    c->GetPad(1)->SetGridy(kTRUE);
    c->GetPad(1)->SetLogy(kTRUE);
    c->cd(1);
  }
  
  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(20);
  h->SetDrawOption("e");
  if (verbose)
    h->Draw("e");
  
  const Double_t rangeStart = 0.5;
  const Double_t rangeEnd   = 2.0;
  const Int_t nPar = 4;
  
  // At first, fit a pure gaussian as reference:
  
  if (verbose)
    printf("Fitting pure gaussian as reference:\n\n");
  
  TFitResultPtr fitResPureGauss = h->Fit("gaus", "S+Q", "same", rangeStart, rangeEnd);
  
  if (verbose) {
    printf("\nchi^2/ndf = %f / %d = %f\n\n", fitResPureGauss->Chi2(), fitResPureGauss->Ndf(),
          fitResPureGauss->Ndf() > 0 ? fitResPureGauss->Chi2() / fitResPureGauss->Ndf() : -1);
    printf("\n\nNow fitting asymmetric shape:\n\n");
  }
  
  TF1* fitFuncAnalytical = new TF1("fitFuncAnalytical", convolutedGaussAnalytical, rangeStart, rangeEnd, nPar);
  
  const Double_t parsAnalytical[nPar] = {
    h->Integral("width"), // amplitude
    h->GetMean(),  // mean
    h->GetRMS(),   // sigma
    lambdaInitialValue // lambda
  };
  
  fitFuncAnalytical->SetParameters(parsAnalytical);
  fitFuncAnalytical->SetParLimits(1, 0.6, 1.4); // restrict mean around unity, where it should be if splines are correct
  fitFuncAnalytical->SetParLimits(2, 0.01, 0.2); // restrict sigma around typical relative resolution of 5-6%
  fitFuncAnalytical->SetParLimits(3, lambdaLowerbound, lambdaUpperBound); // restrict lambda
  
  // Set resolution to be consistent with binning of histogram
  fitFuncAnalytical->SetNpx((rangeEnd - rangeStart) / h->GetXaxis()->GetBinWidth(1) + 1);
  
  fitFuncAnalytical->SetLineColor(kBlue);
  fitFuncAnalytical->SetLineWidth(2);
  
  TStopwatch timer;
//   timer.Start();
  
  TFitResultPtr fitRes = h->Fit(fitFuncAnalytical, "RSM+Q");
  
  // For status codes, see e.g. https://sft.its.cern.ch/jira/browse/ROOT-4445 :
  // 4: MIGRAD error
  // 40 : MINOS error
  // 400: HESSE error
  // 4000: IMPROVE error
  // => 4 is fatal! 400 means that errors are not reliable. All others are not critical at all, e.g. 4000 just means that no better minimum
  // could be found, which doesn't matter, since the fit has obviously converged to the right minimum then!
  if (fitRes->IsValid()) {
    if (verbose) {
      printf("Fit successful! Status: %d\n\n", (Int_t)fitRes);
      for (Int_t i = 0; i < nPar; i++)
        printf("par[%d] = %f +- %f\n", i, fitRes->GetParams()[i], fitRes->GetErrors()[i]);
      
      printf("\nchi^2/ndf = %f / %d = %f\n\n", fitRes->Chi2(), fitRes->Ndf(), fitRes->Ndf() > 0 ? fitRes->Chi2() / fitRes->Ndf() : -1);
      
      printf("\n\n=> Extracted value for lambda is roughly %.1f!\n", fitRes->GetParams()[3]);
    }
  }
  else
    printf("Fit failed! Status = %d\n", (Int_t)fitRes);
  
  if (verbose) {
    timer.Stop();
    printf("\n");
    timer.Print();
    
//     Plot (data - fit) / data
    TH1D* hDiff = new TH1D((TH1D&)*h);
    hDiff->SetTitle("");
    hDiff->SetStats(0);
    hDiff->GetYaxis()->SetTitle("(data - fit) / data");
    hDiff->Add(fitFuncAnalytical, -1);
    hDiff->Divide(h);
    
    c->GetPad(2)->SetGridx(kTRUE);
    c->GetPad(2)->SetGridy(kTRUE);
    c->cd(2);
    
    hDiff->Draw("e");
  }
  
  return fitRes;
}



//_____________________________________________________________________________
Int_t FitAsymmetricShapeADV(Double_t lambdaInitialValue, Double_t lambdaLowerbound, Double_t lambdaUpperBound, 
                            TString path, TString fileName = "bhess_PIDetaAdv.root", TString listName = "bhess_PIDetaAdv")
{
  // etaAdv output as input
  
  TFile* f = TFile::Open(Form("%s/%s", path.Data(), fileName.Data()));
  if (!f)  {
    std::cout << "Failed to open file \"" << Form("%s/%s", path.Data(), fileName.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  THnSparse* hPIDdata = 0x0;
  TList* histList = 0x0;

  histList = (TList*)(f->Get(listName.Data()));
  if (!histList) {
    std::cout << "Failed to load list \"" << listName.Data() << "\"!" << std::endl;
    return -1;
  }
    
  // Extract the data histogram
  hPIDdata = dynamic_cast<THnSparse*>(histList->FindObject("hPIDdataAll"));
  if (!hPIDdata) {
    hPIDdata = dynamic_cast<THnSparse*>(f->Get(Form("%s/hPIDdataAll", listName.Data())));
    if (!hPIDdata)  {
      std::cout << "Failed to load data histo!" << std::endl;
      return -1;
    }
  }
  
  //TODO: CAREFUL: we cannot impose a cluster cut here. Since we have deltaPrime instead of the pull, tracks with different number of clusters
  // (= worse resolution) contribute via a broader shape. Strictly speaky, one then cannot fit the distribution with just one sigma!
  // To avoid this trouble, go a little bit away from MIP. The cluster number is then much more similar for all tracks. 
  // ALSO NOTE: One should not cut on the number of clusters close to MIP, since this introduces a dE/dx bias (the higher the dE/dx, the more likely it is to have a higher number of clusters).
  
  //TODO CAREFUL: Similar as for the clusters, the eta window must be narrow, since the resolution depends on eta!
  
  // TODO maybe take these as parameters of the function!
  // TODO Tune cuts and selections to have sufficient statistics in smallest possible pTPC and eta intervals. See also AliAnalysisTaskPIDetaAdv.*:
  // kMCpid: 1=e, 2=K, 3=#mu, 4=#pi, 5=p, 6=V0+TOF e, 7=V0 e, 8=V0 #pi, 9=V0 p
  // kSelectSpecies: 1=e, 2=K, 3=#pi, 4=p
  Double_t etaLow = 0.5;
  Double_t etaUp = 0.7;
  hPIDdata->GetAxis(kDeltaPrime)->SetRange(0,-1);
  hPIDdata->GetAxis(kMCpid)->SetRange(9,9);
  hPIDdata->GetAxis(kSelectSpecies)->SetRange(4,4);
  hPIDdata->GetAxis(kMultiplicity)->SetRangeUser(0.1, 499);
  // To avoid problems with bin edges, shift the desired values slightly
  hPIDdata->GetAxis(kEta)->SetRangeUser(etaLow + 1e-3, etaUp - 1e-3);
  hPIDdata->GetAxis(kPtpcInner)->SetRangeUser(0.8, 1.0);
  
  TH1D* hProj = (TH1D*)hPIDdata->Projection(kDeltaPrime);

  
  return FitAsymmetricShapeHisto(hProj, lambdaInitialValue, lambdaLowerbound, lambdaUpperBound);
}


//_____________________________________________________________________________
Int_t FitAsymmetricShapeTREE(Double_t lambdaInitialValue, Double_t lambdaLowerbound, Double_t lambdaUpperBound, 
                             TString path, TString fileName = "bhess_PIDetaTree.root", TString treeName = "fTree", Bool_t verbose = kTRUE, Double_t tanThetavalue= 0.0, Int_t tpcSignalNlow = 80, Int_t tpcSignalNUp = 160)
{
  // etaTree output as input
  
  TFile* f = 0x0;
  
  TString pathNameData = Form("%s/%s", path.Data(), fileName.Data());
  
  f = TFile::Open(pathNameData.Data());
  if (!f)  {
    std::cout << "Failed to open file \"" << pathNameData.Data() << "\"!" << std::endl;
    return -1;
  }
  
  // Extract the data tree
  TTree* tree = dynamic_cast<TTree*>(f->Get(treeName.Data()));
  if (!tree) {
    std::cout << "Failed to load data tree!" << std::endl;
    return -1;
  }
  
  tree->SetBranchStatus("*", 0); // Disable all branches
  tree->SetBranchStatus("pTPC", 1);
  tree->SetBranchStatus("dEdx", 1);
  tree->SetBranchStatus("dEdxExpected", 1);
  tree->SetBranchStatus("tanTheta", 1);
  tree->SetBranchStatus("tpcSignalN", 1);
  tree->SetBranchStatus("pidType", 1);
  
  // TODO maybe take these as parameters of the function!
  // TODO Tune cuts and selections to have sufficient statistics in smallest possible pTPC and tanTheta intervals. See also AliAnalysisTaskPIDetaTree.*
  Bool_t useV0s = kTRUE;
  TString pidSelectionString = useV0s ? Form("(pidType == %d || pidType == %d) && pTPC > 1.6 && pTPC < 2.0", kV0idNoTOF, kV0idPlusTOFaccepted)
                                      : Form("pidType == %d && pTPC > 0.84 && pTPC < 0.9", kTPCandTOFid);
  TString selectionString = Form("dEdxExpected > 0 && abs(tanTheta - %.2f) < 0.1 && %s && tpcSignalN > %d && tpcSignalN <= %d", tanThetavalue, 
                                 pidSelectionString.Data(), tpcSignalNlow, tpcSignalNUp);
  

  TH1D* hData  = new TH1D("hData", Form("%s;#Delta';Entries", selectionString.Data()), 200, 0.6, 1.6);
  tree->Project("hData", "dEdx/dEdxExpected", selectionString.Data());
  
  return FitAsymmetricShapeHisto(hData, lambdaInitialValue, lambdaLowerbound, lambdaUpperBound, verbose);
}


//_____________________________________________________________________________
//NOT USED at the moment
Int_t FitAsymmetricShapePIDFIT(Double_t lambdaInitialValue, Double_t lambdaLowerbound, Double_t lambdaUpperBound,
                               TString pathNameData = "finalCuts/MC_pp/7TeV/LHC10f6a/corrected/bhess_PID_Projections_newMethod_idSpectra.root",
                               TString histName = "hDeltaPiMC_0.80_Pt_0.85_species_3")
{
  // Inidiviual histos as input
  
  TFile* f = 0x0;
  
  f = TFile::Open(pathNameData.Data());
  if (!f)  {
    std::cout << "Failed to open file \"" << pathNameData.Data() << "\"!" << std::endl;
    return -1;
  }
  
  // Extract the data hist
  TH1* h = dynamic_cast<TH1*>(f->Get(histName.Data()));
  if (!h) {
    std::cout << "Failed to load data histo!" << std::endl;
    return -1;
  }
  
  return FitAsymmetricShapeHisto(h, lambdaInitialValue, lambdaLowerbound, lambdaUpperBound);
}


//_____________________________________________________________________________
Int_t FitAsymmetricShape(Double_t lambdaInitialValue, Double_t lambdaLowerbound, Double_t lambdaUpperBound, 
                         TString path, TString typeOfData, TString fileName = "", TString listOrTreeName = "", Double_t tanThetavalue = 0.5, Int_t tpcSignalNlow = 80, Int_t tpcSignalNhigh = 160)
{
  if (typeOfData.Contains("ADV", TString::kIgnoreCase))
    return FitAsymmetricShapeADV(lambdaInitialValue, lambdaLowerbound, lambdaUpperBound, path,
                                 fileName.IsNull() ? "bhess_PIDetaAdv.root" : fileName,
                                 listOrTreeName.IsNull() ? "bhess_PIDetaAdv" : listOrTreeName);
  else if (typeOfData.Contains("TREE", TString::kIgnoreCase))
    return FitAsymmetricShapeTREE(lambdaInitialValue, lambdaLowerbound, lambdaUpperBound, path, 
                                  fileName.IsNull() ? "bhess_PIDetaTree.root" : fileName,
                                  listOrTreeName.IsNull() ? "fTree" : listOrTreeName, kTRUE, tanThetavalue, tpcSignalNlow,  tpcSignalNhigh);
  else
    printf("Unknown type of input data: %s!\n", typeOfData.Data());
  
  return -1;
}


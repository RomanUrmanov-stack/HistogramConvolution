#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TH1D.h"

using namespace RooFit;

void roofit_project()
{
   //hadd all histograms if it wasn't already done
   if (gSystem->AccessPathName("K40_ZnCdWO4Ce_gain500_9_70_merged.root"))
   {
      gSystem->Exec("hadd -f K40_ZnCdWO4Ce_gain500_9_70_merged.root K40_ZnCdWO4Ce_gain500_9_70_v1.root K40_ZnCdWO4Ce_gain500_9_70_v2.root K40_ZnCdWO4Ce_gain500_9_70_v3.root");
   }

   //open root files
   TFile* f_pm = TFile::Open("K40.root"); //PM CsI(Tl) spectrum
   TH1F* hist_pm = (TH1F*)f_pm->Get("Integral")->Clone();

   TFile* f_pd = TFile::Open("K40_ZnCdWO4Ce_gain500_9_70_merged.root"); //PD Zn_CdWO4Ce merged spectrum
   TH1F* hist_pd = (TH1F*)f_pd->Get("Amplitude")->Clone();

   //reduce influence of statistical errors and background
   hist_pd->Rebin(8);

   //convolution function parameters
   double mean_conv_user = 0, sigma_conv_user = 6000;

   //first of all we're going to stretch and shift ZnWO4 spectrum
   int n_of_bins = hist_pd->GetXaxis()->GetNbins();
   int x_pm_range = hist_pm->GetXaxis()->GetXmax(), x_pd_range = hist_pd->GetXaxis()->GetXmax();
   int x_shift = 20200; 
   TH1F* new_hist_pd = new TH1F("new hist", "new hist",  n_of_bins, 0, x_pm_range);
   for (int i = 1; i <= n_of_bins; i++)
   {
      double y = hist_pd->GetBinContent(i);
      double x = hist_pd->GetXaxis()->GetBinCenter(i);
      double xnew = x_pm_range/x_pd_range * x + x_shift;
      new_hist_pd->Fill(xnew, y);
   }

   //initial data pdf
   RooRealVar x_pm("x_pm", "x_pm", 0, x_pm_range); 
   RooDataHist data_pm("data init", "intial data set", x_pm, hist_pm);
   RooHistPdf pdf_pm("pdf_pm", "pdf_pm", x_pm, data_pm);

   //let's find position of peak of interest
   int gauss_lower = 94 * pow(10,3), gauss_upper = 111 * pow(10,3); 
   RooRealVar mean_pm("mean_pm", "full absorption peak position", 100000, 0, x_pm_range);
   RooRealVar sigma_pm("sigma_pm", "full absorption peak sigma", 0, 10000);
   RooGaussian g_pm("g_pm","full absorption peak gauss", x_pm, mean_pm, sigma_pm);
   g_pm.fitTo(data_pm, Range(gauss_lower, gauss_upper));

   //convolution gauss
   RooRealVar mean_conv("mean conv", "convlution gauss mean", mean_conv_user, 0, 200 * pow(10,3));
   RooRealVar sigma_conv("sigma conv", "convlution gauss sigma", sigma_conv_user, 0, 10000);

   RooGaussian g_conv("g_conv", "convolution gauss", x_pm, mean_conv, sigma_conv);

   //convolution itself
   RooFFTConvPdf gxg("gxg", "init gauss (x) conv gauss", x_pm, pdf_pm, g_conv);

   //PD spectrum
   RooDataHist data_pd("data res", "low res data set", x_pm, new_hist_pd);
   RooHistPdf pdf_pd("tlr", "tlr", x_pm, data_pd);

   RooPlot* frame = x_pm.frame();
   
   //pdf_pm.plotOn(frame, LineColor(kRed)); //this is for testing convolution results
   //g_conv.plotOn(frame, LineColor(kBlack));
   gxg.plotOn(frame, LineColor(kGreen));
   pdf_pd.plotOn(frame, LineColor(kBlue));
   frame->Draw();

   //inverse transformation
   std::cout<< (mean_pm.getValV() - x_shift) * x_pd_range/x_pm_range << std::endl;
}
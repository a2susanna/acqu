#include <iostream>

void Plot_asymmetries() {

  //  gROOT->SetStyle("Plain");

  const Double_t k = 0.3905;
  const Double_t sintheta = TMath::Sin(81*TMath::DegToRad());

  Double_t asymm[30], asymmerr[30];
  Double_t beampol[30], beampolerr[30];
  Double_t mott[30], motterr[30];
  Double_t counter[30], counter2[30], err[30] = {0};
  Double_t day[22];
  Double_t asymm_minus[22], asymmerr_m[22];
  Double_t asymm_plus[22], asymmerr_p[22];
  Double_t beampol_minus[22], beampol_plus[22];
  Double_t mean[22];
  Double_t chi2r[30];

  Double_t t1, t2, t3, t4, t5;
  Double_t t6, t7, t8, t9, t10, t11, t12;
  Int_t i=0;

  ifstream infile2("MottMeasurements.dat");
  if (infile2.is_open()) {
    while (!infile2.eof()) {
      infile2 >> t6 >> t7 >> t8 >> t9 >> t10 >> t11 >> t12;
      day[i]           = t6;
      asymm_minus[i]   = t7;
      asymmerr_m[i]    = t8;
      beampol_minus[i] = t9;
      asymm_plus[i]    = t10;
      asymmerr_p[i]    = t11;
      beampol_plus[i]  = t12;
      i++;
    }
  }
  infile2.close();
  //   cout << "mott i " << i << endl;

  i = 0;
  ifstream infile("results_26Jan/fitvalues.dat");
  if (infile.is_open()) {
    while (!infile.eof()) {
      infile >> t1 >> t2 >> t3 >> t4 >> t5;
      asymm[i]      = t1;
      asymmerr[i]   = t2;
      beampol[i]    = t3*100;
      beampolerr[i] = t4*100;
      chi2r[i]      = t5;
      i++;
    }
  }
  infile.close();
  //   cout << "fitval i " << i << endl;


  for (Int_t j=0; j<i; j++) {
    if (beampol_minus[j] == 0) {
      mean[j] = (asymm_minus[j] - asymm_plus[j])/2.;
      mott[j] = mean[j]/(k*sintheta);
    }
    else 
      mott[j] = (beampol_minus[j] - beampol_plus[j])/2.;
    motterr[j] = TMath::Sqrt(asymmerr_m[j]*asymmerr_m[j] + asymmerr_p[j]*asymmerr_p[j]) / (2*k*sintheta);
    counter[j] = j+1;
    counter2[j] = j+0.5;
    //    cout << counter[j] << "\t c2 " << counter2[j] << "\t mott " << mott[j] << endl;
  }


  TCanvas * c2 = new TCanvas("c2","Moeller asymmetry", 800, 800);

  TGraphErrors * plot_asymm = new TGraphErrors(i-1, counter, asymm, err, asymmerr);
  plot_asymm->SetMarkerStyle(21);
  plot_asymm->SetMarkerSize(0.7);
  plot_asymm->SetMarkerColor(kBlue);
  plot_asymm->SetLineColor(kBlue);
  c2->cd();
  c2->SetGridx();
  c2->SetGridy();  
  plot_asymm->SetTitle("Moeller asymmetry distribution");
  plot_asymm->GetYaxis()->SetRangeUser(0., 0.05);
  plot_asymm->Draw("AP");
  plot_asymm->GetYaxis()->SetLabelSize(0.03);
  plot_asymm->GetXaxis()->SetTitle("Day");
  plot_asymm->GetXaxis()->SetTitleSize(0.03);
  plot_asymm->GetXaxis()->SetTitleOffset(1.8);


  TCanvas * c3 = new TCanvas("c3","Reduced chi2 distribution", 800, 800);

  TGraph * plot_chi2 = new TGraph(i-1, counter, chi2r);
  plot_chi2->SetMarkerStyle(20);
  plot_chi2->SetMarkerColor(kBlue);
  plot_chi2->SetLineColor(kBlue);
  c3->cd();
  c3->SetGridx();
  c3->SetGridy();  
  plot_chi2->GetYaxis()->SetRangeUser(0., 3);
  plot_chi2->Draw("AP");
  plot_chi2->SetTitle("#chi^{2}_{r} distribution");
  plot_chi2->GetYaxis()->SetLabelSize(0.03);
  plot_chi2->GetYaxis()->SetTitle("#chi^{2}_{r}");
  plot_chi2->GetYaxis()->SetTitleSize(0.03);
  plot_chi2->GetXaxis()->SetTitle("Day");
  plot_chi2->GetXaxis()->SetTitleSize(0.03);
  plot_chi2->GetXaxis()->SetTitleOffset(1.8);


  TCanvas * c = new TCanvas("c","Beam polarisation", 800, 800);
  TGraphErrors * plot_bp = new TGraphErrors(i-1, counter, beampol, err, beampolerr);
  plot_bp->SetMarkerStyle(21);
  plot_bp->SetMarkerSize(0.7);
  plot_bp->SetMarkerColor(kBlue);
  plot_bp->SetLineColor(kBlue);
  //   plot_bp->Draw("AP");

  TGraphErrors * plot_mott = new TGraphErrors(i, counter2, mott, err, motterr);
  plot_mott->SetMarkerStyle(20);
  plot_mott->SetMarkerSize(0.7);
  plot_mott->SetMarkerColor(kRed);
  plot_mott->SetLineColor(kRed);
  //   plot_mott->Draw("AP");


  TMultiGraph *mgres = new TMultiGraph();
  c->cd();
  c->SetGridx();
  c->SetGridy();
  mgres->Add(plot_bp);
  mgres->Add(plot_mott);
  mgres->Draw("AP");
  mgres->GetYaxis()->SetRangeUser(-100,0);
  mgres->SetTitle("Beam polarisation");
  mgres->GetYaxis()->SetTitle("Beam polarisation (%)");
  mgres->GetYaxis()->SetLabelSize(0.03);
  mgres->GetYaxis()->SetTitleSize(0.03);
  mgres->GetYaxis()->SetTitleOffset(1.5);
  mgres->GetXaxis()->SetTitle("Day");
  mgres->GetXaxis()->SetTitleSize(0.03);
  mgres->GetXaxis()->SetTitleOffset(1.8);
  const Int_t nbl = 15;
  const char *binlabels[nbl] = {"30/01","31/01","04/02","06/02","07/02","10/02","12/02","13/02","14/02","16/02","18/02","19/02","20/02","21/02","23/02"};

  mgres->GetXaxis()->SetBinLabel(5,binlabels[0]);
  mgres->GetXaxis()->SetBinLabel(12,binlabels[1]);
  mgres->GetXaxis()->SetBinLabel(18,binlabels[2]);
  mgres->GetXaxis()->SetBinLabel(25,binlabels[3]);
  mgres->GetXaxis()->SetBinLabel(31,binlabels[4]);
  mgres->GetXaxis()->SetBinLabel(38,binlabels[5]);
  mgres->GetXaxis()->SetBinLabel(44,binlabels[6]);
  mgres->GetXaxis()->SetBinLabel(51,binlabels[7]);
  mgres->GetXaxis()->SetBinLabel(57,binlabels[8]);
  mgres->GetXaxis()->SetBinLabel(63,binlabels[9]);
  mgres->GetXaxis()->SetBinLabel(70,binlabels[10]);
  mgres->GetXaxis()->SetBinLabel(76,binlabels[11]);
  mgres->GetXaxis()->SetBinLabel(83,binlabels[12]);
  mgres->GetXaxis()->SetBinLabel(90,binlabels[13]);
  mgres->GetXaxis()->SetBinLabel(96,binlabels[14]);

  plot_asymm->GetXaxis()->SetBinLabel(4,binlabels[0]);
  plot_asymm->GetXaxis()->SetBinLabel(10,binlabels[1]);
  plot_asymm->GetXaxis()->SetBinLabel(17,binlabels[2]);
  plot_asymm->GetXaxis()->SetBinLabel(23,binlabels[3]);
  plot_asymm->GetXaxis()->SetBinLabel(30,binlabels[4]);
  plot_asymm->GetXaxis()->SetBinLabel(36,binlabels[5]);
  plot_asymm->GetXaxis()->SetBinLabel(43,binlabels[6]);
  plot_asymm->GetXaxis()->SetBinLabel(50,binlabels[7]);
  plot_asymm->GetXaxis()->SetBinLabel(56,binlabels[8]);
  plot_asymm->GetXaxis()->SetBinLabel(63,binlabels[9]);
  plot_asymm->GetXaxis()->SetBinLabel(69,binlabels[10]);
  plot_asymm->GetXaxis()->SetBinLabel(76,binlabels[11]);
  plot_asymm->GetXaxis()->SetBinLabel(82,binlabels[12]);
  plot_asymm->GetXaxis()->SetBinLabel(89,binlabels[13]);
  plot_asymm->GetXaxis()->SetBinLabel(95,binlabels[14]);

  plot_chi2->GetXaxis()->SetBinLabel(4,binlabels[0]);
  plot_chi2->GetXaxis()->SetBinLabel(10,binlabels[1]);
  plot_chi2->GetXaxis()->SetBinLabel(17,binlabels[2]);
  plot_chi2->GetXaxis()->SetBinLabel(23,binlabels[3]);
  plot_chi2->GetXaxis()->SetBinLabel(30,binlabels[4]);
  plot_chi2->GetXaxis()->SetBinLabel(36,binlabels[5]);
  plot_chi2->GetXaxis()->SetBinLabel(43,binlabels[6]);
  plot_chi2->GetXaxis()->SetBinLabel(50,binlabels[7]);
  plot_chi2->GetXaxis()->SetBinLabel(56,binlabels[8]);
  plot_chi2->GetXaxis()->SetBinLabel(63,binlabels[9]);
  plot_chi2->GetXaxis()->SetBinLabel(69,binlabels[10]);
  plot_chi2->GetXaxis()->SetBinLabel(76,binlabels[11]);
  plot_chi2->GetXaxis()->SetBinLabel(82,binlabels[12]);
  plot_chi2->GetXaxis()->SetBinLabel(89,binlabels[13]);
  plot_chi2->GetXaxis()->SetBinLabel(95,binlabels[14]);

  for (Int_t ll=0; ll<nbl; ll++) {
    TLine *l = new TLine(counter2[ll], -100, counter2[ll], 0);
    l->SetLineStyle(3);
    c->cd();     l->Draw("same");
    l = new TLine(counter2[ll], 0, counter2[ll], 0.05);
    l->SetLineStyle(3);
    c2->cd();     l->Draw("same");
    l = new TLine(counter2[ll], 0, counter2[ll], 3.);
    l->SetLineStyle(3);
    c3->cd();     l->Draw("same");
  }

  c->cd();
  TLegend *legend = new TLegend(0.5,0.75,0.85,0.85);
  legend->AddEntry(plot_bp,"Moeller measurements", "lp");
  legend->AddEntry(plot_mott, "Mott measurements", "lp");
  legend->SetFillColor(0);
  legend->Draw("same");

  TFile *Output = new TFile("results_26Jan/Polarisation_results.root", "RECREATE");
  c->Write();
  c2->Write();
  c3->Write();
  Output->Write();
  Output->Close();
  
}




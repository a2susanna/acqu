#include <iostream>

const Int_t min_binx = 0;
const Int_t max_binx = 256;
const Int_t width = 5;
const Double_t lowlim = 50;
const Double_t upplim = 200;
//const Double_t Ebeam = 1557.; // MeV
const Double_t Ebeam = 450; // MeV
const Double_t m_ele = 0.511; // MeV
const Double_t foil_polar = 0.0808862;

ofstream fout;
ofstream fout2;

void MoellerPolarimeter(char *file_name) {

  // Creation of the left/right map
  Int_t map_all[5][11][6];
  Int_t map_left[5][11];
  Int_t map_right[5][15];
  Double_t right[4][15];
  Int_t left_offset[4] = {1,1,1,1};
  Int_t right_offset[4] = {16, 11, 3, 15};
  Double_t BeamPol[200];
  Double_t BeamPolErr[200];

  Int_t imax = 0;

  Double_t llower, rlower;
  for (Int_t ivup=0; ivup<4; ivup++) {
    switch (ivup) {
    case 0:
      llower = 128; rlower = 224; break;
    case 1: 
      llower = 96;  rlower = 256; break;
    case 2: 
      llower = 64;  rlower = 288; break;    
    case 3:
      llower = 0;   rlower = 320; break;
    }
    for (Int_t iright=0; iright<14; iright++) 
      map_right[ivup][iright] = rlower + iright + right_offset[ivup];
    for (Int_t ileft=0; ileft<10; ileft++) {
      map_left[ivup][ileft] = llower + ileft + left_offset[ivup];
      for (Int_t iiright=0; iiright<5; iiright++) 
	map_all[ivup][ileft][iiright] = map_right[ivup][9-ileft+iiright];
    }
  }

  // Reading and saving tagger calibration
  Double_t EgammaCal[353];
//   ifstream infile("/home/susanna/acqu_oct14/acqu/acqu_user/data/TagCal1557.dat");
  ifstream infile("/home/susanna/acqu_git_aug13/acqu_user/data-old-mine/Tagger/TagCal450.dat");
 
  if (infile.is_open()) {
    Double_t t1, t2, t3;
    for (Int_t i=1; i<353; i++) {
      infile >> t1 >> t2 >> t3;
      EgammaCal[i-1] = t2;
    }
  }
  infile.close();


  Double_t Esum[5][11][6];
  Double_t content[11][6][3];
  Double_t E_ele[5][11];
  Double_t theta_bar[5][11];
//   Double_t asymmetry[5][11];
//   Double_t errors[5][11];
  Double_t asymmetry[200] = {0.};
  Double_t errors[200] = {0.};
  Double_t a_zz[5][11];
  Double_t integral[6];
  Double_t height[6], pheight[11][6][3], baseline[11][6][3];
  Double_t mean[6];

  for (Int_t nvup=0; nvup<5; nvup++) {
    integral[nvup] = 0.;
    height[nvup] = 0.;
    for (Int_t nleft=0; nleft<10; nleft++) {
      for (Int_t nright=0; nright<6; nright++)
	for (Int_t hel=0; hel<2; hel++) {
	  content[nleft][nright][hel] = 0.;
	  pheight[nleft][nright][hel] = 0.;
	  baseline[nleft][nright][hel] = 0.;
	}
    }
  }


  TCanvas *c[9];
  TCanvas *c_fit[9];
  TCanvas *c_areas[9];

  TFile *infile_root = new TFile(file_name);
  char stringv[10] = {0};
//   strncpy(stringv, file_name+8, 8);
  strncpy(stringv, file_name+18,4);
  cout << stringv << endl;

  TString start_name = "TDC_Vup_Left_Pair_Hel_";
  TString und = "_";

  Int_t ic = 0;
  Int_t counter = 0;
  Int_t index_delta[5][11];
  
  TString asymm_name = "results_29Jan/asymmetry_values_fit_";
  asymm_name += stringv; asymm_name += ".dat";
  fout.open(asymm_name, ios::ate);

  for (Int_t nvup=0; nvup<4; nvup++) {
    fout << "Vuprom # " << nvup << endl;

    for (Int_t nleft=0; nleft<10; nleft++) {
      for (Int_t nright=0; nright<6; nright++)
	for (Int_t hel=0; hel<2; hel++) {
	  content[nleft][nright][hel] = 0.;
	  pheight[nleft][nright][hel] = 0.;
	  baseline[nleft][nright][hel] = 0.;
	}
    }

    for (Int_t hel=0; hel<2; hel++) {  

      TString canvas_name = "Vuprom_"; canvas_name += nvup; 
      canvas_name += "_Helicity_"; canvas_name += hel;
      c[ic] = new TCanvas(canvas_name, canvas_name, 1600, 1000);
      c[ic]->Divide(10, 5);
      TString canvas_fit_name = "Vuprom_"; canvas_fit_name += nvup;
      canvas_fit_name += "_Helicity_"; canvas_fit_name += hel;
      canvas_fit_name += "_fit";
      c_fit[ic] = new TCanvas(canvas_fit_name, canvas_fit_name, 1600, 1000);
      c_fit[ic]->Divide(10, 5);
      TString areas_name = "Integrals_vuprom"; areas_name += nvup;
      areas_name += "_Helicity"; areas_name += hel;
      c_areas[ic] = new TCanvas(areas_name, areas_name, 1600, 1000);
      c_areas[ic]->Divide(5,2);

      for (Int_t nleft=0; nleft<10; nleft++) {
// 	fout << "Left # " << nleft << endl;

	E_ele[nvup][nleft] = Ebeam - EgammaCal[map_left[nvup][nleft]];
	Double_t y = TMath::Sqrt(Ebeam - E_ele[nvup][nleft]);
	Double_t x = TMath::Sqrt(E_ele[nvup][nleft] - m_ele);
	theta_bar[nvup][nleft] = 2*TMath::ATan2(y,x);
	Double_t a_xx = - pow(TMath::Sin(theta_bar[nvup][nleft]),4) / pow((4 - TMath::Sin(theta_bar[nvup][nleft])*TMath::Sin(theta_bar[nvup][nleft])), 2);
  	a_zz[nvup][nleft] = a_xx * (8 - TMath::Sin(theta_bar[nvup][nleft])*TMath::Sin(theta_bar[nvup][nleft]))/(TMath::Sin(theta_bar[nvup][nleft])*TMath::Sin(theta_bar[nvup][nleft]));
 	// 	cout << "############## " << E_ele[nvup][nleft] << "\t" << theta_bar[nvup][nleft] << "\t" << a_zz[nvup][nleft] << endl;
	Double_t Esum_arr[5];
	Double_t max_int = 0;
	Double_t DeltaEsum = 1000;

	for (Int_t nright=0; nright<5; nright++) {

	  Double_t E_ele2 = Ebeam - EgammaCal[map_all[nvup][nleft][nright]];
	  Esum[nvup][nleft][nright] = E_ele[nvup][nleft] + E_ele2;
	  Esum_arr[nright] = Esum[nvup][nleft][nright];

	  integral[nright] = -999;
	  TString hname = start_name;
	  hname += nvup; hname += und; hname += nleft; hname += und;
	  hname += nright; hname += und; hname += hel;
	  
 	  TH1* histo = (TH1*) gROOT->FindObject(hname);
	  Int_t npad = nleft+1 + 10*nright;
 	  c[ic]->cd(npad);
 	  histo->Draw();
 	  height[nright] = histo->GetMaximum();
 	  mean[nright] = histo->GetMaximumBin(20,230,0);
  	  if (height[nright]<10) continue;
 	  TF1 *f = new TF1("f", "pol0", lowlim, upplim);
 	  histo->Fit(f,"Rq");
	  Double_t par0;
	  par0 = f->GetParameter(par0);
	  TF1 *f2 = new TF1("f2", "pol0(0) + gaus(1)", lowlim, upplim);
	  f2->SetParameters(par0, height[nright]-par0, mean[nright], width);
	  histo->Fit(f2, "Rq");
	  Double_t par[4];
	  f2->GetParameters(par);
	  baseline[nleft][nright][hel] = par[0];
 	  TF1 *fpeak = new TF1("fpeak", "gaus", 50, 200);
 	  fpeak->SetParameters(par[1], par[2], par[3]);
	  pheight[nleft][nright][hel] = fpeak->GetMaximum()/baseline[nleft][nright][hel];
	  fout << "L " << nleft << ", R " << nright << ", hel " << hel << "\t max " << fpeak->GetMaximum() << "\t baseline " << baseline[nleft][nright][hel] << "\t h " << pheight[nleft][nright][hel] << endl;

 	  c_fit[ic]->cd(npad);
 	  fpeak->Draw();
 	  integral[nright] = fpeak->Integral(par[2]-20, par[2]+20);
// 	  if (TMath::Abs(integral[nright] - fpeak->Integral(50,200)) > 10)
// 	    integral[nright] = -999;
	  if (integral[nright] > max_int) 
	    max_int = integral[nright];
	  if (TMath::Abs(Esum_arr[nright] - Ebeam) < DeltaEsum) {
	    DeltaEsum = TMath::Abs(Esum_arr[nright] - Ebeam);
	    index_delta[nvup][nleft] = nright;
	  }
//  	  fout << "Right # " << nright << "\t Esum = " << Esum[nvup][nleft][nright] << "\t N = " << integral[nright] << "\t h = " << par[1] << "\tDeltaE = " << TMath::Abs(Esum_arr[nright] - Ebeam) << endl;

	  content[nleft][nright][hel] = integral[nright];
	  if (content[nleft][nright][hel]>imax)
	    imax = content[nleft][nright][hel];
	} // right
	
	TGraph *areas = new TGraph(5, Esum_arr, integral);
	c_areas[ic]->cd(nleft+1);
	c_areas[ic]->cd(nleft+1)->SetGridx();
	c_areas[ic]->cd(nleft+1)->SetGridy();
// 	areas->Draw("AP*");
	areas->Draw("APB");
	areas->SetMarkerSize(1.0);
	areas->SetMarkerStyle(5);
	areas->SetMarkerColor(kBlue);
	areas->SetFillColor(kBlue);
	areas->SetFillStyle(3001);
	areas->GetYaxis()->SetRangeUser(0, max_int+1000);
// 	TLine *l = new TLine(1557, 0, 1557, max_int+1000);
 	TLine *l = new TLine(450, 0, 450, max_int+1000);
	l->SetLineColor(kRed);
	l->SetLineWidth(2);
	l->Draw("same");
      } // nleft

      ic++;
    } // hel

    // CHECK!!!

    for (Int_t nleft=0; nleft<10; nleft++) {
      for (Int_t nright=index_delta[nvup][nleft]; nright < 5; nright++) {
	if ((content[nleft][nright][0] + content[nleft][nright][1]) != 0) {
 	  if (pheight[nleft][nright][0]<1.5 || pheight[nleft][nright][1]<1.5)
 	    continue; 
// 	  if (content[nleft][nright][0]<15000 || content[nleft][nright][1] < 15000)
	  if (content[nleft][nright][0]<0.2*imax || content[nleft][nright][1]<0.2*imax)
	    continue;
	  asymmetry[counter] = TMath::Abs((content[nleft][nright][1] - content[nleft][nright][0]))/(content[nleft][nright][0] + content[nleft][nright][1]);
	  errors[counter]    = 1./TMath::Sqrt(content[nleft][nright][0] + content[nleft][nright][1]);
 	  BeamPol[counter] = asymmetry[counter]/(a_zz[nvup][nleft]*foil_polar*TMath::Sin(65*TMath::DegToRad()));
 	  BeamPolErr[counter] = -errors[counter]/(foil_polar*a_zz[nvup][nleft]*TMath::Sin(65*TMath::DegToRad()));
	  cout << BeamPol[counter] << " +/- " << BeamPolErr[counter] << endl;
	  //	  cout << "errors " << errors[counter] << "\t beampolerr " << BeamPolErr[counter] << "\t or " << TMath::Sqrt(errors[counter]**2+0.00209858**2+0.01) << endl;
	  //	  fout << "######## L " << nleft << "\t R " << nright << "\t" << content[nleft][nright][0] << "\t" << content[nleft][nright][1] <<  "\t th_bar " << theta_bar[nvup][nleft]*TMath::RadToDeg() << "\t a_zz " << a_zz[nvup][nleft] << "\t A = " << asymmetry[counter] << "\t AErr = " << errors[counter] <<  "\t BP = " << BeamPol[counter] << "\t BPErr = " << BeamPolErr[counter] <<  endl;
	  fout << "######## L " << nleft << "\t R " << nright << "\t N0 = " << content[nleft][nright][0] << "\t N1 = " << content[nleft][nright][1] << "\t h0 = " << pheight[nleft][nright][0] << "\t h1 = " << pheight[nleft][nright][1] << "\t A = " << asymmetry[counter] << "\t AErr = " << errors[counter] <<  "\t BP = " << BeamPol[counter] << "\t BPErr = " << BeamPolErr[counter] <<  endl;
	  counter++;
	}
      } // nright
    } // nleft
   
  } // nvup
  
  fout.close();
  //   infile_root->Close();

  Double_t zero[200] = {0.};
  Double_t order[200];
  for (Int_t i=0; i<=counter; i++)
    order[i] = i+1;

  TString fitname = "results_29Jan/fitvalues_";
  fitname += stringv; fitname += ".dat";
  fout2.open(fitname, ios::app);
      
  TGraphErrors * asymm = new TGraphErrors(counter, order, asymmetry, zero, errors);

  TString casymm_name = "Asymmetry_Moeller_";
  casymm_name += stringv;
  TCanvas *casymm = new TCanvas(casymm_name, casymm_name, 1000, 800);
  casymm->SetGridx();
  casymm->SetGridy();
  asymm->SetMarkerStyle(20);
  asymm->SetMarkerSize(0.9) ; 
  asymm->Draw("AP");
  asymm->SetTitle(casymm_name);
  asymm->GetXaxis()->SetTitleSize(0.05);
  asymm->GetXaxis()->SetTitleOffset(0.9);
  asymm->GetYaxis()->SetTitle("Moeller Asymmetry");
  asymm->GetYaxis()->SetTitleSize(0.05);

  TF1 * fasym = new TF1("fasym","pol0", 0, counter);
  asymm->Fit(fasym, "R");
  Double_t par_asym, parerr_asym;
  par_asym = fasym->GetParameter(par_asym);
  parerr_asym = fasym->GetParError(parerr_asym);

  TGraphErrors * beampol = new TGraphErrors(counter, order, BeamPol, zero, BeamPolErr);

  TString cbeampol_name = "Beampol"; cbeampol_name += stringv;
  TCanvas *cbeampol = new TCanvas(cbeampol_name, cbeampol_name, 1000, 800);
  cbeampol->SetGridx();
  cbeampol->SetGridy();
  beampol->SetMarkerStyle(20);
  beampol->SetMarkerSize(0.9) ; 
  beampol->Draw("AP");
  beampol->SetTitle(cbeampol_name);
  beampol->GetXaxis()->SetTitleSize(0.05);
  beampol->GetXaxis()->SetTitleOffset(0.9);
  beampol->GetYaxis()->SetTitle("Beam polarisation");
  beampol->GetYaxis()->SetRangeUser(-1,1);
  beampol->GetYaxis()->SetTitleSize(0.05);

  TF1 * fbeampol = new TF1("fbeampol", "pol0", 0, counter);
  beampol->Fit(fbeampol, "R");
  Double_t par_beampol, parerr_beampol;
  par_beampol = fbeampol->GetParameter(par_beampol);
  parerr_beampol = fbeampol->GetParError(parerr_beampol);
  Double_t chi2 = fbeampol->GetChisquare();
  Int_t ndf = fbeampol->GetNDF();
  fout2 << par_asym << "\t" << parerr_asym << "\t" << par_beampol << "\t" << parerr_beampol << "\t" << chi2/ndf << endl;
  fout2.close();

  cout << "Integral MAX " << imax << endl;

  TString outname = "results_29Jan/Output_asymm_fit_";
  outname += stringv; outname += ".root";
  TFile *Output = new TFile(outname, "RECREATE");
  
  for (Int_t ic=0; ic<8; ic++) {
    c[ic]->Write();
    c_fit[ic]->Write();
    c_areas[ic]->Write();
  }
  casymm->Write();
  cbeampol->Write();

  Output->Write();
  Output->Close();

}

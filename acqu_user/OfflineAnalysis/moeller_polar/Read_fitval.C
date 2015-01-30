#include <iostream>

void Read_fitval() {

  TString file_name;
  Double_t 
    t1 = 0.,
    t2 = 0.,
    t3 = 0.,
    t4 = 0.,
    t5 = 0.;

  ofstream outfile;
  outfile.open("fitvalues.dat", ios::app);

  for (int i=1; i<15; i++) {

    switch (i) {
    case 1:  file_name = "fitvalues_30-31Jan.dat"; break;
    case 2:  file_name = "fitvalues_31-04Feb.dat"; break;
    case 3:  file_name = "fitvalues_04-06Feb.dat"; break;
    case 4:  file_name = "fitvalues_06-07Feb.dat"; break;
    case 5:  file_name = "fitvalues_07-10Feb.dat"; break;
    case 6:  file_name = "fitvalues_10-12Feb.dat"; break;
    case 7:  file_name = "fitvalues_12-13Feb.dat"; break;
    case 8:  file_name = "fitvalues_13-14Feb.dat"; break;
    case 9:  file_name = "fitvalues_14-16Feb.dat"; break;
    case 10: file_name = "fitvalues_16-18Feb.dat"; break;
    case 11: file_name = "fitvalues_18-19Feb.dat"; break;
    case 12: file_name = "fitvalues_19-20Feb.dat"; break;
    case 13: file_name = "fitvalues_20-21Feb.dat"; break;
    case 14: file_name = "fitvalues_21-23Feb.dat"; break;
    }
    
    ifstream infile(file_name);
    if (infile.is_open()) 
      infile >> t1 >> t2 >> t3 >> t4 >> t5;
      
    outfile << t1 << "\t" << t2 << "\t" << t3 << "\t" << t4 << "\t" << t5 << endl;
    
  }

  outfile << endl;
  outfile.close();
}

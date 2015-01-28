#include <iostream>
#include <utility>
#include <vector>
#include <fstream>
      
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"

#include "TH2.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
  
#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;
 
//------------------------------------------------------------------------------

// Analysis macro
#include "EventSelections.C"

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char *appName = "EventSelections";

  if(argc != 3)
  {
    cout << " Usage : " << appName << " input_file output_file" << endl;
    cout << " Usage : " << appName << " input_file output_file" << endl;
    cout << " input_file - input file in ROOT format ('Delphes' tree)" << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  TString inputFile(argv[1]);
  TString outputFile(argv[2]);

//------------------------------------------------------------------------------

// Call the macros main function (in EventSelections.C)

  EventSelections(inputFile,outputFile);

//------------------------------------------------------------------------------

}



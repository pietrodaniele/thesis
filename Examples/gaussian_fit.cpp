#include<cstdio>
#include<cstdlib>
#include<iostream>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TH1F.h"

using namespace std;
using namespace ROOT;

int main() {
    // defining objetcs
    TCanvas *c1 = new TCanvas();
    TH1F *hist = new TH1F("m_yy","m_yy distribution",100,80,140);
    TFile *input = new TFile("../Data/PowhegPy8_NNLOPS_ggH110.root","read"); // reading the 110 Mev file
    TTree *tree = (TTree*)input->Get("myCatNtuple;1"); // loading the TTree file
    int entries = tree->GetEntries(); // number of rows

    // loading the columns
    float m_yy;

    tree->SetBranchAddress("m_yy", &m_yy);

    for(int i=0; i<entries; i++){
      tree->GetEntry(i);
      hist->Fill(m_yy/1000.);
      // cout << m_yy << endl;
    }

    // drawing datas
    hist->Draw();
    c1->cd();

    return 0;
}

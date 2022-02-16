#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit ;

// functions
RooDataSet *read_data(int, RooRealVar&);
vector<double> param_lin(double *, double *, double *, double *, string);
void plot_param(double [4], double [4], string);
void plot_results(string, RooFormulaVar, RooFormulaVar, RooRealVar, RooRealVar, RooRealVar, RooRealVar,
                  RooRealVar, RooDataSet, RooSimultaneous, RooCategory);


// mu and width functions of m_h
// cutFlow > 13 applied


void sim_fit_125removed(){
  // all parameters obtained from previous fit
  double masses[4] = {110,125,130,140};
  double masses_error[4] = {0,0,0,0};
  double mu[4] = {110.122045, 125.120605, 130.119203, 140.128386};
  double mu_error[4] = {0.031419, 0.032085, 0.034132, 0.041798};
  double sigma[4] = {1.683887, 1.788779, 1.823120, 1.906895};
  double sigma_error[4] = {0.037484, 0.037591, 0.040010, 0.048537};
  double a1[4] = {1.737295, 1.723697, 1.739125, 1.748865};
  double a1_error[4] = {0.123130, 0.116492, 0.122213, 0.144540};
  double p1[4] = {4.825106, 4.912847, 4.751721, 4.606459};
  double p1_error[4] = {1.078196, 1.066697, 1.061026, 1.199767};
  double a2[4] = {1.602378, 1.563756, 1.511050, 1.628135};
  double a2_error[4] = {0.145021, 0.130597, 0.126748, 0.178149};
  double p2[4] = {15.929337, 21.855431, 30.468083, 18.320550};
  double p2_error[4] = {9.688985, 16.271775, 28.915310, 15.547188};

  // linear fit results ==> (a,b) from y = a + b*x
  vector<double> mu_par = param_lin(masses,masses_error,mu,mu_error,"mu");
  vector<double> sigma_par = param_lin(masses,masses_error,sigma,sigma_error,"sigma");
  vector<double> a1_par = param_lin(masses,masses_error,a1,a1_error,"a1");
  vector<double> p1_par = param_lin(masses,masses_error,p1,p1_error,"p1");
  vector<double> a2_par = param_lin(masses,masses_error,a2,a2_error,"a2");
  vector<double> p2_par = param_lin(masses,masses_error,p2,p2_error,"p2");

  // creating a single mass range
  RooRealVar m("m","mass",80,170);
  // -------------------------------------------------------------

  // Create observables
  int intervall = 30;
  int mass110 = 110;
  int mass125 = 125;
  int mass130 = 130;
  int mass140 = 140;
  RooRealVar A_mu("A_mu","Bias mu",mu_par[0]);
  RooRealVar B_mu("B_mu","Ang coeff mu",mu_par[1]);
  RooRealVar A_width("A_width","Bias width",sigma_par[0]);
  RooRealVar B_width("B_width","Ang coeff width",sigma_par[1]);

  // dataset
  RooDataSet * dataset110 = read_data(110,m);
  RooDataSet * dataset125 = read_data(125,m);
  RooDataSet * dataset130 = read_data(130,m);
  RooDataSet * dataset140 = read_data(140,m);
  //  110 model
  // variables 110 GeV
  // fit_parameters
  // ---------------------------------------------------------------------------
  RooRealVar m_h_110("m_h_110", "Higgs mass", mass110);
  // ---------------------------------------------------------------------------
  RooFormulaVar mean_110("mu_110","@0 + @1*m_h_110",RooArgList(A_mu,B_mu,m_h_110));
  // ---------------------------------------------------------------------------
  RooFormulaVar width_110("width_110","@0 + @1*m_h_110",RooArgList(A_width,B_width,m_h_110));
  // ---------------------------------------------------------------------------
  RooRealVar a1_110("a1_110","a1",1,0.0001,100.);
  RooRealVar p1_110("p1_110","p1",1,0.0001,100.);
  RooRealVar a2_110("a2_110","a2",1,0.0001,100.);
  RooRealVar p2_110("p2_110","p2",1,0.0001,100.);
  RooCrystalBall dcbPdf_110("dcbPdf_110","DoubleSidedCB",m,mean_110,width_110,a1_110,p1_110,a2_110,p2_110);

  //  125 model
  // variables 125 GeV
  // fit_parameters
  // ---------------------------------------------------------------------------
  RooRealVar m_h_125("m_h_125", "Higgs mass", mass125);
  // ---------------------------------------------------------------------------
  RooFormulaVar mean_125("mu_125","@0 + @1*m_h_125",RooArgList(A_mu,B_mu,m_h_125));
  // ---------------------------------------------------------------------------
  RooFormulaVar width_125("width_125","@0 + @1*m_h_125",RooArgList(A_width,B_width,m_h_125));
  // ---------------------------------------------------------------------------
  RooRealVar a1_125("a1_125","a1",1,0.0001,100.);
  RooRealVar p1_125("p1_125","p1",1,0.0001,100.);
  RooRealVar a2_125("a2_125","a2",1,0.0001,100.);
  RooRealVar p2_125("p2_125","p2",1,0.0001,100.);
  RooCrystalBall dcbPdf_125("dcbPdf_125","DoubleSidedCB",m,mean_125,width_125,a1_125,p1_125,a2_125,p2_125);

  //  130 model
  // variables 130 GeV
  // fit_parameters
  // ---------------------------------------------------------------------------
  RooRealVar m_h_130("m_h_130", "Higgs mass", mass130);
  // ---------------------------------------------------------------------------
  RooFormulaVar mean_130("mu_130","@0 + @1*m_h_130",RooArgList(A_mu,B_mu,m_h_130));
  // ---------------------------------------------------------------------------
  RooFormulaVar width_130("width","@0 + @1*m_h_130",RooArgList(A_width,B_width,m_h_130));
  // ---------------------------------------------------------------------------
  RooRealVar a1_130("a1_130","a1",1,0.0001,100.);
  RooRealVar p1_130("p1_130","p1",1,0.0001,100.);
  RooRealVar a2_130("a2_130","a2",1,0.0001,100.);
  RooRealVar p2_130("p2_130","p2",1,0.0001,100.);
  RooCrystalBall dcbPdf_130("dcbPdf_130","DoubleSidedCB",m,mean_130,width_130,a1_130,p1_130,a2_130,p2_130);

  //  140 model
  // variables 140 GeV
  // fit_parameters
  // ---------------------------------------------------------------------------
  RooRealVar m_h_140("m_h_140", "Higgs mass", mass140);
  // ---------------------------------------------------------------------------
  RooFormulaVar mean_140("mu_140","@0 + @1*m_h_140",RooArgList(A_mu,B_mu,m_h_140));
  // ---------------------------------------------------------------------------
  RooFormulaVar width_140("width","@0 + @1*m_h_140",RooArgList(A_width,B_width,m_h_140));
  // ---------------------------------------------------------------------------
  RooRealVar a1_140("a1_140","a1",1,0.0001,100.);
  RooRealVar p1_140("p1_140","p1",1,0.0001,100.);
  RooRealVar a2_140("a2_140","a2",1,0.0001,100.);
  RooRealVar p2_140("p2_140","p2",1,0.0001,100.);
  RooCrystalBall dcbPdf_140("dcbPdf_140","DoubleSidedCB",m,mean_140,width_140,a1_140,p1_140,a2_140,p2_140);

  // C r e a t e   i n d e x   c a t e g o r y   a n d   j o i n   s a m p l e s
  // ---------------------------------------------------------------------------

  // Define category to distinguish 110,125,130 and 140 GeV samples events
  RooCategory sample("sample","sample") ;
  sample.defineType("myy_110") ;
  sample.defineType("myy_125") ;
  sample.defineType("myy_130") ;
  sample.defineType("myy_140") ;

  // Construct combined dataset in (x,sample)
  RooDataSet combData("combData","combined data",m,Index(sample),Import("myy_110",*dataset110),
                      Import("myy_125",*dataset125),Import("myy_130",*dataset130),Import("myy_140",*dataset140)) ;



  // C o n s t r u c t   a   s i m u l t a n e o u s   p d f   i n   ( x , s a m p l e )
  // -----------------------------------------------------------------------------------

  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;

  // Associate model with the physics state and model_ctl with the control state
  simPdf.addPdf(dcbPdf_110,"myy_110") ;
  //simPdf.addPdf(dcbPdf_125,"myy_125") ;
  simPdf.addPdf(dcbPdf_130,"myy_130") ;
  simPdf.addPdf(dcbPdf_140,"myy_140") ;



  // P e r f o r m   a   s i m u l t a n e o u s   f i t
  // ---------------------------------------------------

  // Perform simultaneous fit of model to data and model_ctl to data_ctl
  simPdf.fitTo(combData) ;

  cout << mean_110.getVal() << endl;
  cout << mean_130.getVal() << endl;


  // P l o t   m o d e l   s l i c e s   o n   d a t a    s l i c e s
  // ----------------------------------------------------------------

  /*// Make a frame for the 110 GeV sample
  RooPlot* frame1 = m.frame(Title("110 GeV m_{h} distribution")) ;
  combData.plotOn(frame1,Cut("sample==sample::myy_110")) ;
  simPdf.plotOn(frame1,Slice(sample,"myy_110"),ProjWData(sample,combData)) ;
  // simPdf.plotOn(frame1,Slice(sample,"myy_110"),Components("dcbPdf_110"),ProjWData(sample,combData),LineStyle(kDashed)) ;

  // Make a frame for the 110 GeV sample
  RooPlot* frame2 = m.frame(Title("125 GeV m_{h} distribution")) ;
  combData.plotOn(frame2,Cut("sample==sample::myy_125")) ;
  simPdf.plotOn(frame2,Slice(sample,"myy_125"),ProjWData(sample,combData)) ;
  // simPdf.plotOn(frame1,Slice(sample,"myy_110"),Components("dcbPdf_110"),ProjWData(sample,combData),LineStyle(kDashed)) ;

  // Make a frame for the 130 GeV sample
  RooPlot* frame3 = m.frame(Title("130 GeV m_{h} distribution")) ;
  combData.plotOn(frame3,Cut("sample==sample::myy_130")) ;
  simPdf.plotOn(frame3,Slice(sample,"myy_130"),ProjWData(sample,combData)) ;
  //simPdf.plotOn(frame2,Slice(sample,"myy_130"),Components("dcbPdf_130"),ProjWData(sample,combData),LineStyle(kDashed)) ;

  // Make a frame for the 130 GeV sample
  RooPlot* frame4 = m.frame(Title("140 GeV m_{h} distribution")) ;
  combData.plotOn(frame4,Cut("sample==sample::myy_140")) ;
  simPdf.plotOn(frame4,Slice(sample,"myy_140"),ProjWData(sample,combData)) ;
  //simPdf.plotOn(frame2,Slice(sample,"myy_130"),Components("dcbPdf_130"),ProjWData(sample,combData),LineStyle(kDashed)) ;
  */

  plot_results("110",mean_110,width_110,a1_110,p1_110,a2_110,p2_110,m,combData,simPdf,sample);
  plot_results("125",mean_125,width_125,a1_125,p1_125,a2_125,p2_125,m,combData,simPdf,sample);
  plot_results("130",mean_130,width_130,a1_130,p1_130,a2_130,p2_130,m,combData,simPdf,sample);
  plot_results("140",mean_140,width_140,a1_140,p1_140,a2_140,p2_140,m,combData,simPdf,sample);


  double mu_sim[4] = {mean_110.getVal(),mean_125.getVal(),mean_130.getVal(),mean_140.getVal()};
  double width_sim[4] = {width_110.getVal(),width_125.getVal(),width_130.getVal(),width_140.getVal()};
  double a1_sim[4] = {a1_110.getVal(),a1_125.getVal(),a1_130.getVal(),a1_140.getVal()};
  double p1_sim[4] = {p1_110.getVal(),p1_125.getVal(),p1_130.getVal(),p1_140.getVal()};
  double a2_sim[4] = {a2_110.getVal(),a2_125.getVal(),a2_130.getVal(),a2_140.getVal()};
  double p2_sim[4] = {p2_110.getVal(),p2_125.getVal(),p2_130.getVal(),p2_140.getVal()};

  plot_param(mu,mu_sim,"mu");
  plot_param(sigma,width_sim,"width");
  plot_param(a1,a1_sim,"a1");
  plot_param(p1,p1_sim,"p1");
  plot_param(a2,a2_sim,"a2");
  plot_param(p2,p2_sim,"p2");

  return;
}

void plot_results(string name, RooFormulaVar mean, RooFormulaVar width, RooRealVar a1, RooRealVar p1, RooRealVar a2, RooRealVar p2,
                  RooRealVar m, RooDataSet combData, RooSimultaneous simPdf, RooCategory sample){
  // plotting
  // adding text
  TPaveText *text_cut  = new TPaveText(0.05,0.88,0.77,0.98,"NDC");
  text_cut->AddText("Ph selection cutFlow > 13 and 125 GeV dataset removed");
  text_cut->SetBorderSize(0);
  text_cut->SetFillStyle(0);
  text_cut->SetTextAlign(12);
  text_cut->SetTextFont(42);
  text_cut->SetTextSize(0.02);

  string intro[6] = {"#mu = ","#sigma = ","a_{1} = ","p_{1} = ","a_{2} = ","p_{2} = "};
  string values[6] = {to_string(mean.getVal()),to_string(width.getVal()),to_string(a1.getVal()),to_string(p1.getVal()),to_string(a2.getVal()),to_string(p2.getVal())};
  string errors[6] = {" "," ",to_string(a1.getAsymErrorHi()),to_string(p1.getAsymErrorHi()),to_string(a2.getAsymErrorHi()),to_string(p2.getAsymErrorHi())};
  string pm = " +- ";

  TPaveText *text_fit  = new TPaveText(0.9, 0.65, 0.7, 0.9,"NDC");
  text_fit->SetFillStyle(0);
  text_fit->SetTextAlign(12);
  text_fit->SetTextFont(42);
  text_fit->SetTextSize(0.025);
  text_fit->AddText("DCB fit:");
  string line;
  for(int j=0; j<6; j++){
    if (j<2){
      line = intro[j]+values[j];
    }else{
      line = intro[j]+values[j]+pm+errors[j];
    }
    text_fit->AddText(line.c_str());
  }

  TCanvas* c1 = new TCanvas("c1","c1",600,400) ;
  c1->cd();
  const char *title_init = "m_{yy} distribution + DCB global fit [";
  const char *title_end = " GeV] with #mu(m_{h}) and #sigma(m_{h})";
  string title_plot = title_init + name + title_end;
  RooPlot* frame = m.frame(Title(title_plot.c_str()));
  const char * name_comb = "sample==sample::myy_";
  string path_comb = name_comb + name;
  combData.plotOn(frame,Name("Data"),Cut(path_comb.c_str())) ;
  const char * name_sim = "myy_";
  string path_sim = name_sim + name;
  //simPdf.plotOn(frame,Name("Fit"),Slice(sample,path_sim.c_str()),ProjWData(sample,combData)) ;
  simPdf.plotOn(frame,Name("Fit"),ProjWData(sample,combData)) ;
  // settings the axis names
  frame->SetXTitle("m_yy [GeV]");
  frame->SetYTitle("Events");

  // legend
  TLegend leg(0.15, 0.8, 0.48, 0.9);
  leg.AddEntry(frame->findObject("Data"), "Data Hist", "lep");
  leg.AddEntry(frame->findObject("Fit"), "DCB Fit", "L");
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.025);

  // drawing
  frame->Draw();
  text_fit->Draw();
  text_cut->Draw();
  leg.DrawClone();


  // saving plot
  const char *save_init = "Plots/myy_";
  const char *save_end = "GeV_DCBfit_cutFlow_global_125removed.pdf";
  string path_plot = save_init + name + save_end;
  c1->SaveAs(path_plot.c_str());

  return;
}

void plot_param(double P_sing[4], double P_all[4], string name){

  TCanvas* c3 = new TCanvas("c3","c3",600, 400);
  double masses[4] = {110,125,130,140};
  TGraph* gr1 = new TGraph(4,masses,P_sing);
  gr1->SetName("Single_fit");
  gr1->SetTitle("Single_fit");
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerColor(2);
  gr1->SetDrawOption("APL");
  gr1->SetLineColor(2);
  gr1->SetLineWidth(4);
  gr1->SetFillStyle(0);

  TGraph* gr2 = new TGraph(4,masses,P_all);
  gr2->SetName("Sim_fit");
  gr2->SetTitle("Sim_fit");
  gr2->SetMarkerStyle(22);
  gr2->SetMarkerColor(3);
  gr2->SetDrawOption("APL");
  gr2->SetLineColor(3);
  gr2->SetLineWidth(4);
  gr2->SetFillStyle(0);

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1);
  mg->Add(gr2);
  const char *title = "Single and global fit comp [no 125 GeV]: ";
  string plot_title = title + name;
  mg->SetTitle(plot_title.c_str());
  mg->GetXaxis()->SetTitle("m_{yy} [GeV]");
  mg->GetYaxis()->SetTitle("Params");
  mg->Draw("ALP");

  c3->BuildLegend();

  const char *repo = "Plots/";
  const char *end = "_param_125removed.pdf";
  string path = repo + name + end;

  c3->SaveAs(path.c_str());
}

RooDataSet *read_data(int MASS, RooRealVar& m){
  // usefull variables
  string name = to_string(MASS);
  int intervall = 30;
  // read variables
  float m_yy, weight;
  int cutFlow;
  // defining objetcs + datas
  const char *repo = "../Data/PowhegPy8_NNLOPS_ggH";
  const char *end = ".root";
  string path = repo + name + end;
  TFile *input = new TFile(path.c_str(),"read"); // reading the 140 GeV file;
  // datas
  TTree *datatree = (TTree*)input->Get("myCatNtuple;1");
  datatree->SetBranchAddress("m_yy", &m_yy);
  datatree->SetBranchAddress("EventWeight", &weight);
  datatree->SetBranchAddress("cutFlow", &cutFlow);

  // RooRealVar m("m_yy","mass",MASS-intervall,MASS+intervall);
  RooRealVar w("weight","weight",-1000,1000);

  RooDataSet * dataset = new RooDataSet("dataset","dataset",RooArgSet(m,w),WeightVar(w));
  for(int i=0;i<datatree->GetEntries();i++){
    datatree->GetEvent(i);
    m.setVal(m_yy*pow(10,-3));
    // photons selection with cutFlow > 13
    if(cutFlow > 13){
      dataset->add(RooArgSet(m),weight*128965.16);
    }
  }
  return dataset;
}


vector<double> param_lin(double *M, double *M_err, double *par, double *par_err, string name){
    vector<double> a_b;

    // preparing the fit
    TF1 * Lfit = new TF1("Linear fit","[0]+x*[1]",100,150);
    Lfit->SetLineColor(kRed);
    Lfit->SetLineStyle(2);
    Lfit->SetParNames("A","B");

    // fit
    TCanvas* c1 = new TCanvas("c1","mu fit",200,10,700,500);
    TGraphErrors* gr = new TGraphErrors(4,M,par,M_err,par_err);
    const char *dist = " distribution with a linear fit";
    string title = name + dist;
    gr->SetTitle(title.c_str());
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->GetXaxis()->SetTitle("m_yy [GeV]");
    const char *ytit = " fit";
    string ytitle = name + ytit;
    gr->GetYaxis()->SetTitle(ytitle.c_str());
    gr->Draw("ALP");
    // fitting
    gr->Fit(Lfit);
    cout << Lfit->GetChisquare() << endl;
    // TLegend
    TLegend* leg = new TLegend(0.15, 0.7, 0.35, 0.8,"Legend");
    leg->AddEntry(gr,"Exp data");
    leg->AddEntry(Lfit,"Fit");
    leg->Draw("SAME");
    // fit results
    TPaveText *text = new TPaveText(0.65, 0.2, 0.85, 0.3,"brNDC");
    string line = "y = "+to_string(Lfit->GetParameter(0))+" + "+to_string(Lfit->GetParameter(1))+"*x";
    text->AddText(line.c_str());
    string line_chi = "ChiQuadro = "+to_string(Lfit->GetChisquare());
    text->AddText(line_chi.c_str());
    text->Draw();
    // saving the plot
    const char *repo = "Plots/";
    const char *end = "_fit.pdf";
    string path = repo + name + end;

    c1->SaveAs(path.c_str());

    // saving linear fit parameters
    a_b.push_back(Lfit->GetParameter(0));
    a_b.push_back(Lfit->GetParameter(1));

    return a_b;
}


#include "TMath.h"
#include "TCanvas.h"
//root -l K.C /////to access IN x2g0 linux
void K(){
  gStyle ->SetOptStat(0);
  TFile file1("KL_flux_p_Pythia_24m.root");//Miks kaon cod//4.02.22 dat

  TFile file2("Cross_hKPXi.root");// co for cross ction
  TH1F *hKPXi = (TH1F*)file2.Get("hKPXi");// tget my cross sction
  hKPXi->GetYaxis()->SetTitle("Counts  [mb]");
  hKPXi->GetXaxis()->SetTitle("p [MeV/c]");


  TH1F *h24p = (TH1F*)file1.Get("h24p");
  // K_L n-> K^+ Xi^-; Xi^- p -> Lambda Lambda
  // Xi^- --> Lambda0 + pi^-
  // Detect K^+ & Lambda0 + pi^-


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Creating histograms here
  TFile fileOutput1("Kaon_Lambda_Output_test.root","recreate");

  TF1* decay_path  = new TF1("decay_path", "TMath::Exp(-x/([0]*[1]*30*0.1639))", 0, 41);
  decay_path->GetYaxis()->SetTitle("P(L)");
  decay_path->GetXaxis()->SetTitle("Length in Chamber [cm]");

  // Angular dependence histograms
  TH2D* h_Kaon_Plus=new TH2D("h_Kaon_Plus","Angular dependence #Theta of K^{+};p [GeV/c]; #Theta [degrees];",200,0,10,200,0,180);
  TH2D* h_Cascade_Minus=new TH2D("h_Cascade_Minus","Angular dependence #Theta of #Xi^{-};p [GeV/c]; #Theta [degrees];",200,0,11,200,0,180);
  TH2D* h_Lamba_Cascade1=new TH2D("h_Lamba_Cascade1","Angular dependence #Theta of #Lambda^{0};p [GeV/c]; #Theta [degrees];",200,0,11,200,0,180);
  //TH2D* h_Lamba_Cascade2=new TH2D("h_Lamba_Cascade2","Angular dependence #Theta of #Lambda^{0};p [GeV/c]; #Theta [degrees];",200,0,11,200,0,180);
  TH2D* h_decay_Proton1=new TH2D("h_decay_Proton1","Angular dependence #Theta of first Proton on #Lambda^{0};p [GeV/c]; #Theta [degrees];",200,0,10,200,0,180);
  TH2D* h_decay_Proton2=new TH2D("h_decay_Proton2","Angular dependence #Theta of Second Proton on #Lambda^{0};p [GeV/c]; #Theta [degrees];",200,0,10,200,0,180);
  TH2D* h_Pi_Minus_Cascade1=new TH2D("h_Pi_Minus_Cascade1","Angular dependence #Theta of first #pi^{-} on #Lambda^{0};p [GeV/c]; #Theta [degrees];",200,0,10,200,0,180);
  TH2D* h_Pi_Minus_Cascade2=new TH2D("h_Pi_Minus_Cascade2","Angular dependence #Theta of Second #pi^{-} on #Lambda^{0};p [GeV/c]; #Theta [degrees];",200,0,10,200,0,180);

  TH1F* h_K_Long=new TH1F("h_K_Long","Momentum of K_{L}; p [Gev/c]; Counts",200,0,5);///Histo for momentum

  TH1F* h_Flux_Cross=new TH1F("h_Flux_Cross","Flux X cross-section; p [Gev/c]; Normalised Yield",12000,0,12);

  TH1F* h_Flux_Cross3=new TH1F("h_Flux_Cross3","Flux X cross-section; p [Gev/c]; Normalised Yield",12000,0,12);

  TH1F* h_Flux_Cross2=new TH1F("h_Flux_Cross2","Flux X cross-section; p [Gev/c]; Normalised Yield",12000,0,12);
  for (size_t i = 0; i < 12001; i++) {
    h_Flux_Cross2->SetBinContent(i,hKPXi->GetBinContent(i));
    h_Flux_Cross3->SetBinContent(i,hKPXi->GetBinContent(i));

  }

  TH1F* h_Flux_Cross2_scal = new TH1F("h_Flux_Cross2_scal","Scaled K_{L} with mean pl ; p [Gev/c]; Normalised Yield",12000,0,12);

  TH1F* h_Flux_Cross2_scal2 = new TH1F("h_Flux_Cross2_scal2","Scaled K_{L} with 40cm l ; p [Gev/c]; Normalised Yield",12000,0,12);




  h_Flux_Cross2_scal2->Multiply(h_Flux_Cross2, h24p);
  h_Flux_Cross2_scal->Multiply(h_Flux_Cross2, h24p);

  h_Flux_Cross2->Multiply(h_Flux_Cross2, h24p);



//ND H24P!!!!!!!!
  TH2D* h_Kaon_Plus_Phi=new TH2D("h_Kaon_Plus_Phi","Angular dependence #Phi of K^{+};p [GeV/c]; #Phi [degrees];",200,0,11,200,-180,180);
  TH2D* h_Cascade_Minus_Phi=new TH2D("h_Cascade_Minus_Phi","Angular dependence #Phi of #Xi^{-};p [GeV/c]; #Phi [degrees];",200,0,11,200,-180,180);
  TH2D* h_Lamba_Cascade1_Phi=new TH2D("h_Lamba_Cascade1_Phi","Angular dependence #Phi of #Lambda^{0};p [GeV/c]; #Phi [degrees];",200,0,11,200,-180,180);
  //TH2D* h_Lamba_Cascade2_Phi=new TH2D("h_Lamba_Cascade2_Phi","Angular dependence #Phi of #Lambda^{0;p [GeV/c]; #Phi [degrees];}",200,0,11,200,-180,180);
  TH2D* h_Pi_Minus_Cascade1_Phi=new TH2D("h_Pi_Minus_Cascade1_Phi","Angular dependence #Phi of #pi^{-} #Lambda^{0};p [GeV/c]; #Phi [degrees];",200,0,11,200,-180,180);
  TH2D* h_Pi_Minus_Cascade2_Phi=new TH2D("h_Pi_Minus_Cascade2_Phi","Angular dependence #Phi of #pi^{-} #Lambda^{0};p [GeV/c]; #Phi [degrees];",200,0,11,200,-180,180);

  TH2D* h_xy=new TH2D("h_xy","Generation of X-Y  Coordinates;x [cm];y [cm]",200,-0.1,0.1,200,-0.1,0.1);

  // Invariant mass histograms
  //TH1F* h_inv_Cascade=new TH1F("h_inv_Lambda","Invariant mass of #Cascade",200,0,2);
  //TH1F* h_inv_Cascade_FD=new TH1F("h_inv_Cascade_FD","Invariant mass of #Cascade after angle cut",200,0,2);



  TH2D* h_t_path = new TH2D("h_t_path", "transverse path length vs Phi; #Phi [rad]; Path Length [m]", 1000, 0, TMath::Pi(), 1000, 0, 0.04);
  TH2D* h_l_p = new TH2D("h_l_p", "longitudinal path length vs Theta; #Theta [rad]; Path Length [m]", 1000, 0, 2, 1000, 0, 0.4);
  TH2D* h_L = new TH2D("h_L", "Grand path length vs Theta; #Theta [rad]; Path Length [m]", 1000, 0, 2/*TMath::Pi()*/, 1000, 0, 0.4);
  TH2D* h_L_D = new TH2D("h_L_D", "Decay path length vs Theta; #Theta [rad]; Path Length [m]", 1000, 0, 2/*TMath::Pi()*/, 1000, 0, 0.4);
//Acceptance AS OFF 6 MAY:
  TH1F* h_KL_all=new TH1F("h_KL_all","Acceptance of total K_{L} Generated; p [Gev/c]; Normalised Yield [particles/s/GeV]",12000,0,12);
  TH1F* h_KL_accp=new TH1F("h_KL_accp","Acceptance of detected K_{L} ; p [Gev/c]; Normalised Yield [particles/s/GeV]",12000,0,12);
  TH1F* h_KL_Accpratio=new TH1F("h_Accpratio","Ratio of detected to Generated K_{L} ; p [Gev/c]; Normalised Yield",12000,0,12);

  TH1F* h_KL_Accpratio1=new TH1F("h_Accpratio1","Ratio of detected to Generated K_{L} ; p [Gev/c]; Acceptance",12000,0,12);


  TH2D* h_Kaon_Plus_AC = new TH2D("h_Kaon_Plus_AC", "Acceptance of K^{+};p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0, 180);
  TH2D* h_Pi_Minus_Cascade1_accp = new TH2D("h_Pi_Minus_Cascade1_accp", "Acceptance of First #pi^{-};p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0,180);
  TH2D* h_Pi_Minus_Cascade2_accp = new TH2D("h_Pi_Minus_Cascade2_accp", "Acceptance of Second #pi^{-};p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0, 180);
  TH2D* h_decay_Proton1_accp = new TH2D("h_decay_Proton1_accp", "Acceptance of First Proton;p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0, 180);
  TH2D* h_decay_Proton2_accp = new TH2D("h_decay_Proton2_accp", "Acceptance of Second Proton;p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0, 180);

  TH2D* h_Kaon_Plus_Con = new TH2D("h_Kaon_Plus_Con", "Constraned K^{+};p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0,180);
  TH2D* h_Pi_Minus_Cascade1_Con = new TH2D("h_Pi_Minus_Cascade1_Con", "Constraned First #pi^{-};p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0,180);
  TH2D* h_Pi_Minus_Cascade2_Con = new TH2D("h_Pi_Minus_Cascade2_Con", "Constraned Second #pi^{-};p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0,180);
  TH2D* h_decay_Proton1_Con = new TH2D("h_decay_Proton1_Con", "Constraned First Proton;p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0,180);
  TH2D* h_decay_Proton2_Con = new TH2D("h_decay_Proton2_Con", "Constraned Second Proton;p [GeV/c]; #Theta [degrees];", 200, 0, 10, 200, 0,180);
//Scale graphs

  //h_Flux_Cross2_scal2=(TH1F*)h_Flux_Cross2->Clone("h_Flux_Cross2_scal2");//k long
  //h_Flux_Cross2_scal=(TH1F*)h_Flux_Cross2->Clone("h_Flux_Cross2_scal");//meanpl
////////////////////////////
  TH1F* h_KL_events = new TH1F("h_KL_events","KL events x accpratio; p [Gev/c]; Normalised Yield",12000,0,12);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Creating particles and terms

  // How many events to simulate and percentage completed
  Int_t nevents=100000;
  Int_t Percentage=nevents/100;

  // Creating TLorentzVectors of particles
  TLorentzVector Beam, Target,Targetp; // Beam and target
  TLorentzVector *Kaon_Plus, *Cascade_Minus; // First vertex particles
  ///add Cascade with proton ?
  TLorentzVector *Lambda_Cascade1,*Lambda_Cascade2, *Pi_Minus_Cascade1, *decay_Proton1; // Second vertex particles
  TLorentzVector *Lambda_Cascade3,*Lambda_Cascade4, *Pi_Minus_Cascade2, *decay_Proton2; //third vertex particles
  // q is momentum transfer from beam to target, used to calculate q^2 weight
  TLorentzVector q;

  // Making Weights
  Double_t Phasespace_Weight_1;
  Double_t Phasespace_Weight_2;

  Double_t qWeight;

  // Setting TLorentzVectors for beam and target in GeV (Px,Py,Pz,M)
  //Mass of k long is 0.493677
  Target.SetXYZM(0,0,0,0.93957); //Mass of Neutron 0.93957
  Targetp.SetXYZM(0,0,0,0.938);

  // Defining initial vertex and masses of particles
  TLorentzVector V1 = Beam + Target;
  Double_t Masses_1[2] = {0.49368,1.321}; // kaon+ and Xi mass
  Double_t Masses_2[2] = {1.11568,1.11568}; // lambda, pi^- (Lambda decay into proton and pi-)// xi mass ,and proton ,then lambda
  Double_t Masses_3[2] = {0.93827,0.139570}; // lambda dcay into proton and pion for 1st Lambda
  Double_t Masses_4[2] = {0.93827,0.1395070}; // lambda dcay inot proton and pion for 2nd Lambda

  TLorentzVector Inv_Cascade; // 4-vector of Cascade from lambda + pion

  // Creating decay vertices
  TGenPhaseSpace Vertex_1, Vertex_2, Vertex_3, Vertex_4, Vertex_5;

  Double_t BeamRandom;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over events
TRandom3* myz=new TRandom3(0);
TRandom3* gnxy=new TRandom3(0);
Double_t z;
Double_t x;
Double_t y;
Double_t T_L = 0.4;
Double_t T_R = 0.03;
Double_t Xi_Phi;
Double_t Xi_Theta;
Double_t Xi_Phi_V;
Double_t alpha;
Double_t t_path;
Double_t l_p;
Double_t z_exit;
Double_t angle_counter;
Double_t L;
Double_t L_Choice;
Double_t vertex_d;
Double_t Xi_Beta;
Double_t Xi_Gamma;
Double_t GDecay;
  // Looping over simulated events
  for (Long64_t i=0;i<nevents;i++){

    // Counter, shows percentage of simulated events completed
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }
    z=myz->Rndm()*0.4;
    //std::cout<<z<<endl;
    Double_t radius=0.04;
    while (radius>0.03) {
      x=gnxy->Gaus(0,0.006);//sigma as k long
      y=gnxy->Gaus(0,0.006);
      radius=TMath::Sqrt((x*x)+(y*y));
    }
    h_xy->Fill(x,y);

    BeamRandom= h24p->GetRandom();
    Beam.SetXYZM(0,0,BeamRandom,0.493677);
    V1=Beam + Target;
    // SetDecay(total energy, no. particles, mass array)
    // Setting the first decay
    if(!Vertex_1.SetDecay(V1,2,Masses_1))continue;

    h_KL_all->Fill(BeamRandom);

    // Generating event and defining the phasespace weight for this decay
    Phasespace_Weight_1=Vertex_1.Generate();

    // Assigning the decay particles from the array above
    Kaon_Plus=Vertex_1.GetDecay(0);
    Cascade_Minus=Vertex_1.GetDecay(1);

    //Lambda=Vertex_1.GetDecay(1); ///////////////These connect to both and assign masses. so 0 i first part of mass array above
    //Kaon_Plus=Vertex_1.GetDecay(2);

    // Defining q^2 weight
    //q=Beam-(TLorentzVector)*Kaon_Plus; // Subtracting e^- out from e^- in
    //qWeight=1/(q.Rho()*q.Rho()); // .Rho() gives momentum of the 4-vector


    // Getting total energy for second decay from Lambda
    TLorentzVector V2 = (TLorentzVector)*Cascade_Minus+Targetp;

    // Setting the second decay
    Vertex_2.SetDecay(V2,2,Masses_2);
    // Defining the phasespace for this decay
    Double_t Phasespace_Weight_2=Vertex_2.Generate();
    Lambda_Cascade1=Vertex_2.GetDecay(0);
    Lambda_Cascade2=Vertex_2.GetDecay(1);//NAM Lambda
//CUT OFF
    TLorentzVector V3 = (TLorentzVector)*Lambda_Cascade1;
    Vertex_3.SetDecay(V3, 2, Masses_3);

    Double_t Phasespace_Weight_3=Vertex_3.Generate();
    decay_Proton1 = Vertex_3.GetDecay(0);
    Pi_Minus_Cascade1= Vertex_3.GetDecay(1);




    TLorentzVector V4 = (TLorentzVector)*Lambda_Cascade2;
    Vertex_4.SetDecay(V4, 2, Masses_4);

    Double_t Phasespace_Weight_4=Vertex_4.Generate();
    decay_Proton2 = Vertex_4.GetDecay(0);
    Pi_Minus_Cascade2 = Vertex_4.GetDecay(1);

    Xi_Theta=Cascade_Minus->Theta();
    Xi_Phi=Cascade_Minus->Phi();
    Xi_Phi_V=TMath::ATan2(y,x);
    vertex_d=TMath::Sqrt(x*x+y*y);
    alpha=TMath::Pi()-Xi_Phi_V+Xi_Phi;
    t_path = vertex_d * TMath::Cos(alpha) + TMath::Sqrt(vertex_d*vertex_d * (TMath::Cos(alpha) * TMath::Cos(alpha)) + T_R*T_R - vertex_d*vertex_d);
    l_p=t_path/TMath::Tan(Xi_Theta);//long path z comp
    z_exit=l_p+z;

    if (z_exit > 0.4) {
        l_p = T_L - z;
        L = l_p / TMath::Cos(Xi_Theta);
        t_path = l_p * TMath::Tan(Xi_Theta);
    }
    else  {
        L = t_path / TMath::Sin(Xi_Theta);
    }
  /*
    else if (z_exit < 0) {
        l_p = z;
        L = l_p / TMath::Cos(Xi_Theta);
        t_path = l_p * TMath::Tan(Xi_Theta);
    }
*/


    Xi_Beta=Cascade_Minus->Beta();
    Xi_Gamma=Cascade_Minus->Gamma();
    decay_path->SetParameter(0,Xi_Beta);
    decay_path->SetParameter(1,Xi_Gamma);
    GDecay=decay_path->GetRandom()/100;
    if (GDecay < L) {
        L_Choice = GDecay;
    }
    else if (GDecay > L) {
        L_Choice = L;
    }

h_L_D->Fill(Xi_Theta, L_Choice);
h_L->Fill(Xi_Theta, L);
h_l_p->Fill(Xi_Theta, l_p);
h_t_path->Fill(Xi_Phi, t_path);

//fill Theta
h_Kaon_Plus->Fill(Kaon_Plus->Rho(),Kaon_Plus->Theta()*TMath::RadToDeg(),Phasespace_Weight_1);
h_Cascade_Minus->Fill(Cascade_Minus->Rho(),Cascade_Minus->Theta()*TMath::RadToDeg(),Phasespace_Weight_1);
h_Lamba_Cascade1->Fill(Lambda_Cascade1->Rho(),Lambda_Cascade1->Theta()*TMath::RadToDeg(),Phasespace_Weight_1*Phasespace_Weight_2);
h_Pi_Minus_Cascade1->Fill(Pi_Minus_Cascade1->Rho(),Pi_Minus_Cascade1->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_3);
h_Pi_Minus_Cascade2->Fill(Pi_Minus_Cascade2->Rho(),Pi_Minus_Cascade2->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_4);
h_decay_Proton1->Fill(decay_Proton1->Rho(),decay_Proton1->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_3);
h_decay_Proton2->Fill(decay_Proton2->Rho(),decay_Proton2->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_4);




h_Kaon_Plus_Phi->Fill(Kaon_Plus->Rho(),Kaon_Plus->Phi()*TMath::RadToDeg(),Phasespace_Weight_1);
h_Cascade_Minus_Phi->Fill(Cascade_Minus->Rho(),Cascade_Minus->Phi()*TMath::RadToDeg(),Phasespace_Weight_1);
h_Lamba_Cascade1_Phi->Fill(Lambda_Cascade1->Rho(),Lambda_Cascade1->Phi()*TMath::RadToDeg(),Phasespace_Weight_1*Phasespace_Weight_2);
h_Pi_Minus_Cascade1_Phi->Fill(Pi_Minus_Cascade1->Rho(),Pi_Minus_Cascade1->Phi()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2);




//add w
if (Kaon_Plus->Rho() > 0.2 && Pi_Minus_Cascade1->Rho() > 0.2 && Pi_Minus_Cascade2->Rho() > 0.2 && decay_Proton1->Rho() > 0.2 && decay_Proton2->Rho() > 0.2) {
        if ((( Pi_Minus_Cascade1->Theta() > 3 * TMath::DegToRad() && Pi_Minus_Cascade1-> Theta() < 15 * TMath::DegToRad()) || (Pi_Minus_Cascade1-> Theta() > 20 * TMath::DegToRad() && Pi_Minus_Cascade1->Theta() < 165 * TMath::DegToRad()))){
            if (((Pi_Minus_Cascade2->Theta() > 3 * TMath::DegToRad() && Pi_Minus_Cascade2-> Theta() < 15 * TMath::DegToRad()) || (Pi_Minus_Cascade2-> Theta() > 20 * TMath::DegToRad() && Pi_Minus_Cascade2->Theta() < 165 * TMath::DegToRad()))) {
                if (((decay_Proton1->Theta() > 3 * TMath::DegToRad() && decay_Proton1->Theta() < 15 * TMath::DegToRad()) || (decay_Proton1->Theta() > 20 * TMath::DegToRad() && decay_Proton1->Theta() < 165 * TMath::DegToRad()))) {
                    if (((decay_Proton2->Theta() > 3 * TMath::DegToRad() && decay_Proton2->Theta() < 15 * TMath::DegToRad()) || (decay_Proton2->Theta() > 20 * TMath::DegToRad() && decay_Proton2->Theta() < 165 * TMath::DegToRad()))) {
                      if (((Kaon_Plus->Theta() > 3 * TMath::DegToRad() && Kaon_Plus->Theta() < 15 * TMath::DegToRad()) || (Kaon_Plus->Theta() > 20 * TMath::DegToRad() && Kaon_Plus->Theta() < 165 * TMath::DegToRad()))) {
                        h_Kaon_Plus_Con->Fill(Kaon_Plus->Rho(), Kaon_Plus->Theta()*TMath::RadToDeg(), Phasespace_Weight_1);
                        h_Pi_Minus_Cascade1_Con->Fill(Pi_Minus_Cascade1->Rho(),Pi_Minus_Cascade1->Theta()*TMath::RadToDeg(), Phasespace_Weight_3);
                        h_Pi_Minus_Cascade2_Con->Fill(Pi_Minus_Cascade2->Rho(), Pi_Minus_Cascade2->Theta()*TMath::RadToDeg(), Phasespace_Weight_4);
                        h_decay_Proton1_Con->Fill(decay_Proton1->Rho(), decay_Proton1->Theta()*TMath::RadToDeg(), Phasespace_Weight_3);
                        h_decay_Proton2_Con->Fill(decay_Proton2->Rho(), decay_Proton2->Theta()*TMath::RadToDeg(), Phasespace_Weight_4);
                        h_KL_accp->Fill(BeamRandom, Phasespace_Weight_2 * Phasespace_Weight_3);

                    }
                }
            }
        }
    }
  }


//kn->k^+xi->LL->p pi........kp->Lpi+,

//CUT OFF
    // Adding 4-vector of proton and pion, use to find invariant mass of lambda with .M()
    //Inv_Cascade = (TLorentzVector)*Lambda_Cascade1 + (TLorentzVector)*Pi_Minus_Cascade;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Filling histograms (x, weights applied)
    //h_inv_Cascade->Fill(Inv_Cascade.M(),Phasespace_Weight_1);


    // Setting Forward Detector acceptance cut for kaon_plus
    //if(Kaon_Plus->Theta()*TMath::DegToRad()>5 && Kaon_Plus->Theta()*TMath::DegToRad()<35){

      // Filling angular distribution histograms (x,y,weights applied)



      // Filling invariant mass histogram
      //h_inv_Cascade_FD->Fill(Inv_Cascade.M(), Phasespace_Weight_1*Phasespace_Weight_2*qWeight);

    //}



}
//bin all 3 sam way
h_KL_Accpratio->Divide(h_KL_accp,h_KL_all);//dc,gn
h_KL_Accpratio1->Rebin(50);
h_KL_all->Rebin(50);
h_KL_accp->Rebin(50);

h_KL_Accpratio1->SetMaximum(1);
h_KL_Accpratio1->SetMinimum(0);

h_KL_Accpratio1->Divide(h_KL_accp,h_KL_all);//dc,gn
h_KL_Accpratio1->Draw("hist");

h_KL_Accpratio->Draw();

h_Kaon_Plus_AC->Divide(h_Kaon_Plus_Con,h_Kaon_Plus);
h_Pi_Minus_Cascade1_accp->Divide(h_Pi_Minus_Cascade1_Con,h_Pi_Minus_Cascade1);//cons w/cuts div / w/o cut
h_Pi_Minus_Cascade2_accp->Divide(h_Pi_Minus_Cascade2_Con,h_Pi_Minus_Cascade2);
h_decay_Proton1_accp->Divide(h_decay_Proton1_Con,h_decay_Proton1);
h_decay_Proton2_accp->Divide(h_decay_Proton2_Con,h_decay_Proton2);


decay_path->Draw("hist");

h_Flux_Cross2->Draw();
Double_t meanpl;
meanpl=h_L_D->GetMean(2);
//h_Flux_Cross2_scal->Scale(meanpl*100*(0.602/2.0141)*0.169/1000);//casca

h_Flux_Cross2_scal->Scale(40*(0.602/2.0141)*0.169/1000000);//casca Flux
h_Flux_Cross2_scal2->Scale(40*(0.602/2.0141)*0.169/1000000);//k long

h_Flux_Cross2_scal->Scale(meanpl*100*30*(0.602/2.0141)*0.169/1000);//casca///ask maya 30


h_Flux_Cross2_scal->Draw("hist");
h_Flux_Cross2_scal2->Draw("hist");

h_Flux_Cross2_scal2->GetYaxis()->SetTitle("Normalised Yield  [particles/s/GeV]");
h_Flux_Cross2_scal->GetYaxis()->SetTitle("Normalised Yield  [particles/s/GeV]");

cout<<"Total #Xi^{-} Yield :  "<<h_Flux_Cross2_scal2->Integral()*100*24*3600<<endl;
cout<<"Total Scattered Yield : "<<h_Flux_Cross2_scal->Integral()*100*24*3600<<endl;
///Scale by 100 24 3600 after 373
//h_Flux_Cross2_scal->Rebin(50);//Scale 1.0/50
//h_Flux_Cross2_scal2->Rebin(50);
//h_KL_Accpratio->Rebin(50);
//h_KL_events->Rebin(50);

h_KL_events->Multiply(h_Flux_Cross2_scal,h_KL_Accpratio);//include in dis as
cout<<"Number of Measured events :  "<<h_KL_events->Integral()*100*24*3600<<endl;

//bin n..scal 1/n.....ti events[particles/s/MeV]...min0
/*
for (int k = 1; k < 12001; k++) {
  if (h_K_Long1->GetBinCenter(k)>1) {
  h_K_Long1->SetBinContent(k, 3 / (h_K_Long1->GetBinCenter(k)* h_K_Long1->GetBinCenter(k)));
  }
  else {
    h_K_Long1->SetBinContent(k,3);
  }
  h_K_Long2->SetBinContent(k, 10);


}
*/


//rho ld2=0.180....

cout<<"Mean Path Length :   "<<meanpl<<endl;
decay_path->Draw();


// Find acceptance by dividing histogram after angle cut by original histogram
// -> Divide (A,B) gives A/B
//h_acceptance_Cascade_FD->Divide(h_inv_Cascade_FD,h_inv_Cascade);

//Writes histograms to file
//fileOutput1.Write();
/*
//Draw histograms
TCanvas* c2 = new TCanvas("c2","stacked hists ",1000,1000);
h_Cascade_Minus->Draw("colz");

TCanvas* c4 = new TCanvas("c4","",1000,1000);
h_Cascade_Minus_Phi->Draw("colz");

TCanvas* c5 = new TCanvas("c5","",1000,1000);
h_Pi_Minus_Cascade->Draw("colz");

TCanvas* c6 = new TCanvas("c6","",1000,1000);
h_Pi_Minus_Cascade_Phi->Draw("colz");
*/


//colz is drawing option (colours the z axis)

//TCanvas* c3 = new TCanvas("c3","stacked hists",800,800);
//c3->cd(1);
//h_Kaon_Plus->DrawClone("colz");

//Writes histograms to file
hKPXi->Write();
decay_path->Write();
h24p->Write();
gStyle->SetOptStat(0);
fileOutput1.Write();

}

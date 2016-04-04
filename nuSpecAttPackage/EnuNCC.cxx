#ifndef ENUNCC_CXX
#define ENUNCC_CXX

#include "EnuNCC.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EnuNCC::EnuNCC( TTree * fOutTree ){
    SetOutTree(fOutTree);
    InitOutTree();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::ImpXsecGraph ( TString XsecFileName , const int Npoints , bool DoPlot){

    double Ev[Npoints] , sigma_over_Ev[Npoints], sigma[Npoints];
    analysis.ReadGraphFromFile( XsecFileName, Npoints, Ev, sigma_over_Ev );
    TGraph * g = new TGraph( Npoints , Ev, sigma_over_Ev  );
    SetXsecGraph( g );

    
    for (int i = 0; i < Npoints; i++)
        sigma[i] = sigma_over_Ev[i] * Ev[i] ;
    TGraph * gE = new TGraph( Npoints , Ev , sigma  );
    SetXsecEGraph( gE );
    
    if(DoPlot) DrawXsecGraph();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::ImpMomentumDist(bool DoPlot){
    // implement momentum distributions: Fermi Gas (FG) and Correlated Fermi Gas (CFG) 
    
//    
//    Double_t        c0 = 4.16 ; // +/- 0.95..
//    Double_t       k_F = 250 ;
//    Double_t       E_F = k_F * k_F / (2. * 940.); // non relativistic
//    Double_t      E_F0 = 12.5 / ( ( pow(2.,2./3) - 1 ) * (3./5.));
//    Double_t      k_F0 = sqrt( 2. * 940. * E_F0 );
//    Double_t  rho2rho0 = pow( ( pow(2.,2./3) - 1 ) * (3./5.) * (E_F/12.5)   , 3./2.);
//    Double_t    lambda = 2.75 ;// +/- 0.25 high-momentum cutoff
//    Double_t       pi2 = pow(TMath::Pi(),2);
//    Double_t        A0 = (3 * pi2 / pow(k_F0,3)) * (1./rho2rho0) * (1 - (1 - (1./lambda)*pow(rho2rho0,1./3.)) * (c0 / pi2 ));
//    SRCk4Tail          = new TF1("tail",Form("%f*%f/(x**4)",c0,k_F),k_F,lambda*k_F0);
//
//    int Nbins = 100;
//    float bin_content , k;
//
//    // Fermi Gas momentum dist.
//    hFG = new TH1F("hFG" , "Fermi Gas momentum dist." , Nbins , 0 , lambda*k_F0 );
//    for (int bin = 0; bin < Nbins ; bin++){
//        k = hFG -> GetXaxis() -> GetBinCenter( bin );
//        bin_content = ( k < k_F ) ? 1 : 0;
//        hFG -> SetBinContent( bin , bin_content );
//    }
//
//    
//    // Correlated Fermi Gas momentum dist.
//    hCFG = new TH1F("hCFG" , "Correlated Fermi Gas momentum dist." , Nbins , 0 , lambda*k_F0 );
//    for (int bin = 0; bin < Nbins ; bin++){
//        k = hCFG -> GetXaxis() -> GetBinCenter( bin );
//        bin_content = ( k < k_F ) ? A0 : SRCk4Tail -> Eval(k);
//        hCFG -> SetBinContent( bin , bin_content );
//    }
//    
//    float Integral = 0 , Tail = 0;
//    for (int bin = 0; bin < Nbins ; bin++){
//        k = hCFG -> GetXaxis() -> GetBinCenter( bin );
//        Integral += k * k * hCFG -> GetBinContent( bin );
//        Tail     += ( k_F < k ) ? k * k * hCFG -> GetBinContent( bin ) : 0 ;
//    }
//    float FracTail = Tail / Integral;//analysis.IntegralH1D(hCFG , 250 , lambda*k_F0 ) / analysis.IntegralH1D(hCFG , 0 , lambda*k_F0 );

    hCFG = calculations.CFGMomentumDist(DoPlot,SRCk4Tail);
    
    if (DoPlot)
        DrawMomentumDist();
    

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::ImpEfluxHisto ( TH1F * h , bool DoPlot ){
    hEflux = new TH1F("hEflux" , "neutrino energy flux" , h->GetNbinsX() , 0 , h->GetXaxis()->GetXmax()/1000. );
    for (int bin = 0 ; bin < h->GetNbinsX(); bin++) {
        hEflux -> SetBinContent(bin , h->GetBinContent(bin));
    }
    if(DoPlot) DrawEfluxGraph();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::ImpEfluxGraph ( TString EfluxFileName , const int Npoints , bool DoPlot){
    
    double Ev[Npoints] , flux_Units[Npoints], flux[Npoints];
    analysis.ReadGraphFromFile( EfluxFileName, Npoints, Ev, flux_Units );
    TGraph * g_Units = new TGraph( Npoints , Ev, flux_Units  );
    SetEfluxGraph( g_Units );
    
    
    for (int i = 0; i < Npoints; i++)
        flux[i] = flux_Units[i] ;
    TGraph * gE = new TGraph( Npoints , Ev , flux  );
    SetEfluxEGraph( gE );
    hEflux = new TH1F("hEflux" , "neutrino energy flux" , 100 , 0 , gEfluxE->GetXaxis()->GetXmax() );
    for (int bin = 0; bin < 100; bin++)
        hEflux -> SetBinContent( bin , gE -> Eval( hEflux->GetXaxis()->GetBinCenter(bin) ) );
    
    if(DoPlot) DrawEfluxGraph();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::InitOutTree(){
    
    OutTree -> Branch("nu"          ,"TLorentzVector"       ,&nu);
    OutTree -> Branch("n"           ,"TLorentzVector"       ,&n);
    OutTree -> Branch("nu_INnRF"    ,"TLorentzVector"       ,&nu_INnRF);    // nu in neutron RestFrame
    OutTree -> Branch("W"           ,"TLorentzVector"       ,&W);
    OutTree -> Branch("mu"          ,"TLorentzVector"       ,&mu);
    OutTree -> Branch("p"           ,"TLorentzVector"       ,&p);
    OutTree -> Branch("prec"        ,"TLorentzVector"       ,&prec);
    OutTree -> Branch("Pcm"         ,"TVector3"             ,&Pcm);
    OutTree -> Branch("Ev_INnRF"    ,&Ev_INnRF              ,"Ev_INnRF/D"); // nu energy in neutron RestFrame
    OutTree -> Branch("XsecWeight"  ,&XsecWeight            ,"XsecWeight/D");    // cross section weight
    OutTree -> Branch("XsecLabFrame",&XsecLabFrame          ,"XsecLabFrame/D");
    OutTree -> Branch("nDirection"  ,&nDirection            ,"nDirection/C");    // cross section weight
    
    
    std::cout << "Initialized Output Tree EnuNCC on " << OutTree -> GetTitle() << std::endl;
    
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::RunInteractions ( TString NuclearModel , TString nuFlux , int Ninteractions , bool DoPrint ){
    for ( int i = 0 ; i < Ninteractions ; i++ ) {
        if (i%(Ninteractions/20)==0) plot.PrintPercentStr((float)i/Ninteractions);

        GenerateNeutrino( nuFlux );
        GenerateNeutron( NuclearModel );
        CalcRestFrameEv();
        GenerateRecoilProton( 0.0 , 0.14 );
        if(DoPrint) PrintDATA(i);
        if(n.P() > 0.25) OutTree -> Fill();
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::GenerateNeutrino( TString nuFlux ){
    if (nuFlux == "BNB flux" || nuFlux == "Ev from mcc6") {
        Ev = hEflux -> GetRandom();
    }
    else if (nuFlux == "monochromatic neutrino 300 MeV"){
        Ev = 0.3;
    }
    nu = TLorentzVector( 0 , 0 , Ev , Ev );
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::GenerateNeutron( TString NuclearModel ){

    
    if (NuclearModel == "nUniform") {
        Pn = rand.Uniform(0,1);
        rand.Sphere(Px,Py,Pz,Pn);
        n = TLorentzVector( Px , Py , Pz , sqrt( Pn*Pn + Mn*Mn ) );
    }

    
    else if (NuclearModel == "nBack") {
        Pn = rand.Uniform(0,1);
        rand.Sphere(Px,Py,Pz,Pn);
        n = TLorentzVector( Px , Py , -fabs(Pz) , sqrt( Pn*Pn + Mn*Mn ) );
    }
    
    else if (NuclearModel == "nForward") {
        Pn = rand.Uniform(0,1);
        rand.Sphere(Px,Py,Pz,Pn);
        n = TLorentzVector( Px , Py , fabs(Pz) , sqrt( Pn*Pn + Mn*Mn ) );
    }

    else if (NuclearModel == "CFG") {
        nInSRC = (rand.Uniform() < 0.20) ? true : false;
        if (!nInSRC){
            Pn = rand.Uniform(0,0.25);
        }
        else if (nInSRC) {
            Pn = SRCk4Tail -> GetRandom() / 1000;
        }
        rand.Sphere(Px,Py,Pz,Pn);
        n = TLorentzVector( Px , Py , Pz , sqrt( Pn*Pn + Mn*Mn ) );
    }
    
    else if (NuclearModel == "CFGnBack") {
        nInSRC = (rand.Uniform() < 0.20) ? true : false;
        if (!nInSRC)
            Pn = rand.Uniform(0,0.25);
        else if (nInSRC)
            Pn = SRCk4Tail -> GetRandom() / 1000;
        rand.Sphere(Px,Py,Pz,Pn);
        n = TLorentzVector( Px , Py , -fabs(Pz) , sqrt( Pn*Pn + Mn*Mn ) );
    }
    
    else if (NuclearModel == "CFGnForward") {
        nInSRC = (rand.Uniform() < 0.20) ? true : false;
        if (!nInSRC)
            Pn = rand.Uniform(0,0.25);
        else if (nInSRC)
            Pn = SRCk4Tail -> GetRandom() / 1000;
        rand.Sphere(Px,Py,Pz,Pn);
        n = TLorentzVector( Px , Py , fabs(Pz) , sqrt( Pn*Pn + Mn*Mn ) );
    }
    
    else if (NuclearModel == "nParallelAntiparallel") {
        int i = (int) (3 * rand.Uniform()) ;
        switch (i) {
            case 1:
                Pn = 0.4;
                break;
            case 2:
                Pn = -0.4;
                break;
            default:
                Pn = 0.0;
                break;
        }
        n = TLorentzVector( 0 , 0 , Pn , sqrt( Pn*Pn + Mn*Mn ) );
    }
    nDirection = (n.Pz() < 0) ? "backward" : "forward";
    
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::GenerateRecoilProton( float mean , float sigma ){
    Pcm = TVector3( rand.Gaus( mean , sigma ) , rand.Gaus( mean , sigma ) , rand.Gaus( mean , sigma ) );
    if (n.P() > 0.25) { // SRC pair above kF
        prec.SetVectM( Pcm - n.Vect()  , Mp );
    }
    else {
        prec.SetVectM( TVector3(0 , 0 , 0) , Mp );
    }
    precDirection = (n.Pz() < 0) ? "backward" : "forward";

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::CalcRestFrameEv(){
    nu_INnRF = nu;    // neutrino in neutron rest frame
    nu_INnRF.Boost( -n.BoostVector() );
    Ev_INnRF = nu_INnRF.E();
    XsecWeight = gXsecE -> Eval ( Ev_INnRF ); // spline interpolation between points
    XsecLabFrame = gXsecE -> Eval ( Ev );
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::PrintDATA(int entry){
    SHOW(entry);
    SHOWTLorentzVector(nu);
    SHOWTLorentzVector(n);
    SHOWTVector3(Pcm);
    SHOWTLorentzVector(prec);
    SHOWTLorentzVector(nu_INnRF);
    SHOW(Ev);
    SHOW(Ev_INnRF);
    SHOW((Ev-Ev_INnRF)/Ev);
    SHOW(XsecWeight);
    SHOW(XsecLabFrame);
    SHOW(XsecWeight/XsecLabFrame);
    PrintLine();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::DrawXsecGraph(){
    
    TCanvas * c = plot.CreateCanvas("Xsec","Divide",2,1);
    
    c -> cd(1);
    plot.SetFrame(gXsec,"#nu - n CCQE Xsec, NUANCE generator assuming M(A) = 1.0 GeV","E#nu [GeV]","#sigma/E#nu [10^{-38} cm^{2} / GeV]");
    gXsec -> Draw("apc");
    
    c -> cd(2);
    plot.SetFrame(gXsecE,"#nu - n CCQE Xsec, NUANCE generator assuming M(A) = 1.0 GeV","E#nu [GeV]","#sigma [10^{-38} cm^{2}]");
    gXsecE -> Draw("apc");
    
    c -> SaveAs("~/Desktop/Xsec.pdf");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::DrawEfluxGraph(){
    
    
    if(gEflux && gEfluxE){
        TCanvas * c = plot.CreateCanvas("Xsec","Divide",2,1);
        c -> cd(1);
        plot.SetFrame(gEflux,"#nu - flux","E#nu [GeV]","#phi / 50 MeV / m ^{2} / 10 ^{6} POT");
        gEflux -> Draw("ap");
        c -> cd(2);
        plot.SetFrame(gEfluxE,"#nu - flux","E#nu [GeV]","flux [a.u.]");
        gEfluxE -> Draw("ap");
        plot.SetFrame(hEflux,"#nu - flux ","E#nu [GeV]","flux [a.u.]",4,38);
        hEflux -> Draw("same HIST");
        c -> SaveAs("~/Desktop/Eflux.pdf");
    }
    else {
        TCanvas * c = plot.CreateCanvas("Xsec");
        plot.SetFrame(hEflux,"#nu - flux ","E#nu [GeV]","flux [a.u.]",4,38);
        hEflux -> Draw("same HIST");
        c -> SaveAs("~/Desktop/Eflux.pdf");
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnuNCC::DrawMomentumDist(){
    
    TCanvas * c = plot.CreateCanvas("MomentumDist","Divide",2,1);
    
    c -> cd(1);
    plot.SetFrame(hFG,"n momentum dist. - Fermi Gas","neutron momentum [GeV/c]","#rho [a.u.]",1,46);
    hFG -> Draw();
    
    c -> cd(2);
    plot.SetFrame(hCFG,"n momentum dist. - Correlated Fermi Gas","neutron momentum [GeV/c]","#rho [a.u.]", 1,46);
    hCFG -> Draw();
    
    c -> SaveAs("~/Desktop/MomentumDist.pdf");
}




#endif

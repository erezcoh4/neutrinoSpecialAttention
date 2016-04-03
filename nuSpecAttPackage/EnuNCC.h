/**
 * \file EnuNCC.h
 *
 * \ingroup nuSpecAttPackage
 * 
 * \brief Class def header for a class EnuNCC
 *
 * @author erezcohen
 */

/** \addtogroup nuSpecAttPackage

    @{*/
#ifndef ENUNCC_H
#define ENUNCC_H

#include <iostream>
#include "MySoftwarePackage/TPlots.h"
#include "MySoftwarePackage/TCalculations.h"
#include "MySoftwarePackage/TAnalysis.h"
#include "TRandom3.h"

/**
   \class EnuNCC
   User defined class EnuNCC ... these comments are used to generate
   doxygen documentation!
 */
class EnuNCC{

public:
    
    TPlots    plot;
    TAnalysis analysis;
    TRandom3  rand;
    
    TTree   * OutTree;
    TGraph  * gXsec , * gXsecE  , * gEflux, * gEfluxE;
    TH1F    * hFG   ,  * hCFG   , * hEflux;
    TF1     * SRCk4Tail;
    
    bool nInSRC;
    
    Double_t Ev          , Pn            , Px       , Py    , Pz    , En;
    Double_t CMenergy    , XsecWeight    , Ev_INnRF , XsecLabFrame ;
    
    
    
    TVector3        Pcm ; // c.m. motion of the pair
    TLorentzVector  nu  , n  , nu_INnRF , W , mu , p , prec;
    TString         nDirection          , precDirection;
    
    /// Default constructor
    EnuNCC(){}
    EnuNCC( TTree * fOutTree );
    
    /// Default destructor
    ~EnuNCC(){}

    
    
    
    // methods
    void      ImpXsecGraph ( TString XsecFileName , int , bool DoPlot = false );
    void     DrawXsecGraph ();

    void   ImpMomentumDist ( bool DoPlot = false );
    void  DrawMomentumDist ();
    
    
    void     ImpEfluxHisto ( TH1F* ,bool DoPlot = false  );
    void     ImpEfluxGraph ( TString , int ,bool DoPlot = false );
    void    DrawEfluxGraph ();
    
    
    // output
    void       InitOutTree ();
    
    // GETs
    TGraph *      GetXsecGraph ( )       {return gXsec;};
    TGraph *     GetXsecEGraph ( )       {return gXsecE;};
    TGraph *          GetEflux ( )       {return gEflux;};
    TH1F *   GetMomentumDistFG ( )       {return hFG;};
    TH1F *  GetMomentumDistCFG ( )       {return hCFG;};
    
    
    // SETs
    void            SetOutTree ( TTree * fOutTree ) {OutTree = fOutTree;};
    void          SetXsecGraph ( TGraph * g )       {gXsec = g;};
    void         SetXsecEGraph ( TGraph * g )       {gXsecE = g;};
    void         SetEfluxGraph ( TGraph * g )       {gEflux = g;};
    void        SetEfluxEGraph ( TGraph * g )       {gEfluxE = g;};
    void     SetMomentumDistFG ( TH1F * h )         {hFG = h;};
    void    SetMomentumDistCFG ( TH1F * h )         {hCFG = h;};
    
    
    
    // interactions
    void       RunInteractions ( TString Nmodel = "CFG" ,  TString nuFlux  = "BNB flux" , int Ninteractions = 1 , bool DoPrint = false );
    void      GenerateNeutrino ( TString nuFlux = "BNB flux" );
    void       GenerateNeutron (TString Nmodel = "CFG");
    void       CalcRestFrameEv ();
    void  GenerateRecoilProton ( float mean = 0.0 , float sigma = 0.0 );

    
    void             PrintDATA (int);
    

    

    
};

#endif
/** @} */ // end of doxygen group 


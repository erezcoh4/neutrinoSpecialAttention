import ROOT
from ROOT import EnuNCC , TAnalysis
from rootpy.interactive import wait
analysis = TAnalysis()

# To run this script:
#   run with some PpPmCut (400/600/800) which is a cut on mcc6 events for |p(p)|+|p(mu)| in order to constrain the Ev
#   run nBack and then nForward, and then use ana_nu_n.py to compare them and look at the asymmetry


DoDraw  = False
NModel  = "CFGnBack"
#NModel = "CFGnForward"

nuFlux = "Ev from mcc6"
#nuFlux = "monochromatic neutrino 300 MeV"

Path = "/Users/erezcohen/Desktop/uboone/AnaFiles/"



PpPmuMin = [0   , 400 , 600 , 800 ]
PpPmuMax = [400 , 600 , 800 , 3000]
hEflux   = []
for i in range(0,len(PpPmuMin)):
    hEflux.append(analysis.GetHistoFromAFile("/Users/erezcohen/Desktop/uBoone/SpecialAttention/Data/MCC6recEventsEv.root"
                                        , "hEflux_%dPpPlusPmu%d"%(int(PpPmuMin[i]),int(PpPmuMax[i]))))
    OutFile = ROOT.TFile(Path+"CCinteractions"+NModel+"_%d_PpPmu_%d"%(int(PpPmuMin[i]),int(PpPmuMax[i]))+".root","recreate")
    OutTree = ROOT.TTree("anaTree","nu-n QE" + "("+NModel + ")" + ", " + "%.0f<p(p)+p(#mu)<%.0f MeV/c )"%(int(PpPmuMin[i]),int(PpPmuMax[i])))
    nuNCC   = EnuNCC( OutTree )
    nuNCC.ImpXsecGraph("/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Xsec.dat" , 70 , DoDraw )    # v-n cross section [A. Schukraft]
    if (DoDraw) : wait()
    nuNCC.ImpMomentumDist( DoDraw )    # neutron momentum distribution
    if (DoDraw) : wait()
    nuNCC.ImpEfluxHisto( hEflux[i] ,  DoDraw) #nuNCC.ImpEfluxGraph( "/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Eflux.dat" , 166 ,  DoDraw )    # BNB energy flux at uboone
    if (DoDraw) : wait()
    nuNCC.RunInteractions( NModel , nuFlux , 100000 , False )       # run interactions and fill output tree
    print "done filling %d events " % OutTree.GetEntries() + "in " + OutTree.GetTitle()
    OutTree.Write()
    OutFile.Close()

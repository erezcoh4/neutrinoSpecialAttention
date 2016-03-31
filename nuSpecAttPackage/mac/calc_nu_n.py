import ROOT
from ROOT import EnuNCC , TAnalysis
from rootpy.interactive import wait

# usage: run with some PpPmCut (400/600/800) which is a cut on mcc6 events for |p(p)|+|p(mu)| in order to constrain the Ev
# run nBack and then nForward, and then use ana_nu_n.py to compare them and look at the asymmetry

analysis = TAnalysis()
Ninteractions = 100000
NModel = "CFGnBack"
#NModel = "CFGnForward"
#nuFlux = "monochromatic neutrino 300 MeV"
PpPmuCut = 800
nuFlux = "mcc6 |p(p)|+|p(mu)|<%d" % PpPmuCut
hEflux = analysis.GetHistoFromAFile("/Users/erezcohen/Desktop/uBoone/SpecialAttention/Data/MCC6recEventsEv.root", "hEflux_PpPlusPmuBelow%d"% PpPmuCut )




DoDraw = False
Path = "/Users/erezcohen/Desktop/uboone/SpecialAttention/AnaFiles/"
OutFile = ROOT.TFile(Path+"CCinteractions"+NModel+"_PpPmuCut%d"%PpPmuCut+".root","recreate")
OutTree = ROOT.TTree("anaTree","sEG for nu-n QE interaction" + "("+NModel + ")" + ", " + nuFlux)


nuNCC   = EnuNCC( OutTree )

nuNCC.ImpXsecGraph("/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Xsec.dat" , 70 , DoDraw )    # neutrino - n cross section by A. Schukraft
print 'generated Xsec'
if (DoDraw) : wait()

nuNCC.ImpMomentumDist( DoDraw )    # neutron momentum distribution
print 'generated MomentumDist'
if (DoDraw) : wait()

nuNCC.ImpEfluxHisto( hEflux ,  DoDraw)
#nuNCC.ImpEfluxGraph( "/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Eflux.dat" , 166 ,  DoDraw )    # BNB energy flux at uboone
print 'generated Eflux'
if (DoDraw) : wait()


nuNCC.RunInteractions( NModel , nuFlux , Ninteractions , False )       # run interactions and fill output tree



print "done filling %d events " % OutTree.GetEntries() + "of " + NModel

OutTree.Write()
OutFile.Close()

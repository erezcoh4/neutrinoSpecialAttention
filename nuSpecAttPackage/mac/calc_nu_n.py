import ROOT
from ROOT import EnuNCC
from rootpy.interactive import wait

Ninteractions = 100000
#NModel = "CFGnBack"
NModel = "CFGnForward"

DoDraw = False
Path = "/Users/erezcohen/Desktop/uboone/SpecialAttention/AnaFiles/"
OutFile = ROOT.TFile(Path+"CCinteractions"+NModel+".root","recreate")
OutTree = ROOT.TTree("anaTree","sEG for nu-n QE interaction" + "("+NModel + ")")


nuNCC   = EnuNCC( OutTree )

nuNCC.ImpXsecGraph("/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Xsec.dat" , 70 , DoDraw )    # neutrino - n cross section by A. Schukraft
print 'generated Xsec'
if (DoDraw) : wait()

nuNCC.ImpMomentumDist( DoDraw )    # neutron momentum distribution
print 'generated MomentumDist'
if (DoDraw) : wait()

nuNCC.ImpEfluxGraph( "/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Eflux.dat" , 166 ,  DoDraw )    # BNB energy flux at uboone
print 'generated Eflux'
if (DoDraw) : wait()


nuNCC.RunInteractions( NModel , Ninteractions , False )       # run interactions and fill output tree



print "done filling %d events " % OutTree.GetEntries() + "of " + NModel

OutTree.Write()
OutFile.Close()

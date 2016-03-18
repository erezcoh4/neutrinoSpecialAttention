import ROOT
from ROOT import EnuNCC



OutFile = ROOT.TFile("$UBOONE/SpecialAttention/AnaFiles/CCinteractions.root","recreate")
OutTree = ROOT.TTree("anaTree","Event genation for neutrino on neutron QE interaction")


nuNCC   = EnuNCC( OutTree )

nuNCC.ImpXsecGraph("/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Xsec.dat" , 70 , False )    # neutrino - n cross section by A. Schukraft

nuNCC.ImpMomentumDist( False )    # neutron momentum distribution

nuNCC.ImpEflux( "/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Eflux.dat" , 70 ,  True )    # neutron momentum distribution


#nuNCC.SetMomentumDist() # get neutron momentum distribution
#nuNCC.SetNuflux()      # get neutrino energy flux
#
#nuNCC.RunInteractions( 100 )       # run interactions and fill output tree



print "done filling %d events" % OutTree.GetEntries()

OutTree.Write()
OutFile.Close()

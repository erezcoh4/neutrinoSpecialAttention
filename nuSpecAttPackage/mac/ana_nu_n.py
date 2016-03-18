import ROOT
from ROOT import TPlots
from ROOT import TAnalysis
from rootpy.interactive import wait

colorWeighted   = 2
colorGen        = 1
Nbins           = 100

NModel = "CFG"

ana = TPlots("/Users/erezcohen/Desktop/uboone/SpecialAttention/AnaFiles/CCinteractions"+NModel+".root" , "anaTree" , NModel )
analysis = TAnalysis()

c = ana.CreateCanvas("n momentum dist.","DivideSquare")

c.cd(1)
cut = ROOT.TCut()
hGen = ana.H1("n.P()",cut,"HIST",Nbins,0,0.7,"n momentum dist. (" + NModel + ")","n momentum [GeV/c]","",colorGen,colorGen)
#analysis.NormalizeHistogram(hGen)
cut = ROOT.TCut("XsecWeight")
hWeighted = ana.H1("n.P()",cut,"HIST same",Nbins,0,0.7,"","","",colorWeighted,colorWeighted)
#analysis.NormalizeHistogram(hWeighted)



c.cd(2)
cut = ROOT.TCut()
ana.H1("XsecWeight",cut,"HIST",Nbins,0,2,"cross-section weighting","E#nu [GeV]","",4,46)



c.cd(3)
cut = ROOT.TCut()
hnu = ana.H1("nu.E()",cut,"HIST",Nbins,0,2,"neutrino energy","E#nu [GeV]","",4,4)
hnu_INnRF = ana.H1("nu_INnRF.E()",cut,"HIST same",Nbins,0,2,"neutrino energy","E#nu [GeV]","",2,2)
#hEv_INnRF = ana.H1("Ev_INnRF",cut,"HIST same",Nbins,0,2,"neutrino energy","E#nu [GeV]","",5,53)
ana.AddLegend("",hnu,"lab frame neutrino",hnu_INnRF,"in neutron rest frame","F");

wait()


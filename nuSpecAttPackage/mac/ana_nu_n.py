import ROOT
from ROOT import TPlots
from ROOT import TAnalysis
from rootpy.interactive import wait
from ROOT import EnuNCC
ROOT.gStyle.SetOptStat(0000)

cutGen = ROOT.TCut()
colorGen        = 1
colorWeighted   = 2
cutW = ROOT.TCut("XsecWeight")

Nbins           = 50

NModel = "CFG"

ana = TPlots("/Users/erezcohen/Desktop/uboone/SpecialAttention/AnaFiles/CCinteractions"+NModel+".root" , "anaTree" , NModel )
analysis = TAnalysis()

c = ana.CreateCanvas("n momentum dist.","Divide",3,3)

c.cd(1)
hGen = ana.H1("n.P()",cutGen,"HIST",Nbins,0,0.7,"n momentum dist. (" + NModel + ")","n momentum [GeV/c]","",colorGen,colorGen)
Integral = 0
Tail = 0
k_F = 0.25
for bin in range (1,Nbins):
    k = hGen.GetXaxis().GetBinCenter( bin )
    Integral = Integral + k * k * hGen.GetBinContent( bin )
    if ( k > k_F ):
        Tail = Tail + k * k * hGen.GetBinContent( bin )

FracTailGen = Tail / Integral;
print "fraction in generated tail is %f" % FracTailGen
analysis.NormalizeHistogram(hGen)

hWeighted = ana.H1("n.P()",cutW,"HIST same",Nbins,0,0.7,"","","",colorWeighted,colorWeighted)
Integral = 0
Tail = 0
for bin in range (1,Nbins):
    k = hWeighted.GetXaxis().GetBinCenter( bin )
    Integral = Integral + k * k * hWeighted.GetBinContent( bin )
    if ( k > k_F ):
        Tail = Tail + k * k * hWeighted.GetBinContent( bin )

FracTailWeighted = Tail / Integral;
print "fraction in weighted tail is %f" % FracTailWeighted
analysis.NormalizeHistogram(hWeighted)
ana.AddLegend("n momentum dist.",hGen,"generated",hWeighted,"weighted","F");



c.cd(2)
ana.H2("nu.E()","XsecWeight",cutGen,"colz",Nbins,0,2,Nbins,0,0.7,"cross-section weighting","E#nu [GeV]","Xsec weight (a.u.)")
nuNCC = EnuNCC()
nuNCC.ImpXsecGraph("/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Xsec.dat" , 70 , False )    # neutrino - n cross section by A. Schukraft
g = nuNCC.GetXsecEGraph();
g.Draw("same c")


c.cd(3)
ana.H1("XsecWeight",cutGen,"HIST",Nbins,0,0.7,"cross-section weighting","Xsec weight (a.u.)","",4,46)



c.cd(4)
hEv = ana.H1("nu.E()",cutGen,"HIST",Nbins,0,2,"#nu energy (lab)","E#nu [GeV]","",colorGen,colorGen)
analysis.NormalizeHistogram(hEv)
hEvW = ana.H1("nu.E()",cutW,"HIST same",Nbins,0,2,"","E#nu [GeV]","",colorWeighted,colorWeighted)
analysis.NormalizeHistogram(hEvW)
ana.AddLegend("Lab Frame E#nu",hEv,"generated",hEvW,"weighted","F");


c.cd(5)
hEv_INnRF = ana.H1("nu_INnRF.E()",cutGen,"HIST",Nbins,0,2,"#nu energy in n RestFrame","E#nu [GeV]","",colorGen,colorGen)
analysis.NormalizeHistogram(hEv_INnRF)
hEvW_INnRF = ana.H1("nu_INnRF.E()",cutW,"HIST same",Nbins,0,2,"","E#nu [GeV]","",colorWeighted,colorWeighted)
analysis.NormalizeHistogram(hEvW_INnRF)
ana.AddLegend("Rest Frame E#nu",hEv_INnRF,"generated",hEvW_INnRF,"weighted","F");




c.cd(6)
ana.H1("nu_INnRF.E() - nu.E()",cutGen,"HIST",Nbins,-0.5,0.5,"boosting #nu","E#nu (Resf Frame) - E#nu (lab) [GeV]","",1,1)


c.cd(7)
ana.H2("n.Pz()","nu_INnRF.Pz() - nu.Pz()",cutGen,"colz",Nbins,-1,1,Nbins,-1,1,"cross-section weighting","p(n)^{z} [GeV/c]","#Delta p(#nu)^{z} [GeV/c]")


c.cd(8)
ana.H2("n.P()","XsecWeight",cutGen,"colz",Nbins,0,1,Nbins,0,0.7,"cross-section weighting","n momentum [GeV/c]","Xsec weight (a.u.)")



wait()


import ROOT
from ROOT import TPlots
from ROOT import TAnalysis
from rootpy.interactive import wait
from ROOT import EnuNCC

colorWeighted   = 2
colorGen        = 1
Nbins           = 50

NModel = "CFG"

ana = TPlots("/Users/erezcohen/Desktop/uboone/SpecialAttention/AnaFiles/CCinteractions"+NModel+".root" , "anaTree" , NModel )
analysis = TAnalysis()

c = ana.CreateCanvas("n momentum dist.","Divide",3,2)

c.cd(1)
cut = ROOT.TCut()
hGen = ana.H1("n.P()",cut,"HIST",Nbins,0,0.7,"n momentum dist. (" + NModel + ")","n momentum [GeV/c]","",colorGen,colorGen)
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
#analysis.NormalizeHistogram(hGen)

cut = ROOT.TCut("XsecWeight")
hWeighted = ana.H1("n.P()",cut,"HIST same",Nbins,0,0.7,"","","",colorWeighted,colorWeighted)
Integral = 0
Tail = 0
for bin in range (1,Nbins):
    k = hWeighted.GetXaxis().GetBinCenter( bin )
    Integral = Integral + k * k * hWeighted.GetBinContent( bin )
    if ( k > k_F ):
        Tail = Tail + k * k * hWeighted.GetBinContent( bin )

FracTailWeighted = Tail / Integral;
print "fraction in weighted tail is %f" % FracTailWeighted
#analysis.NormalizeHistogram(hWeighted)
ana.AddLegend("n momentum dist.",hGen,"generated",hWeighted,"weighted","F");



c.cd(2)
cut = ROOT.TCut()
ana.H2("nu.E()","XsecWeight",cut,"colz",Nbins,0,2,Nbins,0,0.7,"cross-section weighting","E#nu [GeV]","Xsec weight (a.u.)")
nuNCC = EnuNCC()
nuNCC.ImpXsecGraph("/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Xsec.dat" , 70 , False )    # neutrino - n cross section by A. Schukraft
g = nuNCC.GetXsecEGraph();
g.Draw("same c")

c.cd(3)
cut = ROOT.TCut()
ana.H1("XsecWeight",cut,"HIST",Nbins,0,0.7,"cross-section weighting","Xsec weight (a.u.)","",4,46)



c.cd(4)
cut = ROOT.TCut()
hEv = ana.H1("nu.E()",cut,"HIST",Nbins,0,2,"#nu energy (lab)","E#nu [GeV]","",colorGen,colorGen)
cut = ROOT.TCut("XsecWeight")
hEvW = ana.H1("nu.E()",cut,"HIST same",Nbins,0,2,"","E#nu [GeV]","",colorWeighted,colorWeighted)
ana.AddLegend("Lab Frame E#nu",hEv,"generated",hEvW,"weighted","F");


c.cd(5)
cut = ROOT.TCut()
hEv_INnRF = ana.H1("nu_INnRF.E()",cut,"HIST",Nbins,0,2,"#nu energy in n RestFrame","E#nu [GeV]","",colorGen,colorGen)
cut = ROOT.TCut("XsecWeight")
hEvW_INnRF = ana.H1("nu_INnRF.E()",cut,"HIST same",Nbins,0,2,"","E#nu [GeV]","",colorWeighted,colorWeighted)
ana.AddLegend("Rest Frame E#nu",hEv_INnRF,"generated",hEvW_INnRF,"weighted","F");




c.cd(6)
cut = ROOT.TCut()
ana.H1("nu_INnRF.E() - nu.E()",cut,"HIST",Nbins,-0.5,0.5,"boosting #nu","E#nu (Resf Frame) - E#nu (lab) [GeV]","",1,1)
#ana.H1("nu_INnRF.Pz() - nu.Pz()",cut,"HIST",Nbins,-0.5,0.5,"boosting #nu","p(#nu) (Resf Frame) - p(#nu) (lab) [GeV/c]","",1,1)


wait()


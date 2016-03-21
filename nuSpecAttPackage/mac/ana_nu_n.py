import ROOT
from ROOT import TPlots
from ROOT import TAnalysis
from rootpy.interactive import wait
from ROOT import EnuNCC
ROOT.gStyle.SetOptStat(0000)

cutGen = ROOT.TCut()
colorGen        = 1
colorWeighted   = 2
cutW            = ROOT.TCut("XsecWeight")
DoAllVariables  = True
DoAsymmetery    = False
Nbins           = 25

#NModel = "CFG"
NModel = "CFGnBack"

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
    Integral = Integral +  hGen.GetBinContent( bin )
    if ( k > k_F ):
        Tail = Tail +  hGen.GetBinContent( bin )

FracTailGen = Tail / (Integral-Tail);
print "fraction in generated tail is %f" % FracTailGen
analysis.NormalizeHistogram(hGen)

hWeighted = ana.H1("n.P()",cutW,"HIST same",Nbins,0,0.7,"","","",colorWeighted,colorWeighted)
Integral = 0
Tail = 0
for bin in range (1,Nbins):
    k = hWeighted.GetXaxis().GetBinCenter( bin )
    Integral = Integral + hWeighted.GetBinContent( bin )
    if ( k > k_F ):
        Tail = Tail + hWeighted.GetBinContent( bin )

FracTailWeighted = Tail / (Integral-Tail);
print "fraction in weighted tail is %f" % FracTailWeighted
analysis.NormalizeHistogram(hWeighted)
ana.AddLegend("n momentum dist.",hGen,"generated",hWeighted,"weighted","F");



c.cd(2)
ana.H2("nu.E()","XsecWeight",cutGen,"colz",Nbins,0.18,0.7,Nbins,0,0.7,"cross-section weighting","E#nu [GeV]","Xsec weight (a.u.)")
nuNCC = EnuNCC()
nuNCC.ImpXsecGraph("/Users/erezcohen/Desktop/uboone/SpecialAttention/Data/Xsec.dat" , 70 , False )    # neutrino - n cross section by A. Schukraft
g = nuNCC.GetXsecEGraph();
g.Draw("same c")


c.cd(3)
ana.H2("n.Pz()","XsecWeight",cutGen,"colz",Nbins,-0.5,0.5,Nbins,0.01,0.7
       ,"cross-section weighting vs. neutron Pz","p(n) ^{z} [GeV/c]","Xsec weight (a.u.)")



c.cd(4)
hEv = ana.H1("nu.E()",cutGen,"HIST",Nbins,0.18,0.7,"#nu energy (lab)","E#nu [GeV]","",colorGen,colorGen)
analysis.NormalizeHistogram(hEv)
hEvW = ana.H1("nu.E()",cutW,"HIST same",Nbins,0.18,0.7,"","E#nu [GeV]","",colorWeighted,colorWeighted)
analysis.NormalizeHistogram(hEvW)
ana.AddLegend("Lab Frame E#nu",hEv,"generated",hEvW,"weighted","F");


c.cd(5)
hEv_INnRF = ana.H1("nu_INnRF.E()",cutGen,"HIST",Nbins,0.18,0.7,"#nu energy in n RestFrame","E#nu [GeV]","",colorGen,colorGen)
analysis.NormalizeHistogram(hEv_INnRF)
hEvW_INnRF = ana.H1("nu_INnRF.E()",cutW,"HIST same",Nbins,0.18,0.7,"","E#nu [GeV]","",colorWeighted,colorWeighted)
analysis.NormalizeHistogram(hEvW_INnRF)
ana.AddLegend("Rest Frame E#nu",hEv_INnRF,"generated",hEvW_INnRF,"weighted","F");




c.cd(6)
ana.H1("nu_INnRF.E() - nu.E()",cutGen,"HIST",Nbins,-0.5,0.5,"boosting #nu","E#nu (Resf Frame) - E#nu (lab) [GeV]","",1,1)


c.cd(7)
ana.H2("n.Pz()","nu_INnRF.Pz() - nu.Pz()",cutGen,"colz",Nbins,-1,1,Nbins,-1,1,"cross-section weighting","p(n)^{z} [GeV/c]","#Delta p(#nu)^{z} [GeV/c]")


c.cd(8)
ana.H2("n.P()","XsecWeight",cutGen,"colz",Nbins,0,1,Nbins,0,0.7,"cross-section weighting","n momentum [GeV/c]","Xsec weight (a.u.)")


c.cd(9)
ana.H1("(XsecWeight - XsecLabFrame)/XsecWeight"
       ,cutGen,"HIST",Nbins,-1,1,"cc difference (normalized)","[Xsec(RF)-Xsec(lab)]/Xsec(RF)","",46,46)

c.SaveAs("~/Desktop/"+c.GetName()+".pdf")
wait()




if (DoAsymmetery) :

    cutAsym = ROOT.TCut("sign(n.Pz())")
    cAsym   = ana.CreateCanvas("cAsym","Divide",3,1)

    cAsym.cd(1)
    hAsym   = ana.H2("n.P()","nu.E()",cutGen,"colz",Nbins,0,1,Nbins,0.2,0.4,"generated","n momentum [GeV/c]","E#nu (lab) [GeV]")


    cAsym.cd(2)
    hAsym   = ana.H2("n.P()","nu.E()",cutW,"colz",Nbins,0,1,Nbins,0.2,0.4,"weighted","n momentum [GeV/c]","E#nu (lab) [GeV]")


    cAsym.cd(3)
#hAsym   = ana.H2("n.P()","nu.E()",cutAsym,"colz",Nbins,0,1,Nbins,0,1,"asym.","n momentum [GeV/c]","E#nu (lab) [GeV]")


    hAsym   = analysis.Assymetry(ana.GetTree() , "n" , "n" , 25 ,0 , 1 , "nu", 25 ,0.2,0.4 , "XsecWeight" , False )
    ana.SetFrame(hAsym,"(backward - forward) / (backward + forward) asymmetry","n momentum [GeV/c]","E#nu (lab) [GeV]")
    hAsym.Draw("colz")


    cAsym.SaveAs("~/Desktop/"+c.GetName()+".pdf")
    wait()













import ROOT
from ROOT import TPlots
from ROOT import TAnalysis
from rootpy.interactive import wait
from ROOT import EnuNCC
ROOT.gStyle.SetOptStat(0000)


DoAllVariables  = False
DoAsymmetery    = False
DoMomentumDist  = True


cutGen          = ROOT.TCut()
colorGen        = 1
colorWeighted   = 2
cutW            = ROOT.TCut("XsecWeight")
Nbins           = 25
kF              = 0.25

#NModel = "CFG"
#NModel = "CFGnBack"
NModel = "CFGnForward"
PpPmuCut = 800
nuFlux = "mcc6 |p(p)|+|p(mu)|<%d" % PpPmuCut

#ana = TPlots("/Users/erezcohen/Desktop/uboone/SpecialAttention/AnaFiles/CCinteractions"+NModel+".root" , "anaTree" , NModel )
anaBack     = TPlots("/Users/erezcohen/Desktop/uboone/SpecialAttention/AnaFiles/CCinteractions"+"CFGnBack"+"_PpPmuCut%d"%PpPmuCut+".root" , "anaTree" , "CFGnBack" )
anaForward  = TPlots("/Users/erezcohen/Desktop/uboone/SpecialAttention/AnaFiles/CCinteractions"+"CFGnForward"+"_PpPmuCut%d"%PpPmuCut+".root" , "anaTree" , "CFGnForward" )
analysis = TAnalysis()



if DoAllVariables:
    c = ana.CreateCanvas("n momentum dist.","Divide",3,3)

    c.cd(1)
    hGen = ana.H1("n.P()",cutGen,"HIST",Nbins,0,0.7,"n momentum dist. (" + NModel + ")","n momentum [GeV/c]","",colorGen,colorGen)
#    Integral = 0
#    Tail = 0
#    for bin in range (1,Nbins):
#        k = hGen.GetXaxis().GetBinCenter( bin )
#        Integral = Integral +  hGen.GetBinContent( bin )
#        if ( k > k_F ):
#            Tail = Tail +  hGen.GetBinContent( bin )
#    FracTailGen = Tail / (Integral-Tail);
#print "fraction in generated tail is %f" % FracTailGen
    analysis.NormalizeHistogram(hGen)

    hWeighted = ana.H1("n.P()",cutW,"HIST same",Nbins,0,0.7,"","","",colorWeighted,colorWeighted)
#Integral = 0
#Tail = 0
#for bin in range (1,Nbins):
#    k = hWeighted.GetXaxis().GetBinCenter( bin )
#    Integral = Integral + hWeighted.GetBinContent( bin )
#    if ( k > k_F ):
#        Tail = Tail + hWeighted.GetBinContent( bin )
#
#FracTailWeighted = Tail / (Integral-Tail);
#print "fraction in weighted tail is %f" % FracTailWeighted
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
    c.Update()
    #c.SaveAs("~/Desktop/"+c.GetName()+".pdf")
    wait()


if DoAsymmetery:

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


# March - 29
if DoMomentumDist:

    cMomentumDist = anaBack.CreateCanvas("n momentum dist.", "Divide" , 3 , 3)
    
    cMomentumDist.cd(1)
    hGen = anaBack.H1("n.P()",cutGen,"HIST",Nbins,0,0.7,"n momentum dist. (" + NModel + ")","n momentum [GeV/c]","",colorGen,colorGen)
    hWtd = anaBack.H1("n.P()",cutW,"HIST same",Nbins,0,0.7,"","","",colorWeighted,colorWeighted)
    
    
    TailGen = 100 * hGen.Integral(hGen.GetXaxis().FindBin(0.26),hGen.GetNbinsX()) / hGen.Integral()
    TailWtd = 100 * hWtd.Integral(hWtd.GetXaxis().FindBin(0.26),hWtd.GetNbinsX()) / hWtd.Integral()

    kF_bin = hGen.GetXaxis().FindBin(0.25)
    
    hTailGen = ROOT.TH1F("hTailGen","",Nbins,0,0.7)
    anaBack.SetFrame(hTailGen ,"n momentum distribution tail","n momentum [GeV/c]","% of mean field events",colorGen,0)
    hTailWtd = ROOT.TH1F("hTailWtd","",Nbins,0,0.7)
    anaBack.SetFrame(hTailWtd ,"n momentum distribution tail","n momentum [GeV/c]","% of mean field events" ,colorWeighted,0)
    hTlGn_k4 = ROOT.TH1F("hTlGn_k4","",Nbins,0,0.7)
    anaBack.SetFrame(hTlGn_k4 ,"tail times k ^{4}","n momentum [GeV/c]","MF% #times k ^{4} [a.u.]",colorGen,0)
    hTlWd_k4 = ROOT.TH1F("hTlWd_k4","",Nbins,0,0.7)
    anaBack.SetFrame(hTlWd_k4 ,"tail times k ^{4}","n momentum [GeV/c]","MF% #times k ^{4} [a.u.]" ,colorWeighted,0)
    hTailRat = ROOT.TH1F("hTailRat","",Nbins,0,0.7)
    anaBack.SetFrame(hTailRat ,"ratio of increase in tail","n momentum [GeV/c]","% of generated" , 4 , 38)
   
    for bin in range(kF_bin,Nbins):
        k = hGen.GetXaxis().GetBinCenter(bin)
        hTailGen.SetBinContent( bin , 100 * hGen.GetBinContent(bin) / hGen.Integral() )
        hTailWtd.SetBinContent( bin , 100 * hWtd.GetBinContent(bin) /  hWtd.Integral() )
        hTlGn_k4.SetBinContent( bin , 100 * pow( k , 4 ) * hGen.GetBinContent(bin) / hGen.Integral() )
        hTlWd_k4.SetBinContent( bin , 100 * pow( k , 4 ) * hWtd.GetBinContent(bin) /  hWtd.Integral() )
        hTailRat.SetBinContent( bin , 100 * hTailWtd.GetBinContent(bin) / hTailGen.GetBinContent(bin) )
    hTailWtd.Draw()
    hTailGen.Draw("same")
    anaBack.AddLegend("n momentum dist.",hTailGen,"generated (%.1f%%)"%TailGen
                                        ,hTailWtd,"weighted (%.1f%%)"%TailWtd,"l");
    cMomentumDist.cd(2)
    hTlWd_k4.Draw("")
    hTlGn_k4.Draw("same")


    cMomentumDist.cd(3)
    hTailRat.Draw()


    cMomentumDist.cd(4)
    hGenF = anaForward.H1("n.P()",cutGen,"HIST",Nbins,0,0.7,"n momentum dist. (" + NModel + ")","n momentum [GeV/c]","",colorGen,colorGen)
    hWtdF = anaForward.H1("n.P()",cutW,"HIST same",Nbins,0,0.7,"","","",colorWeighted,colorWeighted)
    
    

    TailGenF = 100 * hGenF.Integral(hGenF.GetXaxis().FindBin(0.26),hGenF.GetNbinsX()) / hGenF.Integral()
    TailWtdF = 100 * hWtdF.Integral(hWtdF.GetXaxis().FindBin(0.26),hWtdF.GetNbinsX()) / hWtdF.Integral()


    hTailGenF = ROOT.TH1F("hTailGenF","",Nbins,0,0.7)
    anaForward.SetFrame(hTailGenF ,"n momentum distribution tail","n momentum [GeV/c]","% of mean field events",colorGen,0)
    hTailWtdF = ROOT.TH1F("hTailWtdF","",Nbins,0,0.7)
    anaForward.SetFrame(hTailWtdF ,"n momentum distribution tail","n momentum [GeV/c]","% of mean field events" ,colorWeighted,0)
    hTlGn_k4F = ROOT.TH1F("hTlGn_k4F","",Nbins,0,0.7)
    anaForward.SetFrame(hTlGn_k4F ,"tail times k ^{4}","n momentum [GeV/c]","MF% #times k ^{4} [a.u.]",colorGen,0)
    hTlWd_k4F = ROOT.TH1F("hTlWd_k4F","",Nbins,0,0.7)
    anaForward.SetFrame(hTlWd_k4F ,"tail times k ^{4}","n momentum [GeV/c]","MF% #times k ^{4} [a.u.]" ,colorWeighted,0)
    hTailRatF = ROOT.TH1F("hTailRatF","",Nbins,0,0.7)
    anaForward.SetFrame(hTailRatF ,"ratio of increase in tail","n momentum [GeV/c]","% of generated" , 4 , 38)
    
    for bin in range(kF_bin,Nbins):
        k = hGen.GetXaxis().GetBinCenter(bin)
        hTailGenF.SetBinContent( bin , 100 * hGenF.GetBinContent(bin) / hGenF.Integral() )
        hTailWtdF.SetBinContent( bin , 100 * hWtdF.GetBinContent(bin) / hWtdF.Integral() )
        hTlGn_k4F.SetBinContent( bin , 100 * pow( k , 4 ) * hGenF.GetBinContent(bin) / hGenF.Integral() )
        hTlWd_k4F.SetBinContent( bin , 100 * pow( k , 4 ) * hWtdF.GetBinContent(bin) / hWtdF.Integral() )
        hTailRatF.SetBinContent( bin , 100 * hTailWtdF.GetBinContent(bin) / hTailGenF.GetBinContent(bin) )
    hTailGenF.Draw()
    hTailWtdF.Draw("same")
    anaForward.AddLegend("n momentum dist.",hTailGenF,"generated (%.1f%%)"%TailGenF ,hTailWtdF,"weighted (%.1f%%)"%TailWtdF,"l");
    cMomentumDist.cd(5)
    hTlGn_k4F.Draw("")
    hTlWd_k4F.Draw("same")
    cMomentumDist.cd(6)
    hTailRatF.Draw()

    # neutron before reaction
    cMomentumDist.cd(7)
    hTailFBasymGen = ROOT.TH1F("hTailFBasymGen","",Nbins,0,0.7)
    anaForward.SetFrame(hTailFBasymGen ,"forward/backward asym.","n momentum [GeV/c]","Backward n / Forwad n [%]",colorGen,0)
    hTailFBasymWtd = ROOT.TH1F("hTailFBasymWtd","",Nbins,0,0.7)
    anaForward.SetFrame(hTailFBasymWtd ,"forward/backward asym.","n momentum [GeV/c]","Backward n / Forwad n [%]",colorWeighted,0)
    for bin in range(kF_bin,Nbins):
        hTailFBasymGen.SetBinContent( bin , 100 * hGen.GetBinContent(bin) / hGenF.GetBinContent(bin) )
        hTailFBasymWtd.SetBinContent( bin , 100 * hWtd.GetBinContent(bin) / hWtdF.GetBinContent(bin) )


    hTailFBasymWtd.Draw("")
    hTailFBasymGen.Draw("same")


    # recoil proton
    cMomentumDist.cd(8)
    hGenB_p = anaBack.H1("p.P()",cutGen,"HIST",Nbins,0,0.7,"p momentum dist. (forward)","p momentum [GeV/c]","",colorGen,colorGen)
    hWtdB_p = anaBack.H1("p.P()",cutW,"HIST same",Nbins,0,0.7,"","","",colorWeighted,colorWeighted)
    hGenF_p = anaForward.H1("p.P()",cutGen,"HIST",Nbins,0,0.7,"p momentum dist. (forward)","p momentum [GeV/c]","",colorGen,colorGen)
    hWtdF_p = anaForward.H1("p.P()",cutW,"HIST same",Nbins,0,0.7,"","","",colorWeighted,colorWeighted)
    hTailFBasymGen_p = ROOT.TH1F("hTailFBasymGen_p","",Nbins,0,0.7)
    anaForward.SetFrame(hTailFBasymGen_p ,"forward/backward asym.","p(rec) momentum [GeV/c]","Backward p(rec) / Forwad p(rec) [%]",colorGen,0)
    hTailFBasymWtd_p = ROOT.TH1F("hTailFBasymWtd_p","",Nbins,0,0.7)
    anaForward.SetFrame(hTailFBasymWtd_p ,"forward/backward asym.","p(rec) momentum [GeV/c]","Backward p(rec) / Forwad p(rec) [%]",colorWeighted,0)
    for bin in range(kF_bin,Nbins):
        hTailFBasymGen_p.SetBinContent( bin , 100 * hGenB_p.GetBinContent(bin) / hGenF_p.GetBinContent(bin) )
        hTailFBasymWtd_p.SetBinContent( bin , 100 * hWtdB_p.GetBinContent(bin) / hWtdF_p.GetBinContent(bin) )


    hTailFBasymWtd_p.Draw("")
    hTailFBasymGen_p.Draw("same")

    cMomentumDist.cd(9)
    anaBack.H1("nu.E()",cutGen,"HIST",100,0,2.,"neutrino generated energy","E#nu [GeV]","",1,1)

    cMomentumDist.Update()
    wait()









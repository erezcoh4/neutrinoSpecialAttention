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
cutW            = ROOT.TCut("XsecWeight")
colorGen        = 1
colorWtdB       = 2
colorWtdF       = 4
colorAsym       = 6
Nbins           = 25
kF              = 0.25


PpPmuBin    = 4
PpPmuMin    = [0   , 400 , 600 , 800  , 1000]
PpPmuMax    = [400 , 600 , 800 , 1000 , 3000]
hEflux      = []
anaBack     = []
anaFrwd     = []
Path        = "/Users/erezcohen/Desktop/uboone/SpecialAttention/AnaFiles"
plot        = TPlots()
analysis    = TAnalysis()
for i in range(0,len(PpPmuMin)):
    hEflux.append(analysis.GetHistoFromAFile("/Users/erezcohen/Desktop/uBoone/SpecialAttention/Data/MCC6recEventsEv.root"
                                             , "hEflux_%dPpPlusPmu%d"%(int(PpPmuMin[i]),int(PpPmuMax[i]))))
    anaBack.append(TPlots(Path+"/CCinteractionsCFGnBack_%d_PpPmu_%d"%(int(PpPmuMin[i]),int(PpPmuMax[i]))+".root", "anaTree" , "nB%d"%i ))
    anaFrwd.append(TPlots(Path+"/CCinteractionsCFGnForward_%d_PpPmu_%d"%(int(PpPmuMin[i]),int(PpPmuMax[i]))+".root", "anaTree" , "nF%d"%i ))


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

    anaF    = anaFrwd[PpPmuBin]
    anaB    = anaBack[PpPmuBin]
    Nevents = hEflux[PpPmuBin].Integral()
    
    cMomentumDist = anaB.CreateCanvas("n momentum dist.", "Divide" , 2 , 3)
    
    cMomentumDist.cd(1)
    anaB.H1("nu.E()",cutGen,"HIST",100,0,2.,"neutrino generated energy","E#nu [GeV]","",1,1)

#    cMomentumDist.cd(2)
    hGen_n  = anaB.H1("n.P()",cutGen,"goff",Nbins,0,0.7,"n momentum dist.","n momentum [GeV/c]","",colorGen,colorGen)
    hWtdB_n = anaB.H1("n.P()",cutW,"goff",Nbins,0,0.7,"","","",colorWtdB,colorWtdB)
    hWtdF_n = anaF.H1("n.P()",cutW,"goff",Nbins,0,0.7,"","","",colorWtdF,colorWtdF)
    hGen_p  = anaB.H1("p.P()",cutGen,"goff",Nbins,0,0.7,"p momentum dist.","p momentum [GeV/c]","",colorGen,colorGen)
    hWtdB_p = anaB.H1("p.P()",cutW,"goff",Nbins,0,0.7,"","","",colorWtdB,colorWtdB)
    hWtdF_p = anaF.H1("p.P()",cutW,"goff",Nbins,0,0.7,"","","",colorWtdF,colorWtdF)
    
    
    
    TailGen_n     = 100 * hGen_n.Integral(hGen_n.GetXaxis().FindBin(0.26),hGen_n.GetNbinsX()) / hGen_n.Integral()
    NTailGen_n  = int(TailGen_n*Nevents/100)
    TailGen_p     = 100 * hGen_p.Integral(hGen_p.GetXaxis().FindBin(0.26),hGen_p.GetNbinsX()) / hGen_p.Integral()
    NTailGen_p  = int(TailGen_p*Nevents/100)
    TailWtdB_n    = 100 * hWtdB_n.Integral(hWtdB_n.GetXaxis().FindBin(0.26),hWtdB_n.GetNbinsX()) / hWtdB_n.Integral()
    NTailWtdB_n = int(TailWtdB_n*Nevents/100)
    TailWtdF_n    = 100 * hWtdF_n.Integral(hWtdF_n.GetXaxis().FindBin(0.26),hWtdF_n.GetNbinsX()) / hWtdF_n.Integral()
    NTailWtdF_n = int(TailWtdF_n*Nevents/100)
    print "%d +/- %d events, n > kF %d, p(rec) > kF %d events" % (int(Nevents) , int(ROOT.sqrt(Nevents)) , int(NTailGen_n) , int(NTailGen_p))



    hTailGen_n  = ROOT.TH1F("hTailGen_n","",Nbins,0,0.7)
    plot.SetFrame(hTailGen_n ,"n momentum distribution tail","n momentum [GeV/c]","% of mean field events",colorGen,colorGen)
    hTailWtdB_n = ROOT.TH1F("hTailWtdB_n","",Nbins,0,0.7)
    plot.SetFrame(hTailWtdB_n ,"n momentum distribution tail","n momentum [GeV/c]","% of mean field events" ,colorWtdB,colorWtdB)
    hTailWtdF_n = ROOT.TH1F("hTailWtdF_n","",Nbins,0,0.7)
    plot.SetFrame(hTailWtdF_n ,"n momentum distribution tail","n momentum [GeV/c]","% of mean field events" ,colorWtdF,colorWtdF)
    hTailGen_n_k4 = ROOT.TH1F("hTailGen_n_k4","",Nbins,0,0.7)
    plot.SetFrame(hTailGen_n_k4 ,"tail times k ^{4}","n momentum [GeV/c]","MF% #times k ^{4} [a.u.]",colorGen,colorGen)
    hTailWtdB_n_k4 = ROOT.TH1F("hTailWtdB_n_k4","",Nbins,0,0.7)
    plot.SetFrame(hTailWtdB_n_k4 ,"tail times k ^{4}","n momentum [GeV/c]","MF% #times k ^{4} [a.u.]" ,colorWtdB,colorWtdB)
    hTailWtdF_n_k4 = ROOT.TH1F("hTailWtdF_n_k4","",Nbins,0,0.7)
    plot.SetFrame(hTailWtdF_n_k4 ,"tail times k ^{4}","n momentum [GeV/c]","MF% #times k ^{4} [a.u.]" ,colorWtdF,colorWtdF)
   
    for bin in range(hGen_n.GetXaxis().FindBin(kF),Nbins):
        k = hTailGen_n.GetXaxis().GetBinCenter(bin)
        hTailGen_n.SetBinContent( bin , 100 * hGen_n.GetBinContent(bin) / hGen_n.Integral() )
        hTailWtdB_n.SetBinContent( bin , 100 * hWtdB_n.GetBinContent(bin) /  hWtdB_n.Integral() )
        hTailWtdF_n.SetBinContent( bin , 100 * hWtdF_n.GetBinContent(bin) /  hWtdF_n.Integral() )
        hTailGen_n_k4.SetBinContent( bin , 100 * pow( k , 4 ) * hGen_n.GetBinContent(bin) / hGen_n.Integral() )
        hTailWtdB_n_k4.SetBinContent( bin , 100 * pow( k , 4 ) * hWtdB_n.GetBinContent(bin) / hWtdB_n.Integral() )
        hTailWtdF_n_k4.SetBinContent( bin , 100 * pow( k , 4 ) * hWtdF_n.GetBinContent(bin) / hWtdF_n.Integral() )
    cMomentumDist.cd(3)
    hTailWtdB_n.Draw()
    hTailGen_n.Draw("same")
    hTailWtdF_n.Draw("same")
    anaB.AddLegend("n momentum dist."
                   ,hTailGen_n,"generated (%.1f%% ~ %d+/-%d)"%(TailGen_n,int(NTailGen_n),int(ROOT.sqrt(NTailGen_n)))
                   ,hTailWtdB_n,"Backwards (%.1f%% ~ %d+/-%d)"%(TailWtdB_n,int(NTailWtdB_n),int(ROOT.sqrt(NTailWtdB_n)))
                   ,hTailWtdF_n,"Backwards (%.1f%% ~ %d+/-%d)"%(TailWtdF_n,int(NTailWtdF_n),int(ROOT.sqrt(NTailWtdF_n)))
                   ,"F");

    cMomentumDist.cd(4)
    hTailWtdB_n_k4.Draw()
    hTailGen_n_k4.Draw("same")
    hTailWtdF_n_k4.Draw("same")
 
   
    # neutron before reaction
    cMomentumDist.cd(5)
    hFBasymWtd_n = ROOT.TH1F("hFBasymWtd_n","",Nbins,0,0.7)
    plot.SetFrame(hFBasymWtd_n ,"struck n f/b asym.","n momentum [GeV/c]","Backward n / Forwad n [%]",colorAsym,colorAsym)
    hFBasymGen_n = ROOT.TH1F("hFBasymGen_n","",Nbins,0,0.7)
    for bin in range(hGen_n.GetXaxis().FindBin(kF),Nbins):
        print "bin=%d, hWtdB_n.GetBinContent(bin) =%f hWtdF_n.GetBinContent(bin)=%f "%(bin,hWtdB_n.GetBinContent(bin),hWtdF_n.GetBinContent(bin) )
        hFBasymWtd_n.SetBinContent( bin , 100 * hWtdB_n.GetBinContent(bin) / hWtdF_n.GetBinContent(bin) )
        hFBasymGen_n.SetBinContent( bin , 100 )
    hFBasymWtd_n.Draw()
#    frac = float(hFBasymWtd_n.Integral())/hFBasymGen_n.Integral()
#    plot.Text(0.1 , hFBasymWtd_n.GetMaximum()
#              , "above k_{F} %.1f%% ~ %d+/-%d excess"%(100.*(frac-1),int((frac-1)*NTailGen_n),int(ROOT.sqrt((frac-1)*NTailGen_n))))
#
#    # recoil proton
#    cMomentumDist.cd(6)
#    hFBasymWtd_p = ROOT.TH1F("hFBasymWtd_p","",Nbins,0,0.7)
#    plot.SetFrame(hFBasymWtd_p ,"recoil p f/b asym.","p(rec) momentum [GeV/c]","Backward p(rec) / Forwad p(rec) [%]",colorAsym,colorAsym)
#    hFBasymGen_p = ROOT.TH1F("hFBasymGen_p","",Nbins,0,0.7)
#    for bin in range(hGen_n.GetXaxis().FindBin(kF),Nbins):
#        hFBasymWtd_p.SetBinContent( bin , 100 * hWtdB_p.GetBinContent(bin) / hWtdF_p.GetBinContent(bin) )
#        hFBasymGen_p.SetBinContent( bin , 100 )
#    hFBasymWtd_p.Draw()
#    frac = float(hFBasymWtd_p.Integral())/hFBasymGen_p.Integral()
#    plot.Text(0.1 , hFBasymWtd_p.GetMaximum()
#              , "above k_{F} %.1f%% ~ %d+/-%d excess"%(100.*(frac-1),int((frac-1)*NTailGen_p),int(ROOT.sqrt((frac-1)*NTailGen_p))))
    cMomentumDist.Update()
    wait()
    cMomentumDist.SaveAs("~/Desktop/bin%d.pdf"%PpPmuBin)









# make a histogram of shape variables that use the sphericity tensor

import ROOT

import gc
gc.disable()

print("a")
ROOT.gROOT.SetBatch(1)#don't show graphics
#ROOT.ROOT.EnableImplicitMT()#enables multi-threading
#ROOT.ROOT. DisableImplicitMT()

print("b")
floc = "/store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/"
fnames = [
    "PrivateSamples.SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"]
    #"PrivateSamples.SUEP_2018_mMed-400_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    #"PrivateSamples.SUEP_2018_mMed-750_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    #"PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
#]
print("c")
dFrames = {}
filteredFrames = {}
cHists = {}
dHists = {}
models = {}

trackPtCut=1

for fname in fnames:
    fullname = "root://cmsxrootd.fnal.gov/"+floc+fname
    tname = "TreeMaker2/PreSelection"
    mass = fname.split("_")[2]
    dFrames[mass] = ROOT.ROOT.RDataFrame(tname, fullname)
    print("d1")
    print (dFrames[mass].Count().GetValue())
    # it's more efficient to define ntracks before the loop, right?
    #find the HT the detector "sees" so that we can cut on that for l1 trigger:
    filteredFrames[mass]=dFrames[mass].Define("nTracks", "Tracks.size()") \
        .Define("CutHT", "double cutht=0; for (int i=0; i<Jets.size(); i++) if (Jets[i].Pt()>30 and abs(Jets[i].eta())<2.4) cutht+=Jets[i].Pt(); return cutht") \
        .Filter("CutHT>500") \
        .Define("Momenta", "vector<vector<double>> p; for (int i=0; i<nTracks; i++) {p[i].push_back(Tracks[i].x()); p[i].push_back(Tracks[i].y()); p[i].push_back(Tracks[i].z());} return p;") \
        .Define("Denominator", "double denom=0; for (int i=0; i<nTracks; i++) denom += sqrt(Tracks[i].Mag2()); return denom;") \
        .Define("SphericityTensor", "TMatrixDSym s(3,3); TArrayD array(9); for (int i=0; i<9; i++) array[i]=0; s.SetMatrixArray(array.GetArray()); for (int i=0; i<nTracks; i++) { for (int j=0; j<3; j++) {for (int k=0; k<3; k++) {s[j][k]+= (Momenta[i][j]*Momenta[i][k]/(sqrt(Tracks[i].Mag2())*Denominator));}}} return s;") \
        .Define("EigenVals", "TMatrixDSymEigen eigen(SphericityTensor); return eigen.GetEigenValues();") \
        .Define("Val", "return EigenVals[0]") \
        .Define("C", "return (3*(EigenVals[0]*EigenVals[1]+EigenVals[0]*EigenVals[2]+EigenVals[1]*EigenVals[2]));") \
        .Define("D","return 27*EigenVals[0]*EigenVals[1]*EigenVals[2];")
    #print(filteredFrames[mass].Count().GetValue())
    print("d4")
    models[mass+"C"] = ROOT.RDF.TH1DModel("C"+mass, mass, 50, 0., 1.)
    #cHists[mass] = filteredFrames[mass].Histo1D(models[mass+"C"], "C").Clone("cloneC"+mass)
    models[mass + "D"] = ROOT.RDF.TH1DModel("D" + mass, mass, 50, 0., 1.)
    #dHists[mass] = filteredFrames[mass].Histo1D(models[mass + "D"], "D").Clone("cloneD"+mass)
    cHists[mass] =filteredFrames[mass].Histo1D(models[mass + "C"], "Val").Clone("cloneC"+mass)
    print("d5")
print("e")


print("f")
print(cHists.keys())

#copied from plotHelpers:
can = ROOT.TCanvas("canName", "canTitle")
print("g")
file = ROOT.TFile('hists', 'RECREATE')
print("g1")
cHists["mMed-125"].Write()
print("g2")
'''
i = 0
more_colors = [ROOT.kAzure+2, ROOT.kGreen+2, ROOT.kPink+4, ROOT.kOrange+10, ROOT.kOrange, ROOT.kSpring+7]
print("h")
for mass in cHists.keys():
    print("i")
    cHists[mass].SetLineColor(more_colors[i])
    print("j")
    cHists[mass].SetMarkerColor(more_colors[i])
    print("k")
    i += 1
    print("l")
    i %= len(more_colors)'''

cHists["mMed-125"].SetLineColor(2)
#cHists["mMed-400"].SetLineColor(3)
#cHists["mMed-750"].SetLineColor(4)
#cHists["mMed-1000"].SetLineColor(6)
print("m")

#cHists["mMed-125"].SetMinimum(0)
#cHists["mMed-125"].SetMaximum(1)
cHists["mMed-125"].Draw("hist")
#cHists["mMed-400"].Draw("same")
#cHists["mMed-750"].Draw("same")
#cHists["mMed-1000"].Draw("same")

print("n")

leg = ROOT.TLegend(.66, .64, .8, .88)
print("o")
for mass in cHists.keys():
    leg.AddEntry(cHists[mass], mass, "l")
leg.Draw()
print("p")
can.SaveAs("shapesPlot.pdf")
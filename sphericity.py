# make a histogram of sphericity

import ROOT

print("a")
ROOT.gROOT.SetBatch(1)#don't show graphics
ROOT.ROOT.EnableImplicitMT()#enables multi-threading

print("b")
floc = "/store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/"
fnames = [
    "PrivateSamples.SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-400_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-750_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
]
print("c")
dFrames = {}
cHists = {}
dHists = {}

for fname in fnames:
    fullname = "root://cmsxrootd.fnal.gov/"+floc+fname
    tname = "TreeMaker2/PreSelection"
    mass = fname.split("_")[2]
    dFrames[mass] = ROOT.ROOT.RDataFrame(tname, fullname)
print("d")
#currently not used:
#weights={"mMed-125": 34.8, "mMed-400":5.9, "mMed-750": 0.5, "mMed-1000": 0.17}
#find a way to use this:

trackPtCut=1

for mass in dFrames:
    print("d1")
    # it's more efficient to define ntracks before the loop, right?
    print("d2")
    #find the HT the detector "sees" so that we can cut on that for l1 trigger:
    print("d3")
    filteredFrame=dFrames[mass].Define("nTracks", "Tracks.size()")
    print("d3a")
    filteredFrame=filteredFrame.Define("cutHT", "double cutht=0; for (int i=0; i<Jets.size(); i++) if (Jets[i].Pt()>30 and abs(Jets[i].eta())<2.4) cutht+=Jets[i].Pt(); return cutht")
    print("d3b")
    filteredFrame=filteredFrame.Filter("cutHT>500")
    print("d3c")
    #filteredFrame = filteredFrame.Define("trackMagnitudes", "vector<double> trackmags; for (int i=0; i<nTracks; i++)  trackmags.push_back(sqrt(Tracks[i].Mag2())); return trackmags;")
    filteredFrame = filteredFrame.Define("momenta", "vector<vector<double>> p; for (int i=0; i<nTracks; i++) {p[i].push_back(Tracks[i].x()); p[i].push_back(Tracks[i].y()); p[i].push_back(Tracks[i].z());} return p;")
    filteredFrame = filteredFrame.Define("denominator", "double denom=0; for (int i=0; i<nTracks; i++) denom += sqrt(Tracks[i].Mag2()); return denom;")
    filteredFrame = filteredFrame.Define("sphericityTensor", "TMatrixDSym s(3,3); TArrayD array(9); for (int i=0; i<9; i++) array[i]=0; s.SetMatrixArray(array.GetArray()); for (int i=0; i<nTracks; i++) { for (int j=0; j<3; j++) {for (int k=0; k<3; k++) {s[j][k]+= (momenta[i][j]*momenta[i][k]/(sqrt(Tracks[i].Mag2())*denominator));}}} return s;")
    filteredFrame = filteredFrame.Define("eigenVals", "TMatrixDSymEigen eigen(sphericityTensor); return eigen.GetEigenValues();")
    filteredFrame = filteredFrame.Define("C", "return 3*(eigenVals[0]*eigenVals[1]+eigenVals[0]*eigenVals[2]+eigenVals[1]*eigenVals[2]);")
    filteredFrame = filteredFrame.Define("D","return 27*eigenVals[0]*eigenVals[1]*eigenVals[2];")
    print("d4")
    cHists[mass] = filteredFrame.Histo1D(("C"+mass,mass, 50, 0., 1.), "C")
    dHists[mass] = filteredFrame.Histo1D(("D" + mass, mass, 50, 0., 1.), "D")
    print("d5")
print("e")
#can I put this before the loop so I can combine the loops?

print("f")

#copied from plotHelpers:
can = ROOT.TCanvas("canName", "canTitle")

i = 0
more_colors = [ROOT.kAzure+2, ROOT.kGreen+2, ROOT.kPink+4, ROOT.kOrange+10, ROOT.kOrange, ROOT.kSpring+7]
for mass in cHists.keys():
    cHists[mass].SetLineColor(more_colors[i])
    cHists[mass].SetMarkerColor(more_colors[i])
    i += 1
    i %= len(more_colors)

cHists["mMed-125"].SetMinimum(0)
cHists["mMed-125"].SetMaximum(1)
cHists["mMed-125"].Draw("hist")
cHists["mMed-400"].Draw("same")
cHists["mMed-750"].Draw("same")
cHists["mMed-1000"].Draw("same")

leg = ROOT.TLegend(.66, .64, .8, .88)
for mass in cHists.keys():
    leg.AddEntry(cHists[mass], mass, "l")
leg.Draw()

can.SaveAs("sphericity.pdf")
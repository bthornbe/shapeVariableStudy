# a minimal reproducer for the error in sphericity.py

import ROOT

fullname = "root://cmsxrootd.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/PrivateSamples.SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"

trackPtCut=1

tname = "TreeMaker2/PreSelection"

dFrame = ROOT.ROOT.RDataFrame(tname, fullname)

dFrame2 = dFrame.Define("nTracks", "Tracks.size()")
dFrame3 = dFrame2.Define("CutHT",
            "double cutht=0; for (int i=0; i<Jets.size(); i++) if (Jets[i].Pt()>30 and abs(Jets[i].eta())<2.4) cutht+=Jets[i].Pt(); return cutht")
dFrame4 = dFrame3.Filter("CutHT>500")
dFrame5 = dFrame4.Define("Momenta",
            "vector<vector<double>> p; for (int i=0; i<nTracks; i++) {p[i].push_back(Tracks[i].x()); p[i].push_back(Tracks[i].y()); p[i].push_back(Tracks[i].z());} return p;") \
    dFrame6 = dFrame5.Define("Val", "return (Momenta[0][0])")

model = ROOT.RDF.TH1DModel("Val", "Val", 50, 0., 1.)
cHist = dFrame6.Histo1D(model, "Val")

can = ROOT.TCanvas("canName", "canTitle")
file = ROOT.TFile('reproducerHists', 'RECREATE')
cHist.Write()
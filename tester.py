import ROOT

fullname = "root://cmsxrootd.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/PrivateSamples.SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"

tname = "TreeMaker2/PreSelection"

dFrame = ROOT.ROOT.RDataFrame(tname, fullname).Define("nTracks", "Tracks.size()")

model = ROOT.RDF.TH1DModel("nTracks", "nTracks", 50, 0., 1.)
hist = dFrame.Histo1D(model, "nTracks")

can = ROOT.TCanvas("canName", "canTitle")

hist.Draw()

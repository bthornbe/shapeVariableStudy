# a minimal reproducer for the error in sphericity.py

import ROOT

fullname = "root://cmsxrootd.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/PrivateSamples.SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"

trackPtCut=1

tname = "TreeMaker2/PreSelection"

dFrame = ROOT.ROOT.RDataFrame(tname, fullname)

filteredFrame= dFrame.Define("nTracks", "Tracks.size()") \
    .Define("CutHT",
            "double cutht=0; for (int i=0; i<Jets.size(); i++) if (Jets[i].Pt()>30 and abs(Jets[i].eta())<2.4) cutht+=Jets[i].Pt(); return cutht") \
    .Filter("CutHT>500") \
    .Define("Momenta",
            "vector<vector<double>> p; for (int i=0; i<nTracks; i++) {p[i].push_back(Tracks[i].x()); p[i].push_back(Tracks[i].y()); p[i].push_back(Tracks[i].z());} return p;") \
    .Define("Denominator",
            "double denom=0; for (int i=0; i<nTracks; i++) denom += sqrt(Tracks[i].Mag2()); return denom;") \
    .Define("SphericityTensor",
            "TMatrixDSym s(3,3); TArrayD array(9); for (int i=0; i<9; i++) array[i]=0; s.SetMatrixArray(array.GetArray()); for (int i=0; i<nTracks; i++) { for (int j=0; j<3; j++) {for (int k=0; k<3; k++) {s[j][k]+= (Momenta[i][j]*Momenta[i][k]/(sqrt(Tracks[i].Mag2())*Denominator));}}} return s;") \
    .Define("EigenVals", "TMatrixDSymEigen eigen(SphericityTensor); return eigen.GetEigenValues();") \
    .Define("C", "return (3*(EigenVals[0]*EigenVals[1]+EigenVals[0]*EigenVals[2]+EigenVals[1]*EigenVals[2]));") \

model = ROOT.RDF.TH1DModel("C", "C", 50, 0., 1.)
cHist = filteredFrame.Histo1D(model, "C").Clone("cloneC")

can = ROOT.TCanvas("canName", "canTitle")
file = ROOT.TFile('hists', 'RECREATE')
cHist.Write()
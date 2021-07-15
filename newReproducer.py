import ROOT

fullname = "root://cmsxrootd.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/PrivateSamples.SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"

tname = "TreeMaker2/PreSelection"

mystring= \
'''
TArrayD array(9);   
for (int i=0; i<nTracks; i++) { 
    for (int j=0; j<3; j++) {
        for (int k=0; k<3; k++) {
            array[j+3*k] += Momenta.at(i).at(j)*Momenta.at(i).at(k)/(sqrt(Tracks[i].Mag2())*Denominator);
        }
    }
} 
TMatrixDSym s(3, array.GetArray());
return s;
'''
dFrame= ROOT.ROOT.RDataFrame(tname, fullname).Define("nTracks", "Tracks.size()") \
        .Filter("nTracks > 0") \
        .Define("CutHT", "double cutht=0; for (int i=0; i<Jets.size(); i++) if (Jets[i].Pt()>30 and abs(Jets[i].eta())<2.4) cutht+=Jets[i].Pt(); return cutht") \
        .Filter("CutHT>500") \
        .Define("Momenta", "vector<vector<double>> p; for (int i=0; i<nTracks; i++) {p.emplace_back(); p[i].push_back(Tracks[i].x()); p[i].push_back(Tracks[i].y()); p[i].push_back(Tracks[i].z());} return p;") \
        .Define("Denominator", "double denom=0; for (int i=0; i<nTracks; i++) denom += sqrt(Tracks[i].Mag2()); return denom;") \
        .Define("SphericityTensor", mystring) \
        .Define("Val", "return (SphericityTensor[0][0])") \


model = ROOT.RDF.TH1DModel("Val", "Val", 50, -1000., 1000.)
hist = dFrame.Histo1D(model, "Val")

can = ROOT.TCanvas("canName", "canTitle")
file = ROOT.TFile('reproducerHists', 'RECREATE')
hist.Write()

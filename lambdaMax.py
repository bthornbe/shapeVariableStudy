# make a histogram of the maximum eigenvalue of the sphericity tensor with r = 1
#currently not working
import ROOT
import cms_style

cms_style.setTDRStyle()

ROOT.gROOT.SetBatch(1)  # don't show graphics

floc = "/store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/"
fnames = [
    "PrivateSamples.SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-400_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-750_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"]

dFrames = {}
filteredFrames = {}
hists = {}
models = {}

trackPtCut = 1

# creates the sphericity tensor with r=1, as used for C, D, and LamdaMax:
tensorString1 = \
    '''
    vector<vector<double>> s{{0,0,0},{0,0,0},{0,0,0}};    
    for (int i=0; i<nPassingTracks; i++) 
    { 
        for (int j=0; j<3; j++) 
        {
            for (int k=0; k<3; k++) 
            {
                s.at(j).at(k) += Momenta.at(i).at(j)*Momenta.at(i).at(k)/(sqrt(PassingTracks[i].Mag2())*Denominator);
            }
        }
    } 
    return s;
    '''

eigenString = \
    '''
    double a[9];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            a[i+3*j] = SphericityTensor1.at(i).at(j);
        }
    }
    TMatrixDSym s(3, a);
    TMatrixDSymEigen eigen(s); 
    TVectorD eigenVals = eigen.GetEigenValues();
    vector <double> eigenVals2;
    for (int i = 0; i < 3; i++)
    {
            eigenVals2.push_back(eigenVals[i]);
    }
    return eigenVals2;
    '''

trackPtCut = 1
passTrackString = \
    '''
    auto passTracks = Tracks; 
    passTracks.clear(); 
    for (int i = 0; i < nTracks; i ++) 
    {
        if (Tracks[i].Perp2() > ''' + str(trackPtCut ** 2) + ''' && abs(Tracks[i].Eta()) < 2.5 && Tracks_fromPV0[i] >=2  && Tracks_matchedToPFCandidate[i]) 
        passTracks.push_back(Tracks[i]);
} 
return passTracks;
'''

for fname in fnames:
    fullname = "root://cmsxrootd.fnal.gov/" + floc + fname
    tname = "TreeMaker2/PreSelection"
    mass = fname.split("_")[2]
    dFrames[mass] = ROOT.ROOT.RDataFrame(tname, fullname)

    filteredFrames[mass] = dFrames[mass].Define("nTracks", "Tracks.size()") \
        .Filter("nTracks > 0") \
        .Define("CutHT",
                "double cutht=0; for (int i=0; i<Jets.size(); i++) if (Jets[i].Pt()>30 and abs(Jets[i].eta())<2.4) cutht+=Jets[i].Pt(); return cutht") \
        .Filter("CutHT>500") \
        .Define("PassingTracks", passTrackString) \
        .Define("nPassingTracks", "PassingTracks.size()") \
        .Define("Momenta",
                "vector<vector<double>> p; for (int i=0; i<nPassingTracks; i++) {p.emplace_back(); p[i].push_back(PassingTracks[i].x()); p[i].push_back(PassingTracks[i].y()); p[i].push_back(PassingTracks[i].z());} return p;") \
        .Define("Denominator",
                "double denom=0; for (int i=0; i<nPassingTracks; i++) denom += sqrt(PassingTracks[i].Mag2()); return denom;") \
        .Define("SphericityTensor1", tensorString1) \
        .Define("EigenVals", eigenString) \
        .Define("LambdaMax", "return EigenVals[0];")
    models[mass] = ROOT.RDF.TH1DModel(mass, mass, 50, 0., 1.)
    hists[mass] = filteredFrames[mass].Histo1D(models[mass], "LambdaMax").Clone("cloneL" + mass)

can = ROOT.TCanvas("canName", "canTitle")

hists["mMed-125"].SetLineColor(2)
hists["mMed-400"].SetLineColor(3)
hists["mMed-750"].SetLineColor(4)
hists["mMed-1000"].SetLineColor(6)

hists["mMed-125"].SetMaximum(600)
hists["mMed-125"].Draw("hist")
hists["mMed-400"].Draw("same")
hists["mMed-750"].Draw("same")
hists["mMed-1000"].Draw("same")

leg = ROOT.TLegend(.66, .64, .8, .88)

for mass in hists.keys():
    leg.AddEntry(hists[mass], mass, "l")
leg.Draw()

can.SaveAs("lPlot.pdf")

# make a histogram of sphericity

import ROOT
import cms_style
cms_style.setTDRStyle()


ROOT.gROOT.SetBatch(1) # don't show graphics

# simulated signal:
floc = "/store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/"
fnames = [
    "PrivateSamples.SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-400_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-750_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root"]

# background:
bfloc = "/eos/user/t/tholmes/public/SUEPs/"
bfnames = [
    "Autumn18.QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root"
]

# sphericity tensor for r=2, used for sphericity:
tensorString2 = \
    '''
    vector<vector<double>> s{{0,0,0},{0,0,0},{0,0,0}};    
    for (int i=0; i<nPassingTracks; i++) 
    { 
        for (int j=0; j<3; j++) 
        {
            for (int k=0; k<3; k++) 
            {
                s.at(j).at(k) += Momenta.at(i).at(j)*Momenta.at(i).at(k)/Denominator2;
            }
        }
    } 
    return s;
    '''

eigenString2 = \
'''
double a[9];
for (int i = 0; i < 3; i++)
{
    for (int j = 0; j < 3; j++)
    {
        a[i+3*j] = SphericityTensor2.at(i).at(j);
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
    if (Tracks[i].Perp2() > ''' + str(trackPtCut**2) + ''' && abs(Tracks[i].Eta()) < 2.5 && Tracks_fromPV0[i] >=2  && Tracks_matchedToPFCandidate[i]) 
        passTracks.push_back(Tracks[i]);
} 
return passTracks;
'''

lum=135*1000
xSecArr=[311900,29070,5962,1207,119.9,25.24]
signalXSecs = [34.8, 5.9, 0.5, 0.17]
xSecs = {}
dFrames = {}
filteredFrames = {}
hists = {}
models = {}
trackPtCut = 1
weights = {}

for i, bfname in enumerate(bfnames):
    fullname = "root://cmseos.fnal.gov/ "+bfloc + bfname
    tname = "TreeMaker2/PreSelection"
    range = bfname.split("_")[1]
    dFrames[range] = ROOT.ROOT.RDataFrame("TreeMaker2/PreSelection", bfloc + bfname)
    xSecs[range] = xSecArr[i]
    entries = dFrames[range].Count().Getvalue()
    weights[range] = xSecs[range]*lum/entries

for fname in fnames:
    fullname = "root://cmsxrootd.fnal.gov/"+floc+fname
    tname = "TreeMaker2/PreSelection"
    mass = fname.split("_")[2]
    dFrames[mass] = ROOT.ROOT.RDataFrame(tname, fullname)

for key in dFrames.getkeys:
    # it's more efficient to define ntracks before the loop, right?
    #find the HT the detector "sees" so that we can cut on that for l1 trigger:
    filteredFrames[key]=dFrames[key].Define("nTracks", "Tracks.size()") \
        .Filter("nTracks > 0") \
        .Define("CutHT",
                "double cutht=0; for (int i=0; i<Jets.size(); i++) if (Jets[i].Pt()>30 and abs(Jets[i].eta())<2.4) cutht+=Jets[i].Pt(); return cutht") \
        .Filter("CutHT>500") \
        .Define("PassingTracks", passTrackString) \
        .Define("nPassingTracks", "PassingTracks.size()") \
        .Define("Momenta",
                "vector<vector<double>> p; for (int i=0; i<nPassingTracks; i++) {p.emplace_back(); p[i].push_back(PassingTracks[i].x()); p[i].push_back(PassingTracks[i].y()); p[i].push_back(PassingTracks[i].z());} return p;") \
        .Define("Denominator2", "double denom=0; for (int i=0; i<nPassingTracks; i++) denom += PassingTracks[i].Mag2(); return denom;") \
        .Define("SphericityTensor2", tensorString2) \
        .Define("EigenVals2",eigenString2) \
        .Define("Sphericity", "return (EigenVals2 [1] + EigenVals2 [2])*3/2")# I'm fairly sure GetEigenValues sorts the output from highest to lowest

    models[key + "S"] = ROOT.RDF.TH1DModel("S" + key, key, 50, 0., 1.)
    hists[key] = filteredFrames[key].Histo1D(models[key + "S"], "Sphericity").Clone("cloneS" + key)

can = ROOT.TCanvas("canName", "canTitle")

hists["mMed-125"].SetLineColor(2)
hists["mMed-400"].SetLineColor(3)
hists["mMed-750"].SetLineColor(4)
hists["mMed-1000"].SetLineColor(6)

hists["mMed-125"].SetMaximum(1300)
hists["mMed-125"].Draw("hist")
hists["mMed-400"].Draw("same")
hists["mMed-750"].Draw("same")
hists["mMed-1000"].Draw("same")

leg = ROOT.TLegend(.66, .64, .8, .88)

for key in hists.keys():
    leg.AddEntry(hists[key], key, "l")
leg.Draw()

can.SaveAs("sphericityPlot.pdf")
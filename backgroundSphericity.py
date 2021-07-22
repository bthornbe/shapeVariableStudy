# make a histogram of sphericity

import ROOT
import cms_style
cms_style.setTDRStyle()


ROOT.gROOT.SetBatch(1) # don't show graphics

floc = "/eos/user/t/tholmes/public/SUEPs/"
fnames = [
    "Autumn18.QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root",
    "Autumn18.QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_RA2AnalysisTree.root"
]

lum=135*1000
xsec=[311900,29070,5962,1207,119.9,25.24]
xsec2 = {}

dFrames = {}
filteredFrames = {}
bHist = {} #background histogram
bHist["events"]=ROOT.TH1F("Data", "QCD Histogram", nbins, 0, maxX)
models = {}

trees = {}
for k, fname in enumerate(fnames):
    fullname = floc+fname
    tname = "TreeMaker2/PreSelection"
    Range = fname.split("_")[1]
    trees[Range] = ROOT.ROOT.RDataFrame(tname, fullname)
    xsec2[Range] = xsec[k]

trackPtCut = 1

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

for fname in fnames:
    fullname = "root://cmsxrootd.fnal.gov/"+floc+fname #change for background?

    tname = "TreeMaker2/PreSelection"
    mass = fname.split("_")[2]
    dFrames[mass] = ROOT.ROOT.RDataFrame(tname, fullname)

    # it's more efficient to define ntracks before the loop, right?
    #find the HT the detector "sees" so that we can cut on that for l1 trigger:
    filteredFrames[mass]=dFrames[mass].Define("nTracks", "Tracks.size()") \
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
        .Define("Sphericity", "return (EigenVals2 [1] + EigenVals2 [2])*3/2")# the nTracks cut is probably unnecessary

    models[mass + "S"] = ROOT.RDF.TH1DModel("S" + mass, mass, 50, 0., 1.)
    sHists[mass] = filteredFrames[mass].Histo1D(models[mass + "S"], "Sphericity").Clone("cloneS" + mass)
    bHist["events"].Add()

can = ROOT.TCanvas("canName", "canTitle")

sHists["mMed-125"].SetLineColor(2)
sHists["mMed-400"].SetLineColor(3)
sHists["mMed-750"].SetLineColor(4)
sHists["mMed-1000"].SetLineColor(6)

sHists["mMed-125"].SetMaximum(1300)
sHists["mMed-125"].Draw("hist")
sHists["mMed-400"].Draw("same")
sHists["mMed-750"].Draw("same")
sHists["mMed-1000"].Draw("same")

leg = ROOT.TLegend(.66, .64, .8, .88)

for mass in sHists.keys():
    leg.AddEntry(sHists[mass], mass, "l")
leg.Draw()

can.SaveAs("sphericityPlot.pdf")
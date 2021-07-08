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
hists = {}

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
    #filteredFrame = filteredFrame.Define("trackMagnitudes", "vector<double> trackmags; for (int i=0; i<nTracks; i++)  trackmags.push_back(sqrt(Tracks[i].Perp2())); return trackmags;")
    filteredFrame = filteredFrame.Define("momenta", "double p[nTracks][3]; for (int i=0; i<nTracks; i++) {p[i][0]=Tracks[i].x(); p[i][1]=Tracks[i].y(); p[i][2]=Tracks[i].z();} return p;")
    filteredFrame = filteredFrame.Define("denominator", " double denom=0; for (int i=0; i<nTracks; i++) denom += sqrt(Tracksi[i].Mag2); return denom;")
    filteredFrame = filteredFrame.Define("sphericityTensor", "double s[3][3]={{0,0,0},{0,0,0},{0,0,0}}; for (int i=0; i<nTracks; i++) { for (int j=0; j<3; j++) {for (int k=0; k<3; k++) {s[j][k]+= (momenta[i][j]*momenta[i][k]/(sqrt(Tracksi[i].Mag2)*denominator)}}}; return s;")
    # perp2 is pt squared
    filteredFrame=filteredFrame.Define("filteredTracks", " int ftracks=0; for (int i=0; i<nTracks; i++) {if (TrackPtSquared[i] >" + str(trackPtCut**2) + " && abs(Tracks[i].Eta())<2.5 && Tracks_fromPV0[i]>=2 && Tracks_matchedToPFCandidate[i]) ftracks++;} return ftracks;")
    print("d4")
    hists[mass] = filteredFrame.Histo1D(("filteredTracks"+mass,mass, 50, 0., 500.), "filteredTracks")
    print("d5")
print("e")
#can I put this before the loop so I can combine the loops?

print("f")

#copied from plotHelpers:
can = ROOT.TCanvas("canName", "canTitle")

h_keys = hists.keys()
h_values = hists.values()

maxy = 1.5 * 1e4 * max(hists[mass].GetMaximum() for mass in hists)
miny = 1e-1

i = 0
more_colors = [ROOT.kAzure+2, ROOT.kGreen+2, ROOT.kPink+4, ROOT.kOrange+10, ROOT.kOrange, ROOT.kSpring+7]
for mass in hists:
    hists[mass].SetLineColor(more_colors[i])
    hists[mass].SetMarkerColor(more_colors[i])
    i += 1
    i %= len(more_colors)

hists["mMed-125"].SetMinimum(miny)
hists["mMed-125"].SetMaximum(maxy)
hists["mMed-125"].Draw("hist")
ROOT.gPad.SetLogy(1)
hists["mMed-400"].Draw("same")
hists["mMed-750"].Draw("same")
hists["mMed-1000"].Draw("same")

leg = ROOT.TLegend(.66, .64, .8, .88)
for mass in hists.keys():
    leg.AddEntry(hists[mass], mass, "l")
leg.Draw()

can.SaveAs("rtestplot.pdf")
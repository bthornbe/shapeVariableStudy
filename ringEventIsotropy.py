# creates histogram of the ring event isotropy

import ROOT
import numpy as np
from shapeMethods import from1dhisttoarray
from eventIsotropy.emdVar import _cdist_phicos, emd_Calc
from eventIsotropy.cylGen import ringGen

ROOT.gROOT.SetBatch(1)# don't show graphics
ROOT.ROOT.EnableImplicitMT()# enables multi-threading

floc = "/store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/"
fnames = [
    "PrivateSamples.SUEP_2018_mMed-125_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-400_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-750_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
    "PrivateSamples.SUEP_2018_mMed-1000_mDark-2_temp-2_decay-darkPho_13TeV-pythia8_n-100_0_RA2AnalysisTree.root",
]
# declaring dicts before the loop so we don't lose anything to memory cleaning:
dFrames = {}
# this will be a dictionary of 1d histograms representing the distribution of transverse energy:
hists = {}
# this will be the same data as hists but stored as ndarrays:
arrays = {}
# M will be the distance matrix of the data and uniform distributions, here using the (1 - cos(phi)) metric
M = {}

emdVals = {}
# number of bins:
num = input("number of bins in phi:")

# phi values of uniform(except it's discretized) distribution:
uniformRing = ringGen(num)

# non-normalized weights (normalized in emd_calc):
uniformRingPt=np.array([np.full(len(uniformRing[i]), 1.) for i in range(num)])

for fname in fnames:
    fullname = "root://cmsxrootd.fnal.gov/"+floc+fname
    tname = "TreeMaker2/PreSelection"
    mass = fname.split("_")[2]
    dFrames[mass] = ROOT.ROOT.RDataFrame(tname, fullname)

for mass in dFrames:
    #it's more efficient to define nTracks before the loop, right?
    dFrames[mass]=dFrames[mass].Define("nTracks", "Tracks.size()")
    # find the HT the detector "sees" so that we can cut on that for l1 trigger:
    dFrames[mass]=dFrames[mass].Define("cutHT", "double cutht=0; for (int i=0; i<nTracks; i++) if (Jets[i].Pt()>30 and abs(Jets[i].eta())<2.4) cutht+=Jets[i].Pt(); return cutht")
    filteredFrame=dFrames[mass].Filter("cutHT>500")
    #! describe this better:
    # by creating a histogram of phi weighted by pt, I bin the transverse energy:
    hists[mass] = filteredFrame.Histo1D((mass+"ring", mass, num, 0., 2*np.pi), "phi", "pt")
    arrays[mass] = from1dhisttoarray(hists[mass]) #! this method doesn't work yet
    M[mass] = _cdist_phicos(arrays[mass], uniformRing)
    emdVals[mass] = emd_Calc(arrays[mass], uniformRingPt, M, numItermax=100000000,log=True)

can = ROOT.TCanvas("canName", "canTitle")

i=0
#copied from plotHelpers:
more_colors = [ROOT.kAzure+2, ROOT.kGreen+2, ROOT.kPink+4, ROOT.kOrange+10, ROOT.kOrange, ROOT.kSpring+7]
for mass in hists:
    hists[mass].SetLineColor(more_colors[i])
    hists[mass].SetMarkerColor(more_colors[i])
    i += 1
    i %= len(more_colors)

hists["mMed-125"].Draw()
hists["mMed-400"].Draw("same")
hists["mMed-750"].Draw("same")
hists["mMed-1000"].Draw("same")

can.Draw()
can.SaveAs("rtestplot.pdf")
import ROOT
import cms_style
cms_style.setTDRStyle()

ROOT.gROOT.SetBatch(1) # don't show graphics
# this is a normal comment

#TODO this is a reminder for myself

#? this is a question I want to remember to ask

file = ROOT.TFile('sphericityHists.root', 'READ')

keys = file.GetListOfKeys()

print(keys)

lum=135*1000
xSecs=[311900,29070,5962,1207,119.9,25.24]
signalXSecs = [34.8, 5.9, 0.5, 0.17]
bcknum=10**5
signum=10**4

theme_colors2 = [ROOT.TColor.GetColor("#FFA900"),   # Sunny Orange
                ROOT.TColor.GetColor("#3629AC"),    # Dark Blue
                ROOT.TColor.GetColor("#2FC494"),    # Seafoam
                ROOT.TColor.GetColor("#F65866"),    # Pink
                ROOT.TColor.GetColor("#0E81C4"),    # Light Blue
                ]

hists = {}

for key in keys:
    k=key.GetName()
    histType = k.split("_")[0]
    if histType not in hists:
        hists[histType] = {"sig":{}, "bck":{}}
    if "mMed" in k:
        hists[histType]["sig"][k.split("_")[1]]=file.Get(k)
    else:
        hists[histType]["bck"][k.split("_")[1]]=file.Get(k)

sigs = sorted(hists[histType]["sig"].keys(), key=lambda x:int(x.split("-")[1]))
bcks = sorted(hists[histType]["bck"].keys(), key=lambda x:int(x.split("to")[0].strip("HT")))

for histType in hists:
    can = ROOT.TCanvas("can" + histType, "canTitle")
    if len(bcks) > 0:
        for i, b in enumerate(bcks):
            h = hists[histType]["bck"][b]
            h.Scale(lum * xSecs[i]/bcknum)
            if i > 0:
                hists[histType]["bck"][bcks[0]].Add(h)
        h_bck = hists[histType]["bck"][bcks[0]]
        h_bck.SetFillColor(ROOT.kGray)
        h_bck.SetMarkerSize(0)
        h_bck.SetMinimum(10**4)
        h_bck.SetMaximum(10**11)
        h_bck.GetXaxis().SetTitle(histType)
        h_bck.GetYaxis().SetTitle("Events")
        h_bck.SetTitle("QCD")
        h_bck.Draw("ehist")
    for i, b in enumerate(sigs):
        h = hists[histType]["sig"][b]
        h.Scale(lum * xSecs[i]/signum)
        h.SetMarkerSize(0)
        h.SetLineColor(theme_colors2[i])
        h.SetMinimum(10**4)
        h.SetMaximum(10**11)
        h.Draw("ehist same")
    can.BuildLegend(.76, .64, .9, .88).Draw()
    can.SetLogy()
    l = ROOT.TLatex()
    l.SetNDC()
    l.DrawLatex(.2,.9, "CMS Internal")
    l.DrawLatex(.2, .8, "%d fb^{-1}"%(lum/1000))
    can.SaveAs(histType + ".pdf")

file.Close()

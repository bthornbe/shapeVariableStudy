import ROOT

# this is a normal comment

#TODO this is a reminder for myself

#? this is a question I want to remember to ask

file = ROOT.TFile('sphericityHists.root', 'READ')

keys = file.GetListOfKeys()

print(keys)

hists = {}

for k in keys:
    histType = k.split("_")[0]
    if histType not in hists:
        hists[histType] = {"sig":{}, "bck":{}}
    if "mMed" in k:
        hists[histType]["sig"][k.split("_")[1]]=file.Get(k)
    else:
        hists[histType]["bck"][k.split("_")[1]]=file.Get(k)

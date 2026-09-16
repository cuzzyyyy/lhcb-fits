import ROOT

# wyciszenie logow roofita, zeby nie miec syfu w terminalu
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

file = ROOT.TFile.Open("Bs_DsKstar_magup.root", "READ")
tree = file.Get("LL_KKpi/DecayTree")

mass_var = "KS0_M"
cat_var = "KS0_BKGCAT"

mass_roo = ROOT.RooRealVar(mass_var, "mass [MeV/c2]", 400, 600)
bkgcat_roo = ROOT.RooRealVar(cat_var, "background category", 0, 150)

dataset = ROOT.RooDataSet(
    "data",
    "DecayTree Data",
    ROOT.RooArgSet(mass_roo, bkgcat_roo),
    ROOT.RooFit.Import(tree),
)

categories = {
    "0": ("signal", ROOT.kGreen),
    "10": ("quasi-signal", ROOT.kCyan),
    "20": ("fully reconstructed physical background", ROOT.kBlue),
    "30": ("mis-ID", ROOT.kOrange),
    "40": ("partially reconstructed physical background", ROOT.kBlue),
    "50": ("low-mass background", ROOT.kBlue),
    "60, 63": ("ghost particles - ghost and ghost clone", ROOT.kRed),
    ">63": ("other events", ROOT.kBlack),
}

canvas = ROOT.TCanvas("canvas", "background categories", 1200, 800)
canvas.Divide(4, 2)

frames = []

for index, (key_str, (label, color)) in enumerate(categories.items()):
    if ">" in key_str:
        value = key_str.replace(">", "").strip()
        cut = f"{cat_var} > {value}"
    elif "," in key_str:
        elements = key_str.split(",")
        cut = " || ".join([f"{cat_var} == {el.strip()}" for el in elements])
    else:
        cut = f"{cat_var} == {key_str}"

    reduced_data = dataset.reduce(ROOT.RooFit.Cut(cut))

    hist_name = f"hist_cat_{index}"
    hist = reduced_data.createHistogram(hist_name, mass_roo, ROOT.RooFit.Binning(100))

    hist.SetTitle(label)
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.SetStats(0)

    hist.GetYaxis().SetTitle("Events count")
    hist.GetXaxis().SetTitle("mass [MeV/c^{2}]")

    frames.append(
        hist
    )

    pad = canvas.cd(index + 1)
    pad.SetBottomMargin(0.15)
    pad.SetLeftMargin(0.15)

    hist.Draw("HIST")

canvas.Update()
canvas.SaveAs(f"{mass_var}_bkgcat_LL_roofit.png")
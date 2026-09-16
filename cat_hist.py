import ROOT

file = ROOT.TFile.Open("Bs_DsKstar_magup.root", "READ")
tree = file.Get("DD_KKpi/DecayTree")

mass_var = "B_M"
cat_var = "B_BKGCAT"

min_x = 5000
max_x = 5700
nbins = 100

categories = {
    0: ("signal", ROOT.kGreen),
    10: ("quasi-signal", ROOT.kCyan),
    30: ("mis-ID", ROOT.kOrange),
    "20, 40, 50": ("combined background - physical and low-mass", ROOT.kBlue),
    "60, 63": ("ghost particles - ghost and ghost clone", ROOT.kRed),
    ">63": ("other events", ROOT.kBlack),
}

canvas = ROOT.TCanvas("canvas", "background categories", 1200, 800)
canvas.Divide(3, 2)

histograms = []

for key, (label, color) in categories.items():
    key_str = str(key)

    if ">" in key_str:
        value = key_str.replace(">", "").strip()
        cut = f"{cat_var} > {value}"
    elif "," in key_str:
        elements = key_str.split(",")
        cut = " || ".join([f"{cat_var} == {el.strip()}" for el in elements])
    else:
        cut = f"{cat_var} == {key_str}"

    hist_name = f"h_{color}"

    hist = ROOT.TH1F(
        hist_name, f"{label};mass [MeV/c^2];events count", nbins, min_x, max_x
    )

    tree.Draw(f"{mass_var} >> {hist_name}", cut, "goff")

    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.SetStats(0)

    histograms.append(hist)

for index, hist in enumerate(histograms):
    pad = canvas.cd(index + 1)

    pad.SetBottomMargin(0.15)
    pad.SetLeftMargin(0.15)

    # pad.SetLogy()

    hist.Draw("HIST")

canvas.Update()
canvas.SaveAs(f"{mass_var}_bkgcat_DD.png")
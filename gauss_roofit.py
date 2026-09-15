import ROOT

# wyciszenie logow roofita, zeby nie miec syfu w terminalu
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

file = ROOT.TFile.Open("Bs_DsKstar_magup.root", "READ")

# 1. Definicja ścieżki i automatyczne wyciągnięcie "DD" lub "LL"
tree_path = "DD_KKpi/DecayTree"
track_type = tree_path.split("_")[0]

tree = file.Get(tree_path)

mass_var = "KS0_M"
cat_var = "KS0_BKGCAT"
low_limit = 440
upp_limit = 560


mass_roo = ROOT.RooRealVar(mass_var, "mass [MeV/c2]", low_limit, upp_limit)
bkgcat_roo = ROOT.RooRealVar(cat_var, "background category", 0, 150)

dataset = ROOT.RooDataSet(
    "data",
    "data",
    ROOT.RooArgSet(mass_roo, bkgcat_roo),
    ROOT.RooFit.Import(tree),
)

# odcinanie czesci danych ze wzgledu na kategorie
reduced_data = dataset.reduce(ROOT.RooFit.Cut(f"{cat_var} < 30"))

mean = ROOT.RooRealVar("mean", "mean", 500, 480, 520)
sigma = ROOT.RooRealVar("sigma", "sigma", 5, 0, 10)
gauss = ROOT.RooGaussian("gauss", "gauss pdf", mass_roo, mean, sigma)

# fitowanie
gauss.fitTo(reduced_data, ROOT.RooFit.PrintLevel(-1))

frame = mass_roo.frame(ROOT.RooFit.Title("gaussian fit for data cut (BKG_CAT < 30)"))

reduced_data.plotOn(frame)
gauss.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))

# pull jako obiekt TGraph
pull_graph = frame.pullHist()

# konwersja punktow pulla do klasycznego histogramu
pull_hist = ROOT.TH1F("pull_hist", ";mass [MeV/c^{2}];Pull", 100, low_limit, upp_limit)
x_vals = pull_graph.GetX()
y_vals = pull_graph.GetY()
for i in range(pull_graph.GetN()):
    bin_idx = pull_hist.FindBin(x_vals[i])
    pull_hist.SetBinContent(bin_idx, y_vals[i])

pull_hist.SetFillColor(ROOT.kAzure - 9)
pull_hist.SetLineColor(ROOT.kBlue + 2)
pull_hist.SetStats(0)
pull_hist.GetYaxis().SetRangeUser(-5, 5)
pull_hist.GetYaxis().SetTitleSize(0.08)
pull_hist.GetYaxis().SetTitleOffset(0.5)

canvas = ROOT.TCanvas("canvas", "Fit and Pull", 800, 800)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.SetBottomMargin(0.02)
pad1.Draw()

pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.3)
pad2.Draw()

pad1.SetBottomMargin(0.12)  # robi odstęp na dole górnego panelu
pad2.SetTopMargin(0.05)  # robi odstęp na górze dolnego panelu

pad1.cd()
frame.Draw()

pad2.cd()
pull_hist.Draw("HIST")

# linia zerowa na panelu pulla dla czytelnosci
line = ROOT.TLine(low_limit, 0, upp_limit, 0)
line.SetLineColor(ROOT.kBlack)
line.SetLineStyle(2)
line.Draw("SAME")

# automatyczny zapis wykorzystujący track_type
canvas.SaveAs(f"{mass_var}_gaussian_fit_pull_{track_type}_30.png")

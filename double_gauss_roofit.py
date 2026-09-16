import ROOT

# wyciszenie logow roofita, zeby nie miec syfu w terminalu
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

file = ROOT.TFile.Open("Bs_DsKstar_magup.root", "READ")

# definicja sciezki i automatyczne wyciagniecie "DD" lub "LL"
tree_path = "DD_KKpi/DecayTree"
track_type = tree_path.split("_")[0]

tree = file.Get(tree_path)

mass_var = "D_M"
cat_var = "D_BKGCAT"
low_limit = 1880
upp_limit = 2060


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

# wspolna srednia
mean = ROOT.RooRealVar("mean", "mean", 1900, low_limit, upp_limit)

# rdzen sygnalu
sigma1 = ROOT.RooRealVar("sigma1", "core width", 3, 0.1, 10)
gauss1 = ROOT.RooGaussian("gauss1", "Core Gauss", mass_roo, mean, sigma1)

# ogony
sigma2 = ROOT.RooRealVar("sigma2", "tail width", 20, 10, 40)
gauss2 = ROOT.RooGaussian("gauss2", "Tail Gauss", mass_roo, mean, sigma2)

# 70% sygnału to rdzen, a 30% to szerokie ogony
frac = ROOT.RooRealVar("frac", "core fraction", 0.7, 0.0, 1.0)

# zlozenie obu modeli
double_gauss = ROOT.RooAddPdf(
    "double_gauss",
    "Double Gaussian",
    ROOT.RooArgList(gauss1, gauss2),
    ROOT.RooArgList(frac),
)

# fitowanie
double_gauss.fitTo(reduced_data, ROOT.RooFit.PrintLevel(-1))


frame = mass_roo.frame(ROOT.RooFit.Title("gaussian fit for data cut (BKG_CAT < 30)"))

reduced_data.plotOn(frame)
double_gauss.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))

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
canvas.SaveAs(f"{mass_var}_double_gaussian_fit_pull_{track_type}_30.png")
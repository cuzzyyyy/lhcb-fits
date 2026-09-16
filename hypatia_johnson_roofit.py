import ROOT

# wyciszenie logow roofita, zeby nie miec syfu w terminalu
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

file = ROOT.TFile.Open("Bs_DsKstar_magup.root", "READ")

# definicja sciezki i wyciagniecie "DD" lub "LL" do zapisu nazwy
tree_path = "LL_KKpi/DecayTree"
track_type = tree_path.split("_")[0]

tree = file.Get(tree_path)

# definicja nazw zmiennych, okna masy dla mezonu D i liczby binow
mass_var = "D_M"
cat_var = "D_BKGCAT"
low_limit = 1880
upp_limit = 2060
nbins = 100

# inicjalizacja zmiennych roofitowych z odpowiednimi zakresami fizycznymi
mass_roo = ROOT.RooRealVar(mass_var, "mass [MeV/c2]", low_limit, upp_limit)
bkgcat_roo = ROOT.RooRealVar(cat_var, "background category", 0, 150)

# tworzenie datasetu z drzewa - ladujemy tylko te zmienne, po ktorych realnie tniemy albo fitujemy
dataset = ROOT.RooDataSet(
"data",
"data",
ROOT.RooArgSet(mass_roo, bkgcat_roo),
ROOT.RooFit.Import(tree),
)

# odcinanie czesci danych ze wzgledu na kategorie
reduced_data = dataset.reduce(ROOT.RooFit.Cut(f"{cat_var} < 30"))

# wspolna srednia dla obu modeli
mean = ROOT.RooRealVar("mean", "mean", 1970, 1960, 1990)

# double-sided hypatia - parametry
sigma_h = ROOT.RooRealVar("sigma_h", "sigma Hypatia", 5, 0.1, 20)
lambda_h = ROOT.RooRealVar("lambda_h", "lambda", -2.5, -10, 10)
zeta_h = ROOT.RooRealVar("zeta_h", "zeta", 0.0001, 1e-5, 0.1)
beta_h = ROOT.RooRealVar("beta_h", "beta", 0.0, -0.1, 0.1)
a1_h = ROOT.RooRealVar("a1_h", "a1", 2.0, 1.0, 10.0)
n1_h = ROOT.RooRealVar("n1_h", "n1", 2.0, 1.0, 10.0)
a2_h = ROOT.RooRealVar("a2_h", "a2", 2.0, 1.0, 10.0)
n2_h = ROOT.RooRealVar("n2_h", "n2", 2.0, 1.0, 10.0)

# zamrozenie niektorych parametrow ksztaltu, zeby fit nie zwariowal na surowych danych
lambda_h.setConstant(True)
zeta_h.setConstant(True)
beta_h.setConstant(True)

hypatia = ROOT.RooHypatia2(
"hypatia",
"Hypatia PDF",
mass_roo,
lambda_h,
zeta_h,
beta_h,
sigma_h,
mean,
a1_h,
n1_h,
a2_h,
n2_h,
)

# su johnson - parametry
sigma_j = ROOT.RooRealVar("sigma_j", "sigma Johnson", 10, 0.1, 30)
nu_j = ROOT.RooRealVar("nu_j", "nu", 0.0, -5.0, 5.0)
tau_j = ROOT.RooRealVar("tau_j", "tau", 1.0, 0.1, 10.0)

johnson = ROOT.RooJohnson(
"johnson", "Johnson SU PDF", mass_roo, mean, sigma_j, nu_j, tau_j
)

# skladanie obu rozkladow w jedna pdf
frac = ROOT.RooRealVar("frac", "Hypatia fraction", 0.5, 0.0, 1.0)
signal_model = ROOT.RooAddPdf(
"signal_model",
"Hypatia + Johnson",
ROOT.RooArgList(hypatia, johnson),
ROOT.RooArgList(frac),
)

# fitowanie
signal_model.fitTo(reduced_data, ROOT.RooFit.PrintLevel(-1))

# inicjalizacja ramki do rysowania z odpowiednim tytulem
frame = mass_roo.frame(
ROOT.RooFit.Title("hypatia + johnson fit for data cut (BKG_CAT < 30)")
)

# fitowanie wlasciwe z zapisaniem statusu i nalozenie wynikow na ramke
fit_result = signal_model.fitTo(reduced_data, ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Save())
reduced_data.plotOn(frame, ROOT.RooFit.Binning(nbins))
signal_model.plotOn(
    frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Precision(1e-5)
)

# liczenie chi2/ndf na podstawie wolnych parametrow po ficie
nparams = fit_result.floatParsFinal().getSize()
chi2_ndf = frame.chiSquare(nparams)


# pull jako obiekt tgraph
pull_graph = frame.pullHist()

# konwersja punktow pulla do klasycznego histogramu
pull_hist = ROOT.TH1F("pull_hist", ";mass [MeV/c^{2}];Pull", nbins, low_limit, upp_limit)
x_vals = pull_graph.GetX()
y_vals = pull_graph.GetY()
for i in range(pull_graph.GetN()):
    bin_idx = pull_hist.FindBin(x_vals[i])
    pull_hist.SetBinContent(bin_idx, y_vals[i])

# kosmetyka pulla
pull_hist.SetFillColor(ROOT.kAzure - 9)
pull_hist.SetLineColor(ROOT.kBlue + 2)
pull_hist.SetStats(0)
pull_hist.GetYaxis().SetRangeUser(-5, 5)
pull_hist.GetYaxis().SetTitleSize(0.08)
pull_hist.GetYaxis().SetTitleOffset(0.5)

canvas = ROOT.TCanvas("canvas", "Fit and Pull", 800, 800)

# podzial na dwa pady - gorny na fita (z wiekszym dolnym marginesem na os), dolny na pulla

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad1.SetBottomMargin(0.02)
pad1.Draw()

pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.3)
pad2.Draw()

pad1.SetBottomMargin(0.12) # robi odstęp na dole górnego panelu
pad2.SetTopMargin(0.05) # robi odstęp na górze dolnego panelu

pad1.cd()
frame.Draw()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.DrawLatex(0.15, 0.85, f"#chi^{{2}}/ndf = {chi2_ndf:.2f}")

pad2.cd()
pull_hist.Draw("HIST")

# linia zerowa na panelu pulla dla czytelnosci
line = ROOT.TLine(low_limit, 0, upp_limit, 0)
line.SetLineColor(ROOT.kBlack)
line.SetLineStyle(2)
line.Draw("SAME")

# automatyczny zapis wykorzystujący track_type (DD lub LL w nazwie)
canvas.SaveAs(f"1_{mass_var}_johnson_hypatia_fit_pull_{track_type}_30.png")

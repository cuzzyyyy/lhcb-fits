import ROOT
import numpy as np

# wyciszenie logow roofita, zeby nie miec syfu w terminalu
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

file = ROOT.TFile.Open("Bs_DsKstar_magup.root", "READ")

# definicja sciezki i wyciagniecie "DD" lub "LL" do zapisu nazwy
tree_path = "DD_KKpi/DecayTree"
track_type = tree_path.split("_")[0]

tree = file.Get(tree_path)

# definicja nazw zmiennych dla fita
mass_var = "B_M"
cat_var = "B_BKGCAT"
low_limit = 5000
upp_limit = 5700
nbins = 100

# zmienna po ktorej tniemy na biny (Kstar) - poprawione zakresy
bin_var_name = "B_M_K0Pi"
bin_low = 750
bin_upp = 1100

# inicjalizacja zmiennych roofitowych z odpowiednimi zakresami fizycznymi
mass_roo = ROOT.RooRealVar(mass_var, "mass [MeV/c2]", low_limit, upp_limit)
bkgcat_roo = ROOT.RooRealVar(cat_var, "background category", 0, 150)
bin_var_roo = ROOT.RooRealVar(bin_var_name, "binning mass [MeV/c2]", bin_low, bin_upp)

# tworzenie datasetu z drzewa - ladujemy TYLKO 3 zmienne (bez zdublowanej kategorii)
dataset = ROOT.RooDataSet(
    "data",
    "data",
    ROOT.RooArgSet(mass_roo, bkgcat_roo, bin_var_roo),
    ROOT.RooFit.Import(tree),
)

# wspolna srednia dla obu modeli
mean = ROOT.RooRealVar("mean", "mean", 5350, 5300, 5400)

# double-sided hypatia - parametry
sigma_h = ROOT.RooRealVar("sigma_h", "sigma Hypatia", 5, 0.1, 20)
lambda_h = ROOT.RooRealVar("lambda_h", "lambda", -2.5, -10, -1)

# ustawienie zakresow i zamrozenie zeby nie sypalo bledami o niestabilnosci w terminalu
zeta_h = ROOT.RooRealVar("zeta_h", "zeta", 0.0, 0.0, 1.0)
zeta_h.setConstant(True)
beta_h = ROOT.RooRealVar("beta_h", "beta", 0.0, -1.0, 1.0)
beta_h.setConstant(True)

a1_h = ROOT.RooRealVar("a1_h", "a1", 2.0, 1.0, 10.0)
n1_h = ROOT.RooRealVar("n1_h", "n1", 2.0, 1.0, 10.0)
a2_h = ROOT.RooRealVar("a2_h", "a2", 2.0, 1.0, 10.0)
n2_h = ROOT.RooRealVar("n2_h", "n2", 2.0, 1.0, 10.0)

# zamrozenie niektorych parametrow ksztaltu, zeby fit nie zwariowal na surowych danych
lambda_h.setConstant(True)

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

# najpierw fitujemy calosc zeby model zalapal poprawny ksztalt ogonow i asymetrii
# najpierw fitujemy calosc - podwojne ciecie na czysty sygnal obu czastek
cut_string_global = f"{cat_var} < 30 || {cat_var} == 50"
global_data = dataset.reduce(ROOT.RooFit.Cut(cut_string_global))

signal_model.fitTo(global_data, ROOT.RooFit.PrintLevel(-1))

# betonujemy parametry ogonow i ulamki, w binach pozwalamy plywac tylko sredniej i sigmom
a1_h.setConstant(True)
n1_h.setConstant(True)
a2_h.setConstant(True)
n2_h.setConstant(True)
nu_j.setConstant(True)
tau_j.setConstant(True)
frac.setConstant(True)

# wyliczenie korelacji pearsona zeby sprawdzic czy mozemy bezpiecznie mnozyc pdfy
hist2d = ROOT.TH2F(
    "hist2d", "korelacja", 100, low_limit, upp_limit, 100, bin_low, bin_upp
)

# bierzemy dane bezposrednio z drzewa dla szybkosci, z tym samym warunkiem na kategorie tla
tree.Draw(f"{bin_var_name}:{mass_var} >> hist2d", cut_string_global, "goff")

pearson_corr = hist2d.GetCorrelationFactor()

print("\n" + "=" * 50)
print(
    f"WSPOLCZYNNIK KORELACJI PEARSONA ({mass_var} x {bin_var_name}): {pearson_corr:.4f}"
)
print("=" * 50 + "\n")

# przygotowanie 9 rownych binow dla zmiennej tnacej
bin_edges = np.linspace(bin_low, bin_upp, 10)

# wielkie plotno na 9 wykresow w siatce 3x3
canvas = ROOT.TCanvas("canvas", "Binned Fits", 1800, 1200)
canvas.Divide(3, 3)

# zabezpieczenie obiektow graficznych przed wczesnym wyrzuceniem z pamieci przez garbage collector
keep_alive = []

# glowna petla po binach
for i in range(9):
    bin_min = bin_edges[i]
    bin_max = bin_edges[i + 1]

    # podwojny warunek na sygnał + ramy aktualnego binu
    cut_str = f"({cat_var} < 30 || {cat_var} == 50) && {bin_var_name} >= {bin_min} && {bin_var_name} < {bin_max}"
    # tymczasowy dataset wyciety tylko dla obecnego binu
    binned_data = dataset.reduce(ROOT.RooFit.Cut(cut_str))
    n_events = binned_data.numEntries()

    # przejscie do konkretnego kwadratu na plotnie
    canvas.cd(i + 1)

    # podzial kwadratu na dwa pady - gorny fita, dolny pulla
    pad_fit = ROOT.TPad(f"pad_fit_{i}", "fit", 0, 0.3, 1, 1.0)
    pad_pull = ROOT.TPad(f"pad_pull_{i}", "pull", 0, 0.0, 1, 0.3)

    pad_fit.SetBottomMargin(0.12)
    pad_pull.SetTopMargin(0.05)
    pad_pull.SetBottomMargin(0.3)

    pad_fit.Draw()
    pad_pull.Draw()
    keep_alive.extend([pad_fit, pad_pull])

    pad_fit.cd()

    # ostre ciecie na biny bez statystyki - nie rysujemy ramek, wrzucamy sam tekst na puste biale tlo
    if n_events < 30:
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.08)
        latex.DrawLatex(0.15, 0.5, f"not enough events ({n_events})")
        keep_alive.append(latex)
        continue

    # jesli statystyka jest, budujemy ramke i robimy fita
    frame = mass_roo.frame(
        ROOT.RooFit.Title(f"Bin {i+1}: {bin_min:.0f} - {bin_max:.0f} MeV")
    )
    binned_data.plotOn(frame, ROOT.RooFit.Binning(nbins))

    # pojedyncze wlasciwe fitowanie na zywych parametrach
    fit_result = signal_model.fitTo(
        binned_data, ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Save()
    )
    nparams = fit_result.floatParsFinal().getSize()

    signal_model.plotOn(
        frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Precision(1e-5)
    )
    frame.Draw()

    # wyciaganie wyfitowanych wartosci i ich bledow
    val_mean = mean.getVal()
    err_mean = mean.getError()
    val_sigma_h = sigma_h.getVal()
    err_sigma_h = sigma_h.getError()

    # wrzucenie wartosci chi2/ndf, sredniej i sigmy na ekran
    chi2_ndf = frame.chiSquare(nparams)
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.05)  # lekko mniejszy tekst zeby wszystko weszlo

    latex.DrawLatex(0.15, 0.85, f"#chi^{{2}}/ndf = {chi2_ndf:.2f}")
    latex.DrawLatex(0.15, 0.77, f"#mu = {val_mean:.2f} #pm {err_mean:.2f}")
    latex.DrawLatex(
        0.15, 0.69, f"#sigma_{{H}} = {val_sigma_h:.2f} #pm {err_sigma_h:.2f}"
    )

    keep_alive.extend([frame, latex])

    # rysowanie pulla
    pad_pull.cd()
    pull_graph = frame.pullHist()
    pull_hist = ROOT.TH1F(
        f"pull_hist_{i}", ";mass [MeV/c^{2}];Pull", nbins, low_limit, upp_limit
    )

    # konwersja punktow do klasycznego histogramu
    x_vals = pull_graph.GetX()
    y_vals = pull_graph.GetY()
    for j in range(pull_graph.GetN()):
        bin_idx = pull_hist.FindBin(x_vals[j])
        pull_hist.SetBinContent(bin_idx, y_vals[j])

    # kosmetyka pulla
    pull_hist.SetFillColor(ROOT.kAzure - 9)
    pull_hist.SetLineColor(ROOT.kBlue + 2)
    pull_hist.SetStats(0)
    pull_hist.GetYaxis().SetRangeUser(-5, 5)
    pull_hist.GetYaxis().SetTitleSize(0.1)
    pull_hist.GetYaxis().SetTitleOffset(0.3)
    pull_hist.GetYaxis().SetLabelSize(0.08)
    pull_hist.GetXaxis().SetLabelSize(0.08)
    pull_hist.Draw("HIST")

    # linia zerowa dla czytelnosci
    line = ROOT.TLine(low_limit, 0, upp_limit, 0)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(2)
    line.Draw("SAME")

    keep_alive.extend([pull_graph, pull_hist, line])

# automatyczny zapis calosci po skonczeniu petli
canvas.SaveAs(f"{mass_var}_binned_with_{bin_var_name}_{track_type}.png")

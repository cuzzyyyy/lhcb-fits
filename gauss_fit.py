import ROOT

file = ROOT.TFile.Open("Bs_DsKstar_magup.root", "READ")

if not file or file.IsZombie():
    print("Błąd: Nie można otworzyć pliku!")
    exit()

tree = file.Get("DD_KKpi/DecayTree")

canvas = ROOT.TCanvas("canvas", "Rozklad Masy", 800, 600)

h_mass = ROOT.TH1F("hist_M", "B_M", 100, 5000, 5700)

variable = "B_M"
tree.Draw(f"{variable} >> hist_M")

fit_result = h_mass.Fit("gaus", "SQ")

fit_func = h_mass.GetFunction("gaus")

if fit_func:
    mean = fit_func.GetParameter(1)
    mean_err = fit_func.GetParError(1)
    sigma = fit_func.GetParameter(2)
    sigma_err = fit_func.GetParError(2)

    print("\n--- WYNIKI DOPASOWANIA DLA B ---")
    print(f"Masa (Srednia): {mean:.2f} +/- {mean_err:.2f} MeV/c^2")
    print(f"Rozdzielczosc (Sigma): {sigma:.2f} +/- {sigma_err:.2f} MeV/c^2")

canvas.Draw()
canvas.SaveAs("B_fit_LL.png")

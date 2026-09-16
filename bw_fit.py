import ROOT

file = ROOT.TFile.Open("Bs_DsKstar_magup.root", "READ")

if not file or file.IsZombie():
    print("Błąd: Nie można otworzyć pliku!")
    exit()

tree = file.Get("LL_KKpi/DecayTree")

canvas = ROOT.TCanvas("canvas", "Rozklad Masy", 800, 600)

h_mass = ROOT.TH1F("hist_M", "B_M_K0Pi", 50, 850, 1050)

variable = "B_M_K0Pi"
tree.Draw(f"{variable} >> hist_M")

min_X = 850
max_X = 1050

# [0] to stała (amplituda/wysokość)
# [1] to masa (środek piku)
# [2] to Gamma (szerokość piku)
bw_func = ROOT.TF1("bw_func", "[0] * TMath::BreitWigner(x, [1], [2]) + [3]", min_X, max_X)

# musi byc jakis punkt startowy bo inaczej sie wywali
bw_func.SetParameter(0, h_mass.GetMaximum()) # amplituda - zaczyna od najwyzszego piku
bw_func.SetParameter(1, 892)    # masa - zaczyna od średniej z danych
bw_func.SetParameter(2, 50.0)   # szerokosc rozpadu - strzał
bw_func.SetParameter(3, 3.0)


bw_func.SetParLimits(1, 880, 910)  # środek musi być gdzieś w okolicach piku
bw_func.SetParLimits(2, 20, 100)
bw_func.SetParLimits(3, 0, 20)

# fit
fit_result = h_mass.Fit("bw_func", "SQL")

mean = bw_func.GetParameter(1)
mean_err = bw_func.GetParError(1)
gamma = bw_func.GetParameter(2)
gamma_err = bw_func.GetParError(2)

print("\n--- WYNIKI DOPASOWANIA BREITA-WIGNERA ---")
print(f"masa: {mean:.2f} +/- {mean_err:.2f} MeV/c^2")
print(f"szerokosc rozpadu: {gamma:.2f} +/- {gamma_err:.2f} MeV/c^2")

h_pull = h_mass.Clone("h_pull")
h_pull.Reset()
h_pull.SetTitle("pull")

for i in range(1, h_mass.GetNbinsX()+1):
    data = h_mass.GetBinContent(i)
    fit = bw_func.Eval(h_mass.GetBinCenter(i))
    uncertainty = h_mass.GetBinError(i)

    if uncertainty > 0:
        pull = (data-fit)/uncertainty
        h_pull.SetBinContent(i, pull)

canvas.Clear()
pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.3)
pad1.Draw()
pad2.Draw()

pad1.cd()
h_mass.Draw()

pad2.cd()
h_pull.SetFillColor(ROOT.kBlue)
h_pull.Draw("HIST")

canvas.Draw()
canvas.SaveAs("Kstar_fitandpull_restricted_LL_50.png")
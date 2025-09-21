#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

void plotThrustScatter() {
    // Open your ROOT file
    TFile *f = TFile::Open("output.root");  // change filename if needed
    if (!f || f->IsZombie()) {
        std::cerr << "Error: could not open file!" << std::endl;
        return;
    }

    TTree *tree = (TTree*)f->Get("tree");   // change "tree" if your TTree has a different name
    if (!tree) {
        std::cerr << "Error: tree not found!" << std::endl;
        return;
    }

    // Branch variables
    float thrust_angles[4];
    int mismatch;

    tree->SetBranchAddress("angles", thrust_angles); // assuming you stored as "angles[4]/F"
    tree->SetBranchAddress("mismatch", &mismatch);

    // Prepare containers
    std::vector<double> x_mismatch0, y_mismatch0;
    std::vector<double> x_mismatch1, y_mismatch1;

    // Loop over entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        double x = thrust_angles[3];
        double y = thrust_angles[0];
        if (mismatch == 0) {
            x_mismatch0.push_back(x);
            y_mismatch0.push_back(y);
        } else if (mismatch == 1) {
            x_mismatch1.push_back(x);
            y_mismatch1.push_back(y);
        }
    }

    // Make graphs
    TGraph *g0 = new TGraph(x_mismatch0.size(), x_mismatch0.data(), y_mismatch0.data());
    g0->SetMarkerStyle(20);
    g0->SetMarkerColor(kBlue);

    TGraph *g1 = new TGraph(x_mismatch1.size(), x_mismatch1.data(), y_mismatch1.data());
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kRed);

    // Draw
    TCanvas *c = new TCanvas("c", "Thrust scatter", 800, 600);
    g0->SetTitle("Thrust Angle Scatter;thrust_angles[3];thrust_angles[0]");
    g0->Draw("AP");
    g1->Draw("P SAME");

    // Legend
    auto legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(g0, "mismatch = 0", "p");
    legend->AddEntry(g1, "mismatch = 1", "p");
    legend->Draw();

    c->SaveAs("thrust_scatter.png"); // optional
}

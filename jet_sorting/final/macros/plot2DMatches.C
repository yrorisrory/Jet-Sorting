#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <iostream>

void plot2DMatches() {
    gStyle->SetOptStat(0); // Disable stat box

    const char* folders[] = {"HtHt/", "WbWb/", "ZtZt/"};
    const int nFolders = sizeof(folders)/sizeof(folders[0]);

    const char* tags[] = {"4-1", "5-1", "5-2", "6-1", "6-2",
                          "7-1", "7-2", "7-3", "8-1", "8-2", "8-3"};
    const int nTags = sizeof(tags)/sizeof(tags[0]);

    std::string output_path = "jet_sorting/final/graphs/";
    std::string input_path  = "jet_sorting/final/outputs/";

    // 2D histogram: X = decay product, Y = mass point
    TH2F* h2 = new TH2F("h2", "Normalized Mismatches (scaled to 4000);Decay Product;Mass Point",
                        nFolders, 0, nFolders,
                        nTags, 0, nTags);

    // Label x-axis
    for (int i = 0; i < nFolders; ++i) {
        std::string label = folders[i];
        if (label.back() == '/') label.pop_back();
        h2->GetXaxis()->SetBinLabel(i + 1, label.c_str());
    }
    // Label y-axis
    for (int j = 0; j < nTags; ++j) {
        h2->GetYaxis()->SetBinLabel(j + 1, tags[j]);
    }

    // Loop through datasets
    for (int i = 0; i < nFolders; ++i) {
        for (int j = 0; j < nTags; ++j) {
            std::string filePath = input_path + folders[i] + "clustered_" + tags[j] + ".root";
            TFile* f = TFile::Open(filePath.c_str());
            if (!f || f->IsZombie()) {
                std::cerr << "Error opening file: " << filePath << "\n";
                continue;
            }

            TDirectory* dir = (TDirectory*)f->Get("clusteringAnalyzerAll_nom");
            if (!dir) {
                std::cerr << "Directory 'clusteringAnalyzerAll_nom' not found in " << filePath << "\n";
                f->Close();
                continue;
            }

            TTree* tree = (TTree*)dir->Get("tree_nom");
            if (!tree) {
                std::cerr << "TTree 'tree_nom' not found in " << filePath << "\n";
                f->Close();
                continue;
            }

            int mismatch = -1;
            tree->SetBranchAddress("mismatch", &mismatch);

            int mismatches = 0;
            int matches = 0;
            Long64_t nEntries = tree->GetEntries();
            for (Long64_t entry = 0; entry < nEntries; ++entry) {
                tree->GetEntry(entry);
                if (mismatch == 0) matches++;
                else if (mismatch == 1) mismatches++;
            }

            // Apply formula mismatches * 4000 / (mismatches + matches)
            double value = 0.0;
            if (mismatches + matches > 0) {
                value = (static_cast<double>(mismatches) * 4000.0) / (mismatches + matches);
            }

            h2->SetBinContent(i + 1, j + 1, value);

            f->Close();
        }
    }

    // Draw
    TCanvas* c = new TCanvas("c", "2D Normalized Mismatches", 1000, 800);
    c->SetRightMargin(0.15);
    h2->Draw("COLZ TEXT");

    std::string outputName = output_path + "mismatch_0.40.png";
    c->SaveAs(outputName.c_str());
}

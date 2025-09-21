#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <cstdio>

void plotAllMismatch() {
    gStyle->SetOptStat(0);

    const char* folders[] = {"HtHt/", "WbWb/", "ZtZt/"};
    const int nFolders = sizeof(folders)/sizeof(folders[0]);

    const char* tags[] = {"4-1", "5-1", "5-2", "6-1", "6-2",
                          "7-1", "7-2", "7-3", "8-1", "8-2", "8-3"};
    const int nTags = sizeof(tags)/sizeof(tags[0]);

    const char* input_path  = "jet_sorting/final/outputs/";
    const char* output_path = "jet_sorting/final/graphs/";

    char folder_label[64];
    char filePath[512];
    char outname[512];
    char title[256];

    for (int i = 0; i < nFolders; ++i) {
        // copy folder label without trailing '/'
        std::snprintf(folder_label, sizeof(folder_label), "%s", folders[i]);
        size_t flen = std::strlen(folder_label);
        if (flen > 0 && folder_label[flen-1] == '/') folder_label[flen-1] = '\0';

        for (int j = 0; j < nTags; ++j) {
            std::snprintf(filePath, sizeof(filePath), "%s%sclustered_%s.root",
                          input_path, folders[i], tags[j]);

            TFile* f = TFile::Open(filePath);
            if (!f || f->IsZombie()) {
                std::fprintf(stderr, "Error opening file: %s\n", filePath);
                if (f) f->Close();
                continue;
            }

            TDirectory* dir = (TDirectory*)f->Get("clusteringAnalyzerAll_nom");
            if (!dir) {
                std::fprintf(stderr, "Directory 'clusteringAnalyzerAll_nom' not found in %s\n", filePath);
                f->Close();
                continue;
            }

            TTree* tree = (TTree*)dir->Get("tree_nom");
            if (!tree) {
                std::fprintf(stderr, "TTree 'tree_nom' not found in %s\n", filePath);
                f->Close();
                continue;
            }

            // assume angles stored as a fixed array "angles[4]/F"
            float thrust_angles[4];
            int mismatch = -1;
            tree->SetBranchAddress("thrust_angles", thrust_angles);
            tree->SetBranchAddress("mismatch", &mismatch);

            // Use TGraph and SetPoint to avoid STL vectors
            TGraph* g0 = new TGraph(); // mismatch == 0 -> blue
            TGraph* g1 = new TGraph(); // mismatch == 1 -> red

            Long64_t nEntries = tree->GetEntries();
            for (Long64_t e = 0; e < nEntries; ++e) {
                tree->GetEntry(e);
                double x = thrust_angles[3];
                double y = thrust_angles[0];
                if (mismatch == 0) {
                    g0->SetPoint(g0->GetN(), x, y);
                } else if (mismatch == 1) {
                    g1->SetPoint(g1->GetN(), x, y);
                }
            }

            g0->SetMarkerStyle(20);
            g0->SetMarkerSize(0.5);
            g0->SetMarkerColor(kBlue);

            g1->SetMarkerStyle(20);
            g1->SetMarkerSize(0.5);
            g1->SetMarkerColor(kRed);

            std::snprintf(title, sizeof(title),
                          "Distance from Thrust Axis (%s, %s);Second Chi;First Chi",
                          folder_label, tags[j]);
            g0->SetTitle(title);

            TCanvas* c = new TCanvas("c", "Thrust scatter", 900, 700);
            g0->Draw("AP");
            g1->Draw("P SAME");

            TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.9);
            legend->AddEntry(g0, "All large hard jets correct", "p");
            legend->AddEntry(g1, "At least one hard jet matched to wrong superjet", "p");
            legend->Draw();

            std::snprintf(outname, sizeof(outname), "%s%sscatter_%s.png", output_path, folders[i], tags[j]);
            c->SaveAs(outname);

            // cleanup
            delete legend;
            delete c;
            delete g0;
            delete g1;
            f->Close();
        }
    }
}

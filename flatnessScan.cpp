//
// Build:  clang++ -std=c++17 flatnessScan.cpp Points.cpp `root-config --cflags --libs` -o flatnessScan
//
//------------------------------------------------------------------------------
// File: flatnessScan.cpp
//
// Description:
//   This program analyzes a set of 3D points representing measured surface
//   coordinates and evaluates the "flatness" of the surface by fitting an
//   optimal 3D plane to the ensemble of points and histogramming their
//   deviations from that best-fit plane.
//
//   The program performs the following steps:
//     1. Reads a text file containing one point per line. Each line is expected
//        to contain four values: an integer point index followed by the
//        X, Y, and Z coordinates in millimeters.
//     2. Stores the coordinates in a vector of Point objects (see Points.h).
//     3. Uses the ROOT Minuit2 minimizer to determine the parameters (ax, ay, az)
//        of the best-fit plane defined by the equation:
//
//             ax * X + ay * Y + az * (Z + offset) = 1
//
//        minimizing the total χ² = Σ [ (ax*X + ay*Y + az*(Z+offset) - 1)² / (ax²+ay²+az²) ]
//
//     4. Computes the resulting χ², standard deviation, and plane normal
//        normalization (|a| and 1/|a|), and prints them to the console.
//     5. Fills ROOT histograms for each coordinate (X, Y, Z) and for the
//        deviation of each point from the fitted plane.
//     6. Produces a 2D scatter plot of Y vs. X and displays all histograms
//        and the scatter plot in interactive ROOT canvases.
//     7. Detects whether the data lie on a regular (Nx × Ny) grid. If so,
//        constructs a 2D “flatness map” histogram colored by Z values.
//     8. Writes all histograms and the TGraph to an output ROOT file
//        ("output.root") and displays all results.
//
// Input:
//   A text file (e.g. "points.csv") with four columns per line:
//       index,  X,  Y,  Z
//
// Output:
//   - Console summary of fit results (χ², plane coefficients, flatness).
//   - ROOT file "output.root" containing histograms and scatter plots.
//   - ROOT canvases displaying coordinate distributions and residuals.
//
// Usage example:
//   $ ./flatnessScan my_points.csv
//
// Dependencies:
//   - ROOT framework (TFile, TH1D, TH2D, TGraph, TCanvas, TApplication, Minimizer)
//   - Points.h / Points.cpp for reading input data
//   - GridFinder.h for grid detection
//
// Author: Luciano Ristori
// Version: 1.0
// Date: October 2025
//------------------------------------------------------------------------------

// ROOT and standard headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <limits>
#include <iomanip>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TColor.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "Points.h"
#include "GridFinder.h"

//------------------------------------------------------------------------------
// Program version (update when functionality changes)
const std::string FLATNESS_SCAN_VERSION = "1.1.0 (October 2025)";


using std::cout;
using std::endl;
using std::vector;

std::vector<double> X, Y, Z;
double offset = 400.0;

//------------------------------------------------------------------------------
// Helper classes for formatted console output
//------------------------------------------------------------------------------

class ScientificPrecision {
public:
    ScientificPrecision(std::ostream& os, int precision)
        : os_(os), old_precision_(os.precision()), old_flags_(os.flags()) {
        os_ << std::scientific << std::setprecision(precision);
    }
    ~ScientificPrecision() {
        os_.precision(old_precision_);
        os_.flags(old_flags_);
    }
private:
    std::ostream& os_;
    std::streamsize old_precision_;
    std::ios_base::fmtflags old_flags_;
};

class FloatingPointPrecision {
public:
    FloatingPointPrecision(std::ostream& os, int precision)
        : os_(os), old_precision_(os.precision()) {
        os_ << std::setprecision(precision);
    }
    ~FloatingPointPrecision() {
        os_.precision(old_precision_);
    }
private:
    std::ostream& os_;
    std::streamsize old_precision_;
};

//------------------------------------------------------------------------------
// χ² Function for plane fitting
//------------------------------------------------------------------------------

double chi2Func(const double *x) {
    double ax = x[0], ay = x[1], az = x[2];
    double chi2 = 0.0;

    for (size_t i = 0; i < X.size(); ++i) {
        double delta = ax*X[i] + ay*Y[i] + az*(Z[i] + offset) - 1.0;
        chi2 += delta * delta / (ax*ax + ay*ay + az*az);
    }
    return chi2;
}

//------------------------------------------------------------------------------
// Main program
//------------------------------------------------------------------------------

int main(int argc, char *argv[]) {

// Usage:
//
//   ./flatnessScan input.csv [output.root]
//
//------------------------------------------------------------------------------
// 1. Parse command-line arguments and initialize ROOT application.
//------------------------------------------------------------------------------
////
// If no output file is specified, defaults to "output.root".
// If the name given does not end with ".root", the extension is added automatically.
//

	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " input.csv [output.root]" << std::endl;
		return 1;
	}
	
	std::cout << "\nFlatnessScan version " << FLATNESS_SCAN_VERSION << std::endl;

	std::string filename = argv[1];
	std::string outname = (argc >= 3) ? argv[2] : "output.root";

	// Append ".root" if missing (case-insensitive)
	if (outname.size() < 5 || 
		(outname.substr(outname.size() - 5) != ".root" &&
		 outname.substr(outname.size() - 5) != ".ROOT")) {
		outname += ".root";
	}

	// Initialize ROOT GUI
	TApplication app("app", &argc, argv);
	gROOT->SetBatch(false); // enable GUI
	TH1::AddDirectory(kTRUE);

    // 2. Read 3D points from input file
    
    int n = 4;
    std::vector<Point> points = readPoints(filename, n);
    if (points.empty()) {
        std::cerr << "No valid points found. Exiting." << std::endl;
        return 1;
    }

    for (auto &p : points) {
        X.push_back(p.coords[1]);
        Y.push_back(p.coords[2]);
        Z.push_back(p.coords[3]);
    }

    cout << "Read " << points.size() << " valid points." << endl;

    // 3. Fit a 3D plane using Minuit2
    
    cout << "\nFitting 3D plane..." << endl;
    ROOT::Math::Minimizer* min =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

    min->SetMaxFunctionCalls(1000000);
    min->SetMaxIterations(1000);
    min->SetTolerance(0.001);
    min->SetPrintLevel(0);

    ROOT::Math::Functor f(&chi2Func, 3);
    double step[3] = {0.001, 0.001, 0.001};
    double variable[3] = {0.0, 0.0, 1.0 / offset};
    min->SetFunction(f);

    min->SetVariable(0, "ax", variable[0], step[0]);
    min->SetVariable(1, "ay", variable[1], step[1]);
    min->SetVariable(2, "az", variable[2], step[2]);
    min->Minimize();

    const double *res = min->X();
    const double *err = min->Errors();

    double ax = res[0], ay = res[1], az = res[2];
    double ax_e = err[0], ay_e = err[1], az_e = err[2];
    double minChi2 = min->MinValue();

    {
        ScientificPrecision sp(cout, 2);
        cout << "\n----------------------------------\n";
        cout << "  Plane fit summary\n";
        cout << "  ax = " << ax << " ± " << ax_e << "\n";
        cout << "  ay = " << ay << " ± " << ay_e << "\n";
        cout << "  az = " << az << " ± " << az_e << "\n";
    }

    {
        FloatingPointPrecision fpp(cout, 4);
        cout << "  σ = " << 1000. * sqrt(minChi2 / X.size()) << " µm\n";
        cout << "----------------------------------\n";
    }

    double moda = sqrt(ax*ax + ay*ay + az*az);
    double invModa = 1.0 / moda;
    cout << "\n|a| = " << moda << "   1/|a| = " << invModa << " [mm]" << endl;
    cout << "Offset: " << offset << " [mm]" << endl;

    // 4. Determine coordinate ranges
    
    std::vector<double> mins(n, std::numeric_limits<double>::max());
    std::vector<double> maxs(n, std::numeric_limits<double>::lowest());
    for (const auto &p : points)
        for (int i = 0; i < n; ++i) {
            if (p.coords[i] < mins[i]) mins[i] = p.coords[i];
            if (p.coords[i] > maxs[i]) maxs[i] = p.coords[i];
        }

	// 5. Create histograms for X, Y, Z, and residuals
	//    → Provides coordinate distributions and flatness residuals for visualization

    TFile outfile(outname.c_str(), "RECREATE");

    std::vector<TH1D*> hists;

    for (int i = 0; i < n; ++i) {
        double min = mins[i], max = maxs[i];
        if (min == max) { min -= 0.5; max += 0.5; }
        double margin = 0.5 * (max - min);
        int nBins = static_cast<int>((max - min + 2 * margin) * 1000 + 0.5);

        std::string hname, htitle, xaxis;
        
        if (i == 0) { hname = "hn"; htitle = "Point sequence number"; xaxis = "n"; }
        else if (i == 1) { hname = "hX"; htitle = "X Coordinate Distribution"; xaxis = "X [mm]"; }
        else if (i == 2) { hname = "hY"; htitle = "Y Coordinate Distribution"; xaxis = "Y [mm]"; }
        else if (i == 3) { hname = "hZ"; htitle = "Z Coordinate Distribution"; xaxis = "Z [mm]"; }
        else { hname = "hCoord" + std::to_string(i + 1); htitle = "Coordinate " + std::to_string(i + 1); xaxis = "Value"; }

        auto *h = new TH1D(hname.c_str(), htitle.c_str(),
                           nBins, min - margin, max + margin);
        h->GetXaxis()->SetTitle(xaxis.c_str());
        h->GetYaxis()->SetTitle("Counts");
        hists.push_back(h);

		// For Z coordinate (i == 3), also create a second histogram
    	// to store residuals (deviations from the fitted 3D plane).
        if (i == 3) {
            auto *hDev = new TH1D("hDeviations", "Deviations from 3D Plane Fit",
                                  nBins, min - margin, max + margin);
            hDev->GetXaxis()->SetTitle("Residual [mm]");
            hDev->GetYaxis()->SetTitle("Counts");
            hists.push_back(hDev);
        }
    }

    for (const auto &p : points) {
        for (int i = 0; i < n; ++i)
            hists[i]->Fill(p.coords[i]);
        double delta = (ax*p.coords[1] + ay*p.coords[2] + az*(p.coords[3] + offset) - 1.0) * invModa;
        hists[4]->Fill(delta);
    }
    
    // write code version to histogram file
    
    TNamed versionTag("FlatnessScanVersion", FLATNESS_SCAN_VERSION.c_str());
	versionTag.Write();

    for (auto h : hists) h->Write();

    // 6. 2D Scatter plot of Y vs X
    
    TGraph* g2 = new TGraph(points.size());
    for (size_t i = 0; i < points.size(); ++i)
    g2->SetPoint(i, points[i].coords[1], points[i].coords[2]);
    g2->SetName("g2_xy");
    g2->SetTitle("Y vs X");
    g2->Write();

    // 7. Flatness color map if grid is regular
    
    std::vector<std::pair<double,double>> xy;
    xy.reserve(points.size());
    for (const auto &p : points)
        xy.emplace_back(p.coords[1], p.coords[2]);
        
    // Analyze (X, Y) points to determine if they form a regular Nx×Ny grid.
	// If yes, create a color-coded 2D histogram of Z values — the "flatness map".

    auto grid = GridFinder::analyze(xy);
    TH2D *hZ = nullptr;

    if (grid.regularX && grid.regularY) {
        hZ = new TH2D("hZMap", "Flatness Map;X [mm];Y [mm];Z [mm]",
                      grid.Nx, grid.xMin - grid.dx/2, grid.xMax + grid.dx/2,
                      grid.Ny, grid.yMin - grid.dy/2, grid.yMax + grid.dy/2);

        std::map<std::pair<int,int>, std::vector<double>> bins;

        for (const auto& p : points) {
            int ix = static_cast<int>(std::round((p.coords[1] - grid.xMin) / grid.dx));
            int iy = static_cast<int>(std::round((p.coords[2] - grid.yMin) / grid.dy));
            bins[{ix, iy}].push_back(p.coords[3]);
        }

        for (const auto& [idx, zs] : bins) {
            double zmean = std::accumulate(zs.begin(), zs.end(), 0.0) / zs.size();
            hZ->SetBinContent(idx.first + 1, idx.second + 1, zmean);
        }
		hZ->SetStats(0);  // disables stats box for this histogram
        hZ->Write();
    } else {
        std::cerr << "Warning: points are not on a regular grid — skipping flatness map.\n";
    }

    // 8. Display results
    int canvasWidth = 800, canvasHeight = 600;

    for (size_t i = 0; i < hists.size(); ++i) {
        std::string cname = "cHist_" + std::to_string(i + 1);
        TCanvas *c = new TCanvas(cname.c_str(), hists[i]->GetTitle(), 50 + i * 30, 50 + i * 30,
                                 canvasWidth, canvasHeight);
        c->Connect("Closed()", "TApplication", gApplication, "Terminate()");
        hists[i]->Draw();
        c->Update();
    }

    TCanvas *c2 = new TCanvas("c2", "2D Scatter (Y vs X)", 900, 150, 700, 600);
    c2->Connect("Closed()", "TApplication", gApplication, "Terminate()");
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.8);
    g2->SetMarkerColor(kBlack);
    g2->Draw("AP");
    c2->Update();

    if (hZ) {
        TCanvas *cMap = new TCanvas("cMap", "Flatness Map", 1650, 150, 800, 650);
        gStyle->SetPalette(kBird);
        cMap->SetLeftMargin(0.15);
        cMap->SetRightMargin(0.18);
        cMap->SetBottomMargin(0.12);
        cMap->SetTopMargin(0.08);
        hZ->SetStats(0);
        hZ->GetXaxis()->SetTitleOffset(1.2);
        hZ->GetYaxis()->SetTitleOffset(1.6);
        hZ->Draw("COLZ");
        cMap->Update();
	}
	//------------------------------------------------------------------------------
	// 9. Run ROOT GUI loop
	//------------------------------------------------------------------------------
	//
	// Close the output file before entering the interactive ROOT GUI loop.
	// Canvases remain accessible even after the file is closed.
	//
	
	
	std::cout << "\nHistograms written to " << outname << std::endl;
	std::cout << "\nHit ctrl-c to exit" << std:: endl;
	
		
	
	// detach histograms from file so they survive after outfile.Close()
	for (TH1D* hist : hists) {
		hist->SetDirectory(nullptr);
		hist->Write();
	}	
	// detach and write scatter plot
	g2->Write();

	// write flatness map if present
	if (hZ) {
		hZ->SetDirectory(nullptr);
		hZ->Write();
	}
		outfile.Close();

	// Enter the ROOT GUI event loop — close all canvases or press Ctrl+C to exit.
	
	app.Run();
	
    return 0;
}

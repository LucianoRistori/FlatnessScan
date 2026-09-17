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
//     3b. Fits a second plane -- the ASME Y14.5 "minimum zone" plane -- by
//        minimizing the peak-to-valley span of the same residuals instead of
//        their sum of squares. This is the standard engineering definition
//        of flatness: the distance between the two parallel planes of
//        minimum separation that contain all the measured points.
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
//   - Console summary of fit results (χ², plane coefficients, flatness),
//     mirrored into a matching ".log" text file alongside the ROOT output.
//   - ROOT file (histograms and scatter plots) plus the matching ".log" file.
//   - PDF snapshot of the deviations-from-fit histogram ("..._deviations.pdf").
//   - PDF snapshot of the 2D flatness color map, with its color-scale legend
//     ("..._flatnessMap.pdf"; only produced when the (X,Y) points form a
//     regular grid).
//   - One-page presentation-ready summary image ("..._summary.png", 16:9)
//     combining the best-fit plane residual RMS, the GD&T ASME Y14.5
//     flatness result, the deviations histogram, and the 2D color map with
//     its legend -- meant to be dropped straight into a PowerPoint or
//     Keynote slide.
//   - ROOT canvases displaying coordinate distributions and residuals.
//   - Unless an explicit output name is given, each run gets its own fresh
//     subfolder named after the input file, with a sequential number on the
//     FOLDER (not the files): "myPoints.csv" -> "myPoints_1/myPoints.root"
//     + "myPoints_1/myPoints.log" + "myPoints_1/myPoints_deviations.pdf"
//     + "myPoints_1/myPoints_flatnessMap.pdf" + "myPoints_1/myPoints_summary.png"
//     on the first run, "myPoints_2/..." on the next, and so on.
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

#define FLATNESSSCAN_VERSION "v1.0"


// ROOT and standard headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>
#include <cstdlib>
#include <limits>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <set>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TBox.h"
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
const std::string FLATNESS_SCAN_VERSION = "1.5.7 (August 2026)";


using std::cout;
using std::endl;
using std::vector;

std::vector<double> X, Y, Z;
double offset = 400.0;

#include <cstdlib>
#include <stdexcept>



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
// TeeBuf: a streambuf that duplicates every character written to it into two
// underlying streambufs. Used to mirror std::cout / std::cerr into a log
// file in addition to the console, without touching every individual
// cout/cerr statement elsewhere in the program.
//------------------------------------------------------------------------------

class TeeBuf : public std::streambuf {
public:
    TeeBuf(std::streambuf* sb1, std::streambuf* sb2) : sb1_(sb1), sb2_(sb2) {}
protected:
    int overflow(int c) override {
        if (c == EOF) return !EOF;
        int r1 = sb1_->sputc(static_cast<char>(c));
        int r2 = sb2_->sputc(static_cast<char>(c));
        return (r1 == EOF || r2 == EOF) ? EOF : c;
    }
    int sync() override {
        int r1 = sb1_->pubsync();
        int r2 = sb2_->pubsync();
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }
private:
    std::streambuf* sb1_;
    std::streambuf* sb2_;
};

//------------------------------------------------------------------------------
// χ² Function for plane fitting (least-squares plane)
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
// Minimum-zone (Chebyshev) objective for the standard ASME Y14.5 flatness
// definition. For a candidate plane, this is the peak-to-valley span of the
// signed perpendicular distances of all points from that plane. Minimizing
// this over the plane parameters finds the orientation for which the two
// parallel planes containing all the points are closest together; that
// minimum separation is the reported flatness value.
//------------------------------------------------------------------------------

double minimaxFunc(const double *x) {
    double ax = x[0], ay = x[1], az = x[2];
    double moda = std::sqrt(ax*ax + ay*ay + az*az);
    if (moda < 1e-12) return std::numeric_limits<double>::max();

    double dmax = std::numeric_limits<double>::lowest();
    double dmin = std::numeric_limits<double>::max();
    for (size_t i = 0; i < X.size(); ++i) {
        double delta = (ax*X[i] + ay*Y[i] + az*(Z[i] + offset) - 1.0) / moda;
        if (delta > dmax) dmax = delta;
        if (delta < dmin) dmin = delta;
    }
    return dmax - dmin;
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
// 1. Parse command-line arguments, set up the log file, and initialize the
//    ROOT application.
//------------------------------------------------------------------------------
////
// If no output file is specified, defaults to "output.root".
// If the name given does not end with ".root", the extension is added automatically.
//

	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " input.csv [output.root]" << std::endl;
		return 1;
	}

	std::string filename = argv[1];

	// Default output basename tracks the input file: strip any directory
	// path and extension from the input filename. An explicit output name
	// on the command line (argv[2]) still overrides everything below.
	//
	// Every run is collected into its own fresh subfolder named after the
	// input file, created alongside it -- the SEQUENTIAL NUMBER is on the
	// folder, not on the files inside it (e.g. input "myPoints.csv"
	// produces "myPoints_1/myPoints.root" + "myPoints_1/myPoints.log" on
	// the first run, "myPoints_2/myPoints.root" + "myPoints_2/myPoints.log"
	// on the next, and so on).
	std::string inputBase = filename;
	std::string inputDir;  // directory the input file lives in, including trailing slash (may be empty)
	size_t slashPos = inputBase.find_last_of("/\\");
	if (slashPos != std::string::npos) {
		inputDir = inputBase.substr(0, slashPos + 1);
		inputBase = inputBase.substr(slashPos + 1);
	}
	size_t dotPos = inputBase.find_last_of('.');
	if (dotPos != std::string::npos) inputBase = inputBase.substr(0, dotPos);

	// outDir holds the directory all outputs for this run land in (the
	// fresh numbered subfolder in the default case, or whatever directory
	// the explicit output name points at). Used later to place the PDF
	// snapshots of the deviations histogram and the flatness map.
	std::string outname;
	std::string outDir;
	if (argc >= 3) {
		outname = argv[2];
		size_t outSlash = outname.find_last_of("/\\");
		outDir = (outSlash != std::string::npos) ? outname.substr(0, outSlash + 1) : "";
	} else {
		int runN = 1;
		while (true) {
			outDir = inputDir + inputBase + "_" + std::to_string(runN) + "/";
			if (!std::filesystem::exists(outDir)) break;
			++runN;
		}
		std::filesystem::create_directories(outDir);
		outname = outDir + inputBase + ".root";
	}
	// Append ".root" if missing (case-insensitive)
	if (outname.size() < 5 ||
		(outname.substr(outname.size() - 5) != ".root" &&
		 outname.substr(outname.size() - 5) != ".ROOT")) {
		outname += ".root";
	}

	// Derive a matching log file name (same base name, ".log" extension) and
	// mirror everything written to std::cout / std::cerr into it, so every
	// run leaves a complete text record alongside the ROOT output file.
	std::string logname = outname.substr(0, outname.size() - 5) + ".log";
	std::ofstream logFile(logname.c_str());
	if (!logFile) {
		std::cerr << "Warning: could not open log file \"" << logname
		          << "\" -- continuing without a log file." << std::endl;
	}
	TeeBuf teeOutBuf(std::cout.rdbuf(), logFile.rdbuf());
	TeeBuf teeErrBuf(std::cerr.rdbuf(), logFile.rdbuf());
	std::streambuf* origCoutBuf = std::cout.rdbuf(&teeOutBuf);
	std::streambuf* origCerrBuf = std::cerr.rdbuf(&teeErrBuf);

	cout << "\n====================================\n";
	cout << " FlatnessScan " << FLATNESSSCAN_VERSION << " — Luciano Ristori\n";
	cout << " Built: " << __DATE__ << " " << __TIME__ << endl;
	cout << "====================================\n";
	cout << "Log file: " << logname << endl;

	// Initialize ROOT GUI
	TApplication app("app", &argc, argv);
	gROOT->SetBatch(false); // enable GUI
	TH1::AddDirectory(kTRUE);

    // 2. Read 3D points from input file
    
	cout << "Data file: " << filename << endl;
    int n = 3;
    std::vector<Point> points = readPoints(filename, n);
    if (points.empty()) {
        std::cerr << "No valid points found. Exiting." << std::endl;
        return 1;
    }

    for (auto &p : points) {
        X.push_back(p.coords[0]);
        Y.push_back(p.coords[1]);
        Z.push_back(p.coords[2]);
    }

    cout << "Read " << points.size() << " valid points." << endl;

    // 3. Fit a 3D plane using Minuit2 (least-squares plane)

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

    // 3b. Fit the ASME Y14.5 "minimum zone" plane and report the standard
    //     flatness value. This is a second, independent minimization of the
    //     same three plane parameters, seeded from the least-squares result,
    //     but minimizing the peak-to-valley span of the residuals instead of
    //     their sum of squares. No plots are produced for this -- the result
    //     is a single number, printed to the console (and, via the tee above,
    //     to the log file).

    cout << "\nFitting minimum-zone plane (ASME Y14.5 flatness definition)..." << endl;
    ROOT::Math::Minimizer* minMZ =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");

    minMZ->SetMaxFunctionCalls(1000000);
    minMZ->SetMaxIterations(1000);
    minMZ->SetTolerance(0.001);
    minMZ->SetPrintLevel(0);

    ROOT::Math::Functor fMZ(&minimaxFunc, 3);
    double stepMZ[3]     = {0.0001, 0.0001, 0.0001};
    double variableMZ[3] = {ax, ay, az};   // seed from the least-squares plane
    minMZ->SetFunction(fMZ);

    minMZ->SetVariable(0, "ax", variableMZ[0], stepMZ[0]);
    minMZ->SetVariable(1, "ay", variableMZ[1], stepMZ[1]);
    minMZ->SetVariable(2, "az", variableMZ[2], stepMZ[2]);
    minMZ->Minimize();

    double flatness = minMZ->MinValue();  // peak-to-valley span at the optimum [mm]

    {
        FloatingPointPrecision fpp(cout, 4);
        cout << "\n----------------------------------\n";
        cout << "  Flatness (ASME Y14.5 minimum-zone method)\n";
        cout << "  Flatness = " << flatness << " mm  (" << flatness * 1000.0 << " µm)\n";
        cout << "----------------------------------\n";
    }

    // 4. Determine coordinate ranges

    std::vector<double> mins(n, std::numeric_limits<double>::max());
    std::vector<double> maxs(n, std::numeric_limits<double>::lowest());
    for (const auto &p : points)
        for (int i = 0; i < n; ++i) {
            if (p.coords[i] < mins[i]) mins[i] = p.coords[i];
            if (p.coords[i] > maxs[i]) maxs[i] = p.coords[i];
        }

        cout<<"ranges done" << endl;

    // Precompute the fit residuals (perpendicular distance of each point
    // from the best-fit plane) up front, and their own min/max. The
    // deviations histogram below is sized from THIS range, not from the
    // Z coordinate's own range: Z can sit far from zero (e.g. a CMM
    // reading around 6-7 mm), while the residuals are tiny deviations
    // centered near zero, so reusing the Z axis range left every residual
    // landing outside it -- a histogram that was 100% underflow.
    std::vector<double> devVals;
    devVals.reserve(points.size());
    double devMin = std::numeric_limits<double>::max();
    double devMax = std::numeric_limits<double>::lowest();
    for (const auto &p : points) {
        double d = (ax*p.coords[0] + ay*p.coords[1] + az*(p.coords[2] + offset) - 1.0) * invModa;
        devVals.push_back(d);
        if (d < devMin) devMin = d;
        if (d > devMax) devMax = d;
    }

	// 5. Create histograms for X, Y, Z, and residuals
	//    → Provides coordinate distributions and flatness residuals for visualization

    TFile outfile(outname.c_str(), "RECREATE");

    std::vector<TH1D*> hists;

    cout << "n = " << n << endl;

    // Bin width for the Z-coordinate and residual/deviation histograms:
    // the default formula below (used for X and Y) works out to a 1 micron
    // bin, which is far finer than useful for a sub-millimeter surface
    // measurement -- 5 microns per bin instead for those two.
    const double coarseBinWidth = 0.005;  // mm (5 microns)

    for (int i = 0; i < n; ++i) {
        double min = mins[i], max = maxs[i];
        if (min == max) { min -= 0.5; max += 0.5; }
        double margin = 0.5 * (max - min);
        int nBins = static_cast<int>((max - min + 2 * margin) * 1000 + 0.5);
        if (i == 2) {
            nBins = std::max(1, static_cast<int>((max - min + 2 * margin) / coarseBinWidth + 0.5));
        }

        std::string hname, htitle, xaxis;

        if (i == 0) { hname = "hX"; htitle = "X Coordinate Distribution"; xaxis = "X [mm]"; }
        else if (i == 1) { hname = "hY"; htitle = "Y Coordinate Distribution"; xaxis = "Y [mm]"; }
        else if (i == 2) { hname = "hZ"; htitle = "Z Coordinate Distribution"; xaxis = "Z [mm]"; }
        else { hname = "hCoord" + std::to_string(i + 1); htitle = "Coordinate " + std::to_string(i + 1); xaxis = "Value"; }

        auto *h = new TH1D(hname.c_str(), htitle.c_str(),
                           nBins, min - margin, max + margin);
        h->GetXaxis()->SetTitle(xaxis.c_str());
        h->GetYaxis()->SetTitle("Counts");
        hists.push_back(h);

		// For Z coordinate (i == 2), also create a second histogram to
    	// store residuals (deviations from the fitted 3D plane) -- using
    	// their OWN range (devMin/devMax, computed above), not the Z
    	// coordinate's range used for "h" above: the two are unrelated
    	// quantities and can sit nowhere near each other numerically.
        if (i == 2) {
            double devMinL = devMin, devMaxL = devMax;
            if (devMinL == devMaxL) { devMinL -= 0.5; devMaxL += 0.5; }
            double devMargin = 0.5 * (devMaxL - devMinL);
            int nBinsDev = std::max(1, static_cast<int>((devMaxL - devMinL + 2 * devMargin) / coarseBinWidth + 0.5));
            auto *hDev = new TH1D("hDeviations", "Deviations from 3D Plane Fit",
                                  nBinsDev, devMinL - devMargin, devMaxL + devMargin);
            hDev->GetXaxis()->SetTitle("Residual [mm]");
            hDev->GetYaxis()->SetTitle("Counts");
            hists.push_back(hDev);
        }
    }

    cout << "histograms done" << endl;

    for (size_t idx = 0; idx < points.size(); ++idx) {
        const auto &p = points[idx];
        for (int i = 0; i < n; ++i)
            hists[i]->Fill(p.coords[i]);
        	hists[3]->Fill(devVals[idx]);
    }

    // write code version to histogram file

    TNamed versionTag("FlatnessScanVersion", FLATNESS_SCAN_VERSION.c_str());
	versionTag.Write();

    for (auto h : hists) h->Write();

    // 6. 2D Scatter plot of Y vs X

    TGraph* g2 = new TGraph(points.size());
    for (size_t i = 0; i < points.size(); ++i)
    g2->SetPoint(i, points[i].coords[0], points[i].coords[1]);
    g2->SetName("g2_xy");
    g2->SetTitle("Y vs X");
    g2->Write();

    // 7. Flatness color map if grid is regular

    std::vector<std::pair<double,double>> xy;
    xy.reserve(points.size());
    for (const auto &p : points)
        xy.emplace_back(p.coords[0], p.coords[1]);

    // Analyze (X, Y) points to determine if they form a regular Nx×Ny grid.
	// If yes, create a color-coded 2D histogram of Z values — the "flatness map".

    auto grid = GridFinder::analyze(xy);
    TH2D *hZ = nullptr;

    // Grid cells that actually received at least one point -- GridFinder
    // only requires the overall set of X's and set of Y's to be evenly
    // spaced (see GridFinder.h), NOT that every row/column is complete, so
    // a ragged/triangular scan shape is accepted and the Nx x Ny bounding
    // box below will contain cells nothing ever landed in. Declared here
    // (rather than inside the "if" block below, where the fill loop lives)
    // so both places that draw hZ -- the standalone map and the summary
    // slide -- can mark those cells the same way.
    std::set<std::pair<int,int>> measuredCells;

    // Paint every un-measured cell solid black, on top of whatever hZ's own
    // "COLZ" draw already put there. hZ leaves those cells at their default
    // content of 0, which (now that the map is colored by deviation, not
    // raw Z) sits in the middle of the real color scale -- indistinguishable
    // from "measured and genuinely flat" without this. Must be called with
    // the correct pad already current (cd()'d) and after hZ->Draw("COLZ").
    auto drawUnmeasuredMask = [&]() {
        if (!hZ) return;
        for (int ix = 0; ix < grid.Nx; ++ix) {
            for (int iy = 0; iy < grid.Ny; ++iy) {
                if (measuredCells.count({ix, iy})) continue;
                double x1 = grid.xMin + (ix - 0.5) * grid.dx;
                double x2 = grid.xMin + (ix + 0.5) * grid.dx;
                double y1 = grid.yMin + (iy - 0.5) * grid.dy;
                double y2 = grid.yMin + (iy + 0.5) * grid.dy;
                TBox *box = new TBox(x1, y1, x2, y2);
                box->SetFillColor(kBlack);
                box->SetFillStyle(1001);
                box->SetLineWidth(0);
                box->Draw("same");
            }
        }
    };

    if (grid.regularX && grid.regularY) {
        // Color-code by DEVIATION from the fitted plane (devVals, computed
        // above), not by the raw Z coordinate -- a "flatness map" should
        // show how far each point sits from the reference plane, and the
        // raw Z coordinate carries whatever baseline the CMM measured from
        // (which can be far from zero), unrelated to that.
        hZ = new TH2D("hZMap", "Flatness Map;X [mm];Y [mm];Deviation [mm]",
                      grid.Nx, grid.xMin - grid.dx/2, grid.xMax + grid.dx/2,
                      grid.Ny, grid.yMin - grid.dy/2, grid.yMax + grid.dy/2);

        std::map<std::pair<int,int>, std::vector<double>> bins;

        for (size_t i = 0; i < points.size(); ++i) {
            const auto& p = points[i];
            int ix = static_cast<int>(std::round((p.coords[0] - grid.xMin) / grid.dx));
            int iy = static_cast<int>(std::round((p.coords[1] - grid.yMin) / grid.dy));
            bins[{ix, iy}].push_back(devVals[i]);
            measuredCells.insert({ix, iy});
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

    std::string deviationsPdfPath, flatnessMapPdfPath;

    for (size_t i = 0; i < hists.size(); ++i) {
        std::string cname = "cHist_" + std::to_string(i + 1);
        TCanvas *c = new TCanvas(cname.c_str(), hists[i]->GetTitle(), 50 + i * 30, 50 + i * 30,
                                 canvasWidth, canvasHeight);
        c->Connect("Closed()", "TApplication", gApplication, "Terminate()");
        hists[i]->Draw();
        c->Update();

        // Save a PDF snapshot of the deviations-from-fit histogram.
        if (std::string(hists[i]->GetName()) == "hDeviations") {
            deviationsPdfPath = outDir + inputBase + "_deviations.pdf";
            c->SaveAs(deviationsPdfPath.c_str());
        }
    }

    TCanvas *c2 = new TCanvas("c2", "2D Scatter (Y vs X)", 900, 150, 700, 600);
    c2->Connect("Closed()", "TApplication", gApplication, "Terminate()");
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.8);
    g2->SetMarkerColor(kBlack);
    g2->Draw("AP");
    c2->Update();

    if (hZ) {


    double xrange = grid.xMax - grid.xMin;
	double yrange = grid.yMax - grid.yMin;

	int width = 800;
	int height = 800;

	if(xrange > yrange) height = static_cast<int>(width * yrange/xrange);
		else            width = static_cast<int>(height * xrange/yrange);

	TCanvas *cMap = new TCanvas("cMap","Flatness Map",
                            1650,150,width,height);



        //TCanvas *cMap = new TCanvas("cMap", "Flatness Map", 1650, 150, 800, 650);
        gStyle->SetPalette(kBird);
        cMap->SetLeftMargin(0.15);
        cMap->SetRightMargin(0.18);
        cMap->SetBottomMargin(0.12);
        cMap->SetTopMargin(0.08);
        hZ->SetStats(0);
        hZ->GetXaxis()->SetTitleOffset(1.2);
        hZ->GetYaxis()->SetTitleOffset(1.6);

        gPad->SetFixedAspectRatio();

        hZ->Draw("COLZ");
        drawUnmeasuredMask();
        cMap->Update();

        // Save a PDF snapshot of the 2D flatness map, including its color
        // legend (drawn automatically by the "COLZ" draw option above).
        flatnessMapPdfPath = outDir + inputBase + "_flatnessMap.pdf";
        cMap->SaveAs(flatnessMapPdfPath.c_str());

        outfile.cd();
    	cMap->Write("cMap");
	}

	//------------------------------------------------------------------------------
	// 8b. One-page summary slide -- the key results in a single image sized
	//     for a presentation (16:9), ready to drop into PowerPoint or Keynote:
	//       - best-fit plane residual RMS
	//       - GD&T ASME Y14.5 flatness result
	//       - deviations-from-fit histogram
	//       - 2D flatness color map with its legend
	//     Saved as a PNG (a plain raster image imports cleanly into both
	//     PowerPoint and Keynote with no conversion surprises) at a
	//     resolution well above screen size so text and plots stay crisp.
	//------------------------------------------------------------------------------
	{
		std::string summaryPngPath = outDir + inputBase + "_summary.png";

		TCanvas *cSummary = new TCanvas("cSummary", "FlatnessScan Summary", 1900, 100, 1920, 1080);
		cSummary->SetFillColor(kWhite);
		cSummary->cd();

		// Shared left margin: the title, the histogram, and the map all line
		// up against this same left edge.
		const double leftNdc = 0.03;

		// --- Title: a single left-aligned line along the very top. Text
		//     uses a "precision 3" ROOT font (the trailing "3" in 63/43
		//     below), which is a FIXED PIXEL SIZE independent of pad size --
		//     unlike the previous TPaveText auto-fit boxes, this makes the
		//     size directly and predictably tunable (just edit the pixel
		//     numbers below), which is what "reduce it by ~2 points" needs.
		//
		// Note: these strings are drawn on the ROOT canvas itself, not
		// printed to the terminal, so they use ROOT's own "#name" text-symbol
		// codes (#sigma, #mum, ...) instead of raw UTF-8 characters -- ROOT's
		// text rendering does not reliably interpret multi-byte UTF-8 and
		// mangles it. Plain std::cout output elsewhere in this program is
		// unaffected and keeps using the real UTF-8 characters, which render
		// correctly in a terminal.
		// Heap-allocated (not a stack local): cSummary stays open for
		// interactive viewing during the GUI event loop below and may get
		// repainted, so anything added to its pads must outlive this block --
		// same reasoning as every other ROOT primitive in this file, which is
		// why none of them are ever explicitly deleted.
		TLatex *latTitle = new TLatex();
		latTitle->SetNDC();
		latTitle->SetTextFont(63);   // Helvetica bold, fixed pixel size
		latTitle->SetTextSize(30);   // px -- reduced from the previous auto-fit size
		latTitle->SetTextAlign(13);  // left, top
		latTitle->DrawLatex(leftNdc, 0.985, ("FlatnessScan Summary  --  " + inputBase).c_str());

		// --- Deviations-from-fit histogram, left-aligned under the title. ---
		const double histX1 = leftNdc, histX2 = leftNdc + 0.40;
		const double histRowY1 = 0.46, histRowY2 = 0.93;

		TPad *padDev = new TPad("padDev", "padDev", histX1, histRowY1, histX2, histRowY2);
		padDev->SetFillColor(kWhite);
		padDev->Draw();
		padDev->cd();
		for (auto h : hists) {
			if (std::string(h->GetName()) == "hDeviations") {
				h->Draw();
				break;
			}
		}

		cSummary->cd();

		// --- Results text, immediately to the right of the histogram, in a
		//     smaller fixed pixel size than the title. ---
		std::ostringstream ossRms, ossFlat;
		ossRms.setf(std::ios::fixed);   ossRms   << std::setprecision(2);
		ossFlat.setf(std::ios::fixed);  ossFlat  << std::setprecision(2);
		ossRms  << "Best-fit plane residual RMS (#sigma):  " << (1000.0 * sqrt(minChi2 / X.size())) << " #mum";
		ossFlat << "GD&T ASME Y14.5 Flatness:  " << (flatness * 1000.0) << " #mum";

		const double resultsX = histX2 + 0.03;

		TLatex *latRms = new TLatex();
		latRms->SetNDC();
		latRms->SetTextFont(43);   // Helvetica regular, fixed pixel size
		latRms->SetTextSize(20);  // px -- smaller than the title
		latRms->SetTextAlign(12); // left, vcenter
		latRms->DrawLatex(resultsX, 0.72, ossRms.str().c_str());

		TLatex *latFlat = new TLatex();
		latFlat->SetNDC();
		latFlat->SetTextFont(43);
		latFlat->SetTextSize(20);
		latFlat->SetTextAlign(12);
		latFlat->DrawLatex(resultsX, 0.64, ossFlat.str().c_str());

		cSummary->cd();

		// --- 2D flatness color map with its legend, stacked below the
		//     histogram and left-aligned to the same left edge. Unlike the
		//     histogram, this plot's aspect ratio carries real information
		//     (it mirrors the physical X/Y extent of the scanned surface),
		//     so instead of just filling its allotted region it is sized to
		//     match the data's (xrange / yrange) aspect ratio -- exactly
		//     like the standalone flatness-map PDF above -- then placed flush
		//     against the left border and vertically centered within the
		//     region below the histogram. Only produced when the points form
		//     a regular grid. The freed-up vertical space from shrinking the
		//     title block above goes entirely into this region and the
		//     histogram row above, growing both by the same proportion so
		//     their relative sizing and alignment stay the same as before.
		if (hZ) {
			const double canvasWpx = 1920.0, canvasHpx = 1080.0;
			// Region reserved for the map: left-aligned, bottom 46% of the canvas.
			const double regionY1Ndc = 0.0, regionY2Ndc = 0.46;
			const double regionWpx = canvasWpx;  // available width to fit within
			const double regionHpx = (regionY2Ndc - regionY1Ndc) * canvasHpx;

			double mapXrange = grid.xMax - grid.xMin;
			double mapYrange = grid.yMax - grid.yMin;
			double dataAspect = mapXrange / mapYrange;  // width / height

			// "Contain" fit: shrink to the region's width or height,
			// whichever binds first, preserving dataAspect exactly.
			const double fitMargin = 0.92;  // small margin so it doesn't touch the region's edges
			double availWpx = regionWpx * fitMargin;
			double availHpx = regionHpx * fitMargin;
			double mapWpx, mapHpx;
			if (availWpx / availHpx > dataAspect) {
				mapHpx = availHpx;
				mapWpx = availHpx * dataAspect;
			} else {
				mapWpx = availWpx;
				mapHpx = availWpx / dataAspect;
			}

			// Flush against the left border (not centered), vertically
			// centered within the reserved region, then converted back to
			// NDC coordinates on the full canvas.
			double mapX1px = leftNdc * canvasWpx;
			double mapY1px = regionY1Ndc * canvasHpx + (regionHpx - mapHpx) / 2.0;

			double padX1 = mapX1px / canvasWpx;
			double padX2 = (mapX1px + mapWpx) / canvasWpx;
			double padY1 = mapY1px / canvasHpx;
			double padY2 = (mapY1px + mapHpx) / canvasHpx;

			TPad *padMap = new TPad("padMap", "padMap", padX1, padY1, padX2, padY2);
			padMap->SetFillColor(kWhite);
			padMap->Draw();
			padMap->cd();
			padMap->SetLeftMargin(0.15);
			padMap->SetRightMargin(0.18);
			padMap->SetBottomMargin(0.12);
			padMap->SetTopMargin(0.08);
			gPad->SetFixedAspectRatio();
			hZ->Draw("COLZ");
			drawUnmeasuredMask();
		} else {
			cSummary->cd();
			TLatex *latNoMap = new TLatex();
			latNoMap->SetNDC();
			latNoMap->SetTextFont(43);
			latNoMap->SetTextSize(20);
			latNoMap->SetTextAlign(12);
			latNoMap->DrawLatex(leftNdc, 0.23, "No regular grid -- flatness map unavailable");
		}

		cSummary->cd();
		cSummary->Update();
		cSummary->SaveAs(summaryPngPath.c_str());
		std::cout << "Summary slide (PNG) written to " << summaryPngPath << std::endl;
	}

	//------------------------------------------------------------------------------
	// 9. Run ROOT GUI loop
	//------------------------------------------------------------------------------
	//
	// Close the output file before entering the interactive ROOT GUI loop.
	// Canvases remain accessible even after the file is closed.
	//


	std::cout << "\nHistograms written to " << outname << std::endl;
	std::cout << "Log written to " << logname << std::endl;
	if (!deviationsPdfPath.empty())
		std::cout << "Deviations histogram (PDF) written to " << deviationsPdfPath << std::endl;
	if (!flatnessMapPdfPath.empty())
		std::cout << "Flatness map (PDF) written to " << flatnessMapPdfPath << std::endl;
	std::cout << "\nHit ctrl-c to exit" << std:: endl;



	// detach histograms from file so they survive after outfile.Close()
	for (TH1D* hist : hists) {
		hist->SetDirectory(nullptr);
		//hist->Write();
	}
	// detach and write scatter plot
	g2->Write();

	// write flatness map if present
	if (hZ) {
		hZ->SetDirectory(nullptr);
		//hZ->Write();
	}
		outfile.Close();

	// Enter the ROOT GUI event loop — close all canvases or press Ctrl+C to exit.

	app.Run();

	// Restore the original stream buffers (reached only if all canvases are
	// closed normally rather than via Ctrl+C).
	std::cout.rdbuf(origCoutBuf);
	std::cerr.rdbuf(origCerrBuf);

    return 0;
}

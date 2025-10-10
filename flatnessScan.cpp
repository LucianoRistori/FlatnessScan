
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
//     7. Writes all histograms and the TGraph to an output ROOT file
//        ("output.root").
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
//   - ROOT framework (TFile, TH1D, TGraph, TCanvas, TApplication, Minimizer)
//   - Points.h / Points.cpp for reading input data
//
// Author: Luciano Ristori
// Version: 1.0
// Date: October 2025
//------------------------------------------------------------------------------

// This program depends on the ROOT framework for histogramming, fitting,
// and visualization. It also uses the 'Points' module for reading 3D data files.

//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <limits>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TGraph.h"


#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "Points.h"

// This program depends on the ROOT framework for histogramming, fitting,
// and visualization. It also uses the 'Points' module for reading 3D data files.


using std::cout;
using  std::endl;
using  std::vector;


// Global vectors storing coordinate components, used by the minimization function.
// These are filled in main() after reading the input points.


//std::vector<Point> points;
vector<double> X, Y, Z;

double offset = 400.;

////////////////////////////////////////////////////////////////////////////////////
// RAII helper class to temporarily set scientific notation and precision for console output.

class ScientificPrecision {

public:

	// Constructor saves current format and applies scientific formatting.
    ScientificPrecision(std::ostream& os, int precision)
        : os_(os), old_precision_(os.precision()), old_flags_(os.flags()) {
        os_ << std::scientific << std::setprecision(precision);
    }
    
    // Destructor restores original stream settings automatically.
    ~ScientificPrecision() {
        os_.precision(old_precision_);
        os_.flags(old_flags_);
    }

    ScientificPrecision(const ScientificPrecision&) = delete;
    ScientificPrecision& operator=(const ScientificPrecision&) = delete;

private:
    std::ostream& os_;
    std::streamsize old_precision_;
    std::ios_base::fmtflags old_flags_;
};
/*

// Usage example

int main() {
    double value = 3.14159265359;

    std::cout << "Default notation: " << value << std::endl;

    {
        ScientificPrecision sp(std::cout, 4);
        std::cout << "Scientific notation with precision 4: " << value << std::endl;
    }

    std::cout << "Default notation restored: " << value << std::endl;

    return 0;
}
*/

////////////////////////////////////////////////////////////////////////////////////
// RAII helper class to temporarily set floating point precision for console output.

class FloatingPointPrecision {

public:

	// Constructor saves current format and applies floating point precision. 
    FloatingPointPrecision(std::ostream& os, int precision)
        : os_(os), old_precision_(os.precision()) {
        os_ << std::setprecision(precision);
    }
    
    // Destructor restores original stream settings automatically.
    ~FloatingPointPrecision() {
        os_.precision(old_precision_);
    }

    FloatingPointPrecision(const FloatingPointPrecision&) = delete;
    FloatingPointPrecision& operator=(const FloatingPointPrecision&) = delete;

private:
    std::ostream& os_;
    std::streamsize old_precision_;
};

/*

// usage example

int main() {
    double value = 3.14159265359;

    std::cout << "Default precision: " << value << std::endl;

    {
        FloatingPointPrecision fpp(std::cout, 4);
        std::cout << "Precision 4: " << value << std::endl;
    }

    std::cout << "Default precision restored: " << value << std::endl;

    return 0;
}
*/


////////////////////////////////////////////////////////////////////////////////////
// chi2Func()
// Computes the chi-square value for the given plane parameters (ax, ay, az).
// The function is minimized by the ROOT Minuit2 optimizer to find the best-fit plane.

double chi2Func(const double *x){
	 	   	       	
	 double ax = x[0]; 
	 double ay = x[1];
	 double az = x[2]; 
	 
	double chi2 = 0.;
	
	// For each point, compute the signed distance from the trial plane and accumulate chi².

	for(int i = 0; i != (int)X.size(); ++i) {	
		double delta = ax*X[i]+ay*Y[i]+az*(Z[i]+offset) - 1.;
		chi2 += delta*delta/(ax*ax+ay*ay+az*az);
	}
	
	return chi2;

}// end chi2Func

////////////////////////////////////////////////////////////////////////////////////




int main(int argc, char *argv[]) {


// 1. Parse command-line arguments and initialize ROOT application.
 
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.csv" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

 	TApplication app("app", &argc, argv);  // Initialize ROOT GUI 	
    gROOT->SetBatch(false); // enable GUI
    TH1::AddDirectory(kTRUE);


// 2. Read 3D points from input file (see Points.cpp for parser details).

    int n = 4; // >= 4 numbers expected on each line
    
    std::vector<Point> points = readPoints(filename, n);

    if (points.empty()) {
        std::cerr << "No valid points found. Exiting." << std::endl;
        return 1;
    }
    
// 3. Fill global coordinate vectors (X, Y, Z) used in the fit function.
  
    for(int iPoint = 0; iPoint != (int)points.size(); ++iPoint){
    	X.push_back(points[iPoint].coords[1]);
    	Y.push_back(points[iPoint].coords[2]);
    	Z.push_back(points[iPoint].coords[3]);   
    }
    std::cout << "Read " << points.size() << " valid points with "
              << n << " coordinates each." << std::endl;
              
              

// 4. Configure and run ROOT's Minuit2 minimizer to fit a 3D plane.

	cout << "\n Fitting 3D plane..." << endl;

	// create minimizer giving a name and a name (optionally) for the specific
	// algorithm
	// possible choices are: 
	//     minName                  algoName
	// Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
	//  Minuit2                     Fumili2
	//  Fumili
	//  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS, 
	//                              BFGS2, SteepestDescent
	//  GSLMultiFit
	//   GSLSimAn xxx problems
	//   Genetic

	const char * minName = "Minuit2";
	const char *algoName = "";

	ROOT::Math::Minimizer* min = 
	  ROOT::Math::Factory::CreateMinimizer(minName, algoName);
	
		// set tolerance , etc...
		min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
		min->SetMaxIterations(1000);  // for GSL 
		min->SetTolerance(0.001);
		min->SetPrintLevel(0);

		// create function wrapper for minmizer
		// a IMultiGenFunction type 
		ROOT::Math::Functor f(&chi2Func,3); 
		double step[3] = {0.001,0.001,0.001};
		
		
		double variable[3] = {0.,0.,1./offset};
		
		
		 min->SetFunction(f);
					
		 
		// Set the free variables to be minimized!

		min->SetVariable(0,"ax",variable[0], step[0]);
		min->SetVariable(1,"ay",variable[1], step[1]);
		min->SetVariable(2,"az",variable[2], step[2]);
		
		 // do the minimization 
				
		min->Minimize(); 
		
// 5. Retrieve fit results and compute derived quantities (σ, |a|, etc.).

		const double *res = min->X();
		const double *err = min->Errors();
				
		// return parameters f
		
		double minChi2 = min->MinValue(); 
		double ax = res[0];
		double ay = res[1];
		double az = res[2];
		
		double ax_e = err[0];
		double ay_e = err[1];
		double az_e = err[2];
		
		// print results
		
		{
		ScientificPrecision sp(cout, 2);	
		
		std::cout << "\n----------------------------------\n";
		std::cout << "  Plane fit summary\n";
		std::cout << "  ax = " << ax << " ± " << ax_e << "\n";
		std::cout << "  ay = " << ay << " ± " << ay_e << "\n";
		std::cout << "  az = " << az << " ± " << az_e << "\n";
		}
		
		{
        FloatingPointPrecision fpp(std::cout, 4);
		std::cout << "  σ = " << 1000.*sqrt(minChi2 / X.size()) << " µm\n";
		std::cout << "----------------------------------\n";
		}
		
		// Modulus of a
		
		double moda = sqrt(ax*ax+ay*ay+az*az);
		
		cout << "\nModulus of a: " << moda << endl;
		
		// Inverse modulus
		
		double invModa = 1./moda;
		
		cout << "Inverse modulus of a: " << invModa << " [mm]" << endl;
		
		cout << "Offset: " << offset << " [mm]" << endl;
		
		
//////////////////////////////////////////////////////////////////////////////////////				

// 6. Determine coordinate ranges for histogram axis scaling.

    std::vector<double> mins(n, std::numeric_limits<double>::max());
    std::vector<double> maxs(n, std::numeric_limits<double>::lowest());

    for (const auto &p : points) {
        for (int i = 0; i < n; ++i) {
            if (p.coords[i] < mins[i]) mins[i] = p.coords[i];
            if (p.coords[i] > maxs[i]) maxs[i] = p.coords[i];
        }
    }
    
    
// 7. Create and fill histograms for X, Y, Z, and residuals from plane fit.

    // Create output ROOT file
    
    std::string outname = "output.root";
    TFile outfile(outname.c_str(), "RECREATE");

    // Automatically determine histogram binning and range for each coordinate.
	// Adds a 50% margin to ensure all points fall comfortably inside the histogram range.
    
    std::vector<TH1D*> hists;
    for (int i = 0; i < n; ++i) {
        double min = mins[i];
        double max = maxs[i];
        if (min == max) {
            min -= 0.5;
            max += 0.5;
        }
        double margin = 0.5 * (max - min); // 50% margin
        int nBins = static_cast<int>((max - min + 2*margin) * 1000 + 0.5);  // round to nearest int
        std::string hname = "hCoord" + std::to_string(i+1);
        std::string htitle = "Coordinate " + std::to_string(i+1);
        hists.push_back(new TH1D(hname.c_str(), htitle.c_str(),
                                 nBins, min - margin, max + margin));
        if(i == 3){// this is Z
        	hname = "Deviations";
        	htitle = "Deviations from 3D plane fit";
        	hists.push_back(new TH1D(hname.c_str(), htitle.c_str(),
                                 nBins, min - margin, max + margin));                               
        }
    }
    
	cout << "nHists = " << hists.size() << endl;

    // Fill histograms
    for (const auto &p : points) {
        for (int i = 0; i < n; ++i) {
            hists[i]->Fill(p.coords[i]);
        }
        // Compute perpendicular deviation (signed distance) of each point from the best-fit plane.
		double delta = (ax*p.coords[1]+ay*p.coords[2]+az*(p.coords[3]+offset) -1.)*invModa;
        hists[4]->Fill(delta);
    }

    // Write to file
    for (auto h : hists) {
        h->Write();
    }
    

///////////////////////////////////////
/// Display histograms
//////////////////////////////////////

	// One canvas per 1D histogram
	
	n = (int)hists.size();
	
	std::vector<TCanvas*> canvases;
	int canvasWidth = 800;  
	int canvasHeight = 600;  

	for (int i = 0; i < n; ++i) {
		std::string cname = "c1D_" + std::to_string(i+1);
		std::string ctitle = "Histogram of coordinate " + std::to_string(i+1);

		int xPos = 50 + (i * 40);
		int yPos = 50 + (i * 40);

		TCanvas* c = new TCanvas(cname.c_str(), ctitle.c_str(),
								 xPos, yPos, canvasWidth, canvasHeight);
		
		c->Connect("Closed()", "TApplication", gApplication, "Terminate()");

		hists[i]->Draw();
		c->Update();

		canvases.push_back(c);
	}
	
// 8. Create a 2D scatter plot of Y vs. X for visualization.
//      2D scatter plot using TGraph (safe for many points)
//
	// Make a TGraph with one point per entry
	TGraph* g2 = new TGraph(points.size());
	for (size_t i = 0; i < points.size(); ++i) {
		if (points[i].coords.size() >= 3) {  // safety check
			g2->SetPoint(i, points[i].coords[1], points[i].coords[2]);
		}
	}

	// Optional: associate with current ROOT file for saving
	g2->SetName("g2_xy");
	g2->SetTitle("coords[2] vs coords[1]");
	g2->Write(); // write to file
	
		
// 9. Display histograms and scatter plots in ROOT canvases, and write all objects to file.

	// Display in a separate canvas
	TCanvas* c2 = new TCanvas("c2", "2D Scatter", 600, 600);
	c2->Connect("Closed()", "TApplication", gApplication, "Terminate()"); 
	c2->SetWindowPosition(200, 200);

	// Marker settings
	g2->SetMarkerStyle(20);  // filled circle
	g2->SetMarkerSize(0.8);
	g2->SetMarkerColor(kBlack);

	// Draw the points with axes
	g2->Draw("AP"); // A = draw axes, P = points
	c2->Update();

// 10. Keep GUI open until user closes the ROOT windows.
//    	actually you need to cntrl-c

    app.Run();  
    
    
// End of program. Histograms and graphs are written to 'output.root' for later inspection.
   
    outfile.Close();
    std::cout << "Histograms written to " << outname << std::endl;

    return 0;
    
}// end of main

// End of program. Histograms and graphs are written to 'output.root' for later inspection.



/*
 * testGridFinder.cpp
 *
 * ------------------------------------------------------------
 * Purpose:
 *     Standalone test and demonstration program for GridFinder.h.
 *     Reads (n, X, Y, Z) data from a text or CSV file using
 *     readPoints() and checks whether the (X,Y) coordinates
 *     correspond to a regular rectangular grid.
 *
 * Workflow:
 *     1. Input file lines must contain four numeric values:
 *            n, X, Y, Z
 *        where n is a point index (ignored in analysis).
 *
 *     2. The program calls readPoints(filename, 4)
 *        to read all points into a std::vector<Point>.
 *
 *     3. It extracts the second and third coordinates
 *        (X = coords[1], Y = coords[2]) into a vector of pairs (X,Y).
 *
 *     4. The vector is passed to:
 *            GridFinder::analyze(xy, tolFrac, presenceFrac, mergeFrac)
 *
 *        which returns a GridFinder::Result containing:
 *            Nx, Ny, dx, dy, regularX/Y, missingPoints
 *
 * Command line usage:
 *        ./testGridFinder <pointsFile>
 *
 * Example:
 *        ./testGridFinder ASTRAL_GRANITE_VISION_FLATNESS_004.csv
 *
 * Example output:
 *        Analyzed file: ASTRAL_GRANITE_VISION_FLATNESS_004.csv
 *        Number of points: 100
 *        Nx=10 Ny=10  dx=31.1111  dy=19.4444
 *        Regular X: true   Regular Y: true
 *        Missing grid points: 0
 *
 * Notes:
 *     - Lines containing text headers or wrong column counts
 *       are skipped automatically by readPoints().
 *     - GridFinder tolerates small floating-point differences
 *       (e.g., 29.444 vs 29.445) via mergeStepFraction.
 *     - The output provides a quick diagnostic of grid
 *       completeness and uniformity before deeper analysis
 *       (e.g., flatness or Z-based surface studies).
 *
 * ------------------------------------------------------------
 * Author:    Luciano Ristori (with ChatGPT assistance)
 * Created:   Oct 2025
 * License:   MIT / open use
 * ------------------------------------------------------------
 */

#include <iostream>
#include <vector>
#include <utility>
#include "Points.h"
#include "GridFinder.h"

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <pointsFile>" << std::endl;
        return 1;
    }

    std::string fileName = argv[1];
    int nCols = 4;  // each line has: n, X, Y, Z

    // --- Read all points ---
    std::vector<Point> points = readPoints(fileName, nCols);
    if (points.empty()) {
        std::cerr << "Error: no points read from file " << fileName << std::endl;
        return 1;
    }

    // --- Extract (X,Y) pairs from coords[1] and coords[2] ---
    std::vector<std::pair<double,double>> xy;
    xy.reserve(points.size());
    for (const auto &p : points) {
        if (p.coords.size() >= 3)
            xy.emplace_back(p.coords[1], p.coords[2]);
    }

    // --- Analyze grid ---
    GridFinder::Result res = GridFinder::analyze(xy, 0.05, 0.2, 0.10);

    // --- Print results ---
    std::cout << "Analyzed file: " << fileName << std::endl;
    std::cout << "Number of points: " << points.size() << std::endl;
    std::cout << "Nx=" << res.Nx << " Ny=" << res.Ny
              << "  dx=" << res.dx << "  dy=" << res.dy << std::endl;
    std::cout << "Regular X: " << std::boolalpha << res.regularX
              << "   Regular Y: " << res.regularY << std::endl;
    std::cout << "Missing grid points: " << res.missingPoints << std::endl;

    return 0;
}

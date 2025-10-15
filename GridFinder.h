/*
 * GridFinder.h
 *
 * ------------------------------------------------------------
 * Purpose:
 *     Provides utilities to verify whether a collection of (X,Y)
 *     points lies on a regular rectangular grid.
 *
 * Overview:
 *     - Determines the number of distinct X and Y coordinates (Nx, Ny)
 *     - Computes average grid steps (dx, dy)
 *     - Checks grid regularity (uniform spacing within tolerance)
 *     - Estimates how many grid points are missing
 *
 * Usage:
 *     #include "GridFinder.h"
 *
 *     std::vector<std::pair<double,double>> xy = { {10,10}, {10,30}, {40,10}, {40,30} };
 *     GridFinder::Result res = GridFinder::analyze(xy);
 *
 *     if (res.regularX && res.regularY)
 *         std::cout << "Regular " << res.nx << "x" << res.ny << " grid with dx="
 *                   << res.dx << " dy=" << res.dy << std::endl;
 *
 * Data model:
 *     The input is a vector of (X,Y) pairs.  The algorithm:
 *       1. Extracts all unique X and Y coordinates.
 *       2. Merges coordinates that differ by less than a
 *          user-defined fraction of the grid step (mergeStepFraction).
 *       3. Computes mean and spread of consecutive spacings.
 *       4. Declares the grid "regular" if the fractional spread is below
 *          toleranceFraction (default 5%).
 *       5. Scans all (X,Y) combinations to detect missing grid points.
 *
 * Adjustable parameters:
 *     toleranceFraction       - allowed deviation from uniform spacing
 *     presenceEpsilonFraction - proximity threshold for missing points
 *     mergeStepFraction       - fractional step difference used to merge
 *                               nearly identical coordinates
 *
 * Typical defaults:
 *     toleranceFraction       = 0.05   (5% spacing uniformity)
 *     presenceEpsilonFraction = 0.2    (20% distance tolerance)
 *     mergeStepFraction       = 0.10   (10% merge tolerance)
 *
 * Output:
 *     Returned GridFinder::Result structure contains:
 *         nx, ny          → number of distinct grid positions in X and Y
 *         dx, dy          → average step size
 *         regularX, regularY → spacing uniformity flags
 *         missingPoints   → count of missing grid intersections
 *
 * Example:
 *     // A perfect 10x10 grid with minor rounding noise
 *     std::vector<std::pair<double,double>> xy = loadPoints(...);
 *     auto res = GridFinder::analyze(xy, 0.05, 0.2, 0.1);
 *     // → Nx=10 Ny=10 dx≈31.11 dy≈19.44 regularX=Y=true missing=0
 *
 * ------------------------------------------------------------
 * Author:    Luciano Ristori (with ChatGPT assistance)
 * Created:   Oct 2025
 * License:   MIT / open use
 * ------------------------------------------------------------
 */
 


#ifndef GRID_FINDER_H
#define GRID_FINDER_H

#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace GridFinder {

struct Result {
    bool regularX = false;
    bool regularY = false;
    double dx = 0.0;
    double dy = 0.0;
    int Nx = 0;
    int Ny = 0;
    double xMin = 0.0;
    double xMax = 0.0;
    double yMin = 0.0;
    double yMax = 0.0;
    int missingPoints = 0;
};

// ------------------------------------------------------------
// analyze()
//   Checks if a set of (X,Y) points form a regular rectangular grid.
// ------------------------------------------------------------
inline Result analyze(
    const std::vector<std::pair<double,double>> &points,
    double toleranceFraction = 0.05,        // allowed deviation in spacing (~5%)
    double presenceEpsilonFraction = 0.2,   // proximity for missing-point detection
    double mergeStepFraction = 0.10         // merge threshold = 10% of mean step
)
{
    Result result;
    if (points.size() < 4) return result;

    // --- Extract X and Y coordinates ---
    std::vector<double> xs, ys;
    xs.reserve(points.size());
    ys.reserve(points.size());
    for (auto &p : points) {
        xs.push_back(p.first);
        ys.push_back(p.second);
    }
    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end());

    // --- Estimate average step for merge tolerance ---
    auto estimateMeanStep = [](const std::vector<double> &v) {
        if (v.size() < 2) return 0.0;
        std::vector<double> d;
        for (size_t i = 1; i < v.size(); ++i)
            d.push_back(v[i] - v[i-1]);
        return std::accumulate(d.begin(), d.end(), 0.0) / d.size();
    };
    double dxGuess = estimateMeanStep(xs);
    double dyGuess = estimateMeanStep(ys);

    // --- Merge close coordinates ---
    auto mergeClose = [&](std::vector<double> &v, double step) {
        std::vector<double> merged;
        if (v.empty()) return merged;
        merged.push_back(v[0]);
        double eps = std::fabs(step) * mergeStepFraction;
        for (double val : v) {
            if (std::fabs(val - merged.back()) > eps)
                merged.push_back(val);
        }
        return merged;
    };
    xs = mergeClose(xs, dxGuess);
    ys = mergeClose(ys, dyGuess);

	result.Nx = xs.size();
	result.Ny = ys.size();
	if (result.Nx < 2 || result.Ny < 2) return result;

    // --- Compute step statistics ---
    auto spacingStats = [](const std::vector<double> &v) {
        std::vector<double> d;
        for (size_t i = 1; i < v.size(); ++i)
            d.push_back(v[i] - v[i-1]);
        double mean = std::accumulate(d.begin(), d.end(), 0.0) / d.size();
        double spread = 0.0;
        for (double di : d)
            spread += std::fabs(di - mean);
        return std::make_pair(mean, spread / d.size());
    };

    auto [dx, dxSpread] = spacingStats(xs);
    auto [dy, dySpread] = spacingStats(ys);
    result.dx = dx;
    result.dy = dy;

    result.regularX = (dxSpread / dx < toleranceFraction);
    result.regularY = (dySpread / dy < toleranceFraction);

    // --- Count missing grid points ---
    double eps = std::min(dx, dy) * presenceEpsilonFraction;
    size_t missing = 0;
    for (double x : xs) {
        for (double y : ys) {
            bool found = false;
            for (auto &p : points) {
                if (std::fabs(p.first - x) < eps &&
                    std::fabs(p.second - y) < eps) {
                    found = true;
                    break;
                }
            }
            if (!found) ++missing;
        }
    }
    result.missingPoints = missing;
    
    result.xMin = *std::min_element(xs.begin(), xs.end());
	result.xMax = *std::max_element(xs.begin(), xs.end());
	result.yMin = *std::min_element(ys.begin(), ys.end());
	result.yMax = *std::max_element(ys.begin(), ys.end());
	result.Nx = std::round((result.xMax - result.xMin) / result.dx) + 1;
	result.Ny = std::round((result.yMax - result.yMin) / result.dy) + 1;


    return result;
}

} // namespace GridFinder

#endif // GRID_FINDER_H

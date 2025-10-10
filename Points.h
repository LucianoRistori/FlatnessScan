#ifndef POINTS_H
#define POINTS_H
/*
 * Points.h
 *
 * Provides the Point struct and readPoints() function for reading N-dimensional coordinates
 * from a text or CSV file. Designed to be reusable in multiple programs.
 *
 * -------------------------------------------
 * Workflow / Usage:
 *
 * 1. Include this header and compile with Points.cpp:
 *
 *      #include "Points.h"
 *      g++ -std=c++17 -O2 -Wall main.cpp Points.cpp `root-config --cflags --libs` -o myProgram
 *
 * 2. Read points from a file:
 *
 *      int n = 3; // number of coordinates per point
 *      std::vector<Point> points = readPoints("points.csv", n);
 *
 * 3. Process points:
 *
 *      for (const auto& p : points) {
 *          for (int i = 0; i < n; ++i) {
 *              std::cout << p.coords[i] << " ";
 *          }
 *          std::cout << std::endl;
 *      }
 *
 * 4. Optional: create histograms with ROOT (1D per coordinate) or a 2D scatter plot:
 *
 *      // Using makeHists.cpp workflow:
 *      ./makeHists points.csv 3
 *
 *      This will:
 *          - Read valid points from the file
 *          - Create one TH1D histogram per coordinate
 *          - Create a TGraph of coords[1] vs coords[0]
 *          - Save all objects to output.root
 *          - Display each plot in its own GUI window
 *
 * -------------------------------------------
 * File format:
 *
 *  - Text or CSV file with one point per line
 *  - Coordinates separated by spaces or commas
 *  - Lines with fewer numbers than expected (n) are skipped
 *
 * Example file (3 coordinates per line):
 *
 *   1.23, 4.56, 7.89
 *   2.34, 5.67, 8.90
 *   3.45, 6.78, 9.01
 *
 * Example minimal C++ program:
 *
 *   #include "Points.h"
 *   #include <iostream>
 *
 *   int main() {
 *       int n = 3;
 *       auto points = readPoints("points.csv", n);
 *       for (const auto& p : points) {
 *           for (int i = 0; i < n; ++i) std::cout << p.coords[i] << " ";
 *           std::cout << std::endl;
 *       }
 *   }
 *
 */



#include <string>
#include <vector>

// A simple struct to hold one point with n coordinates
struct Point {
    std::vector<double> coords;
    Point(int n) : coords(n, 0.0) {}
};

// Reads a file of points, each line containing n coordinates.
// Returns a vector of Point objects.
std::vector<Point> readPoints(const std::string &filename, int n);

#endif // POINTS_H

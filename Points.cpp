//------------------------------------------------------------------------------
// File: Points.cpp
//
// Description:
//   Utility function to read 3D (or N-dimensional) points from a text or CSV file.
//   Each line is expected to contain at least 'n' numerical values.
//   Lines with fewer than 'n' values are skipped with a warning.
//   Both space- and comma-separated files are supported.
//
// Behavior:
//   - On missing or unreadable file → exits with error message.
//   - On malformed lines → prints warning and skips the line.
//   - Returns a vector<Point> with 'n' coordinates per point.
//
//------------------------------------------------------------------------------
// more info in Points.h
//------------------------------------------------------------------------------

#include "Points.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

std::vector<Point> readPoints(const std::string &filename, int n) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        exit(1);
    }

    std::vector<Point> points;
    std::string line;
    int lineNum = 0;
    
	// begin loop that reads lines///////////////////////
	
    while (std::getline(infile, line)) {
        ++lineNum;

        // Replace commas with spaces to allow reading 
        //both CSV and space-separated formats.
        
        for (char &c : line) {
            if (c == ',') c = ' ';
        }

        std::istringstream iss(line);
        std::vector<double> numbers;
        double val;

        while (iss >> val) {
            numbers.push_back(val);
        }
        
        // Skip lines that do not have enough numbers (e.g. missing coordinates).
		// This avoids exceptions and allows the program to continue with valid data only.

        if (numbers.size() < static_cast<size_t>(n)) {
            std::cerr << "Warning: line " << lineNum
                      << " has only " << numbers.size()
                      << " numbers, expected " << n << ". Skipping." << std::endl;
            continue;
        }
        
        // Create a Point object with n coordinates and copy the parsed values into it.
        Point p(n);
        for (int i = 0; i < n; ++i) {
            p.coords[i] = numbers[i];
        }
        points.push_back(p);
        
    }// end loop on lines
    
    // Return the vector of successfully read points.
	// If no valid lines were found, the vector will be empty.

    return points;
}

# FlatnessScan

**FlatnessScan** is a C++17 program for analyzing surface flatness from 3D measurement data.  
It reads a CSV-like file of point coordinates, fits an optimal plane to the data, evaluates deviations,  
and visualizes the results as histograms and color maps using the ROOT framework.

---

## ✳️ Overview

The program performs the following tasks:

1. **Reads** a text file with four columns:  
   `index, X, Y, Z` — where `X, Y, Z` are coordinates in millimeters.
   Additional numbers on the same line are ignored.
   Any text is ignored
2. **Fits** a best-fit plane using ROOT’s Minuit2 minimizer.
3. **Computes** the flatness standard deviation and plane coefficients.
4. **Produces**:
   - 1D histograms for X, Y, Z, and residuals from the plane fit.  
   - A 2D scatter plot of Y vs X.  
   - A **colored flatness map** (Z as color) if the (X,Y) points form a regular grid.
5. **Saves** all histograms and graphs to a ROOT file.

---

## 🧩 Input Format

Example of a valid input file:
(only the first four columns are actually used)

STATS
PLANK SCAN
07/16/25 14:39:46 
Layer3ScanTopV1.RTN

1,10,10,0.001,0.004,0.009,1
2,41.111,10,0.006,0.004,0.009,1
3,72.222,10,0.024,0.004,0.009,1
4,103.333,10,-0.006,0.004,0.009,1
5,134.444,10,-0.03,0.004,0.009,1
6,165.556,10,-0.01,0.004,0.009,1
7,196.667,10,0.003,0.004,0.009,1
8,227.778,10,-0.011,0.004,0.009,1
9,258.889,10,0.033,0.004,0.009,1
.......

# FlatnessScan v2.1.0

**Author:** Luciano Ristori  
**Date:** October 2025  

---

## Overview

`FlatnessScan` analyzes CMM (Coordinate Measurement Machine) surface measurements to determine the flatness of a scanned surface.  
It fits a plane to the measured points, computes residuals (ΔZ, ΔR), and generates histograms and 2D color maps using ROOT.

This version introduces integration with the shared `common` module, which now provides unified point handling and optional label support.

---

## Features

- Reads 3D point data files (`X Y Z` or `label X Y Z`, with optional CSV format)
- Fits a plane using Minuit2 minimization
- Computes residuals and flatness statistics (mean, sigma)
- Produces ROOT histograms and 2D color maps
- Compatible with labeled point data via `common v1.2.1`

---

## Dependencies

- **common module:** v1.2.1  
  Provides the `Points` structure and `readPoints()` function with label and CSV support  
  (Tag: `v1.2.1` in the `common` repository)

- **ROOT Framework:** v6.30+  
  Required for histogramming and visualization  
  [https://root.cern](https://root.cern)

- **C++17 compiler (clang++)**

---

## Build Instructions

```bash
make clean
make

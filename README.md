# FlatnessScan

**Version:** 1.1.0 — *October 2025*  
**Author:** Luciano Ristori  

---

## 🧭 Overview

**FlatnessScan** is a ROOT-based analysis tool for evaluating the flatness of measured 3D surfaces.  
It reads point coordinate data, fits an optimal plane using ROOT’s Minuit2 minimizer, and produces detailed histograms, residuals, and 2D visualizations of the surface flatness.

The program was designed to help characterize and visualize mechanical surface deviations with micrometer precision.

---

## 🚀 Features

- **Robust plane fitting**
  - Uses ROOT’s Minuit2 minimizer to determine optimal plane coefficients.
  - Reports χ², standard deviation, and plane-normal normalization.

- **Improved histogram handling**
  - Histograms remain interactive and stable after the ROOT file is closed.
  - Fully functional on macOS and Linux.

- **Automatic grid detection**
  - Integrated `GridFinder` module identifies regular Nx×Ny point grids.
  - Generates a color-coded 2D “Flatness Map” when the grid is regular.

- **Flexible output**
  - Optional output file name (adds `.root` automatically if missing).
  - Saves all histograms, scatter plots, and metadata tags.

- **Readable, modular source code**
  - Extensive comments and clear separation between reading, fitting, plotting, and saving.

---

## 🧪 Example Usage

```bash
./flatnessScan my_points.csv results.root
```

**Input format:**
```
index   X   Y   Z
```
Each line contains a point index and its X, Y, Z coordinates in millimeters.

---

## 🛠️ Building the Program

You need the **ROOT framework** (https://root.cern/) installed and configured.  
Then compile with:

```bash
clang++ -std=c++17 flatnessScan.cpp Points.cpp `root-config --cflags --libs` -o flatnessScan
```

or simply run:

```bash
make
```

(Provided that the included `Makefile` is configured for your ROOT installation.)

---

## 📦 Repository Contents

| File | Description |
|------|--------------|
| `flatnessScan.cpp` | Main program |
| `Points.cpp` / `Points.h` | Point reading utilities |
| `GridFinder.h` | Grid detection logic |
| `Makefile` | Build configuration |
| `README.md` | This documentation |
| Example `.csv` files | Sample datasets |

---

## 📊 Output

- ROOT file containing:
  - X, Y, Z coordinate histograms  
  - Plane deviation histogram  
  - Y–vs–X scatter plot  
  - Optional 2D “Flatness Map”  
- Console summary of fit parameters and flatness statistics  
- Interactive canvases (close windows or press Ctrl-C to exit)

---

## 🏷️ Version History

| Version | Date | Notes |
|----------|------|-------|
| **1.1.0** | October 2025 | First stable release with GridFinder integration and improved histogram handling. |

---

## 📜 License

© 2025 Luciano Ristori  
This project is provided for research and educational use.  
You may modify and redistribute it freely, provided attribution is maintained.

---

*For more details, see the [Releases](https://github.com/LucianoRistori/FlatnessScan/releases) section on GitHub.*

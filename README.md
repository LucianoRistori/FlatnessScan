# FlatnessScan

**FlatnessScan** is a C++ program for analyzing geometric or detector-related data to study surface or alignment flatness.  
It reads point data from input files, performs statistical analysis, and can produce ROOT histograms or summary outputs.

---

## 🧩 Features

- Reads input files containing 2D or 3D point data.  
- Performs flatness or deviation scans.  
- Optionally generates ROOT histograms for visualization.  
- Includes built-in protection against malformed input.  
- Modular structure (e.g., `Point.cpp`, `flatnessScan.cpp`).  

---

## ⚙️ Build Instructions

```bash
make flatnessScan
```
`

## 🧠 Usage Example

```bash
./flatnessScan input_points.txt
```

This reads the file `input_points.txt`, performs a flatness scan, and prints a short summary.  
If ROOT is enabled, histograms may be written to `flatness.root`.

---

## 🧰 Dependencies

- **ROOT** (https://root.cern)  
- **C++17 or later**  
- Tested on macOS 13+  

---

## 🧑‍💻 Author

**Luciano Ristori**  
GitHub: [LucianoRistori](https://github.com/LucianoRistori)

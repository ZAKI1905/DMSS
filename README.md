# DMSS (Dark Matter Solar Signals)
This program was used in "[Constraining Time Dependent Dark Matter Signals from the Sun](https://arxiv.org/abs/2109.11662)".

## Info 


If you use this package please cite:
```
  Constraining Time Dependent Dark Matter Signals from the Sun
  Mohammadreza Zakeri and Yu-Feng Zhou
  https://arxiv.org/abs/2109.11662
```

For BibTeX style citation use:
```
@article{Zakeri:2021cur,
    author = "Zakeri, Mohammadreza and Zhou, Yu-Feng",
    title = "{Constraining Time Dependent Dark Matter Signals from the Sun}",
    eprint = "2109.11662",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    month = "9",
    year = "2021"
}
```

## Setup
You need [Root](https://root.cern/), Minuit2 (included in Root), [GSL](https://www.gnu.org/software/gsl), and [OpenMP](https://www.openmp.org). [Confind](https://github.com/ZAKI1905/CONFIND) & Zaki libraries are already included inside dependencies (Only for Mac OS). Confind generates the contour plots and it uses OpenMP, but the main code for DMSS doesn't directly use OpenMP. The CMake file will use find_package() to find these packages. If find_package() fails, open CMakeLists.txt in the main folder and edit the cmake path for root. Compile as follows:
  1. Create a build folder (e.g. mkdir build).
  2. Change directory to the build folder (e.g. cd build)
  3. Run cmake (e.g. cmake ..)
  4. Run make (or make -j) 
  5. The executables should be created inside the build directory.


## References
Check CMakeLists.txt.

GSL:
```
M. Galassi et al, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078.
http://www.gnu.org/software/gsl/
```

Root:

See https://root.cern/about/license for Root's license.

```
Rene Brun and Fons Rademakers,
ROOT - An Object Oriented Data Analysis Framework,
Proceedings AIHENP'96 Workshop, Lausanne, Sep. 1996, Nucl. Inst. & Meth. in Phys. Res. A 389 (1997) 81-86. See also http://root.cern.ch/
```

Minuit2:
```
F. James, M. Roos,
Minuit - a system for function minimization and analysis of the parameter errors and correlations,
Computer Physics Communications,
Volume 10, Issue 6,
1975,
Pages 343-367,
ISSN 0010-4655,
https://doi.org/10.1016/0010-4655(75)90039-9.
(https://www.sciencedirect.com/science/article/pii/0010465575900399)
```

OpenMP:
```
The OpenMP API specification for parallel programming
https://www.openmp.org/
```

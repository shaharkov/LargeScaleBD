Large-Scale Bounded Distortion Mappings
====

A Matlab implementation of the paper ["Large-Scale Bounded Distortion Mappings"](https://shaharkov.github.io/projects/LargeScaleBD_lowRes.pdf).

----
The class `SolverProjectorBD.m` implements the algorithm described in the paper for computing bounded-distortion mappings via a Gauss-Newton-like approach.

Three examples are included in this package:
- `example_ProjBD_Mesh2D.m` demonstrates the computation a bounded-distortion mapping of a 2D triangular mesh.
- `example_ProjBD_Mesh3D.m` demonstrates the computation a bounded-distortion mapping of a 3D tetrahedral mesh.
- `example_ProjBD_Surface.m` demonstrates the computation a bounded-distortion mapping of a 3D surface triangular mesh.
The input to all examples are randomly generated mappings that have high distortions and flipped elements.

**Compatibility:** The code was tested with Matlab (2014b). The code depends on two `mex` files, Windows (x64) binaries are included; they are compiled with Intel C++ Composer XE 2016 with Microsoft Visual Studio 2013; compilation requires [Eigen](http://eigen.tuxfamily.org/). The source code is provided under the `mex/` folder; run `compileAllMex.m` to compile all mex files (only tested under windows).

**Disclaimer:**
The code is provided as-is for academic use only and without any guarantees. Please contact the author to report any bugs. 
Written by [Shahar Kovalsky](http://www.wisdom.weizmann.ac.il/~shaharko/).

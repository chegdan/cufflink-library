# News #
  * Cufflink has been tested and is known to work with the newest versions of Thrust and Cusp, as well as CUDA 4.2 (May 26<sup>th</sup>, 2012)
  * Be sure to check the [issue tracker](http://code.google.com/p/cufflink-library/issues/list) to see open issues and new fixes
  * Finally, cufflink is now released and ready for [download](http://code.google.com/p/cufflink-library/source/checkout) via [svn](http://svnbook.red-bean.com/en/1.4/index.html).  Thanks to all that helped get it this far (October 17<sup>th</sup>, 2011)
  * The [How to Contribute](HowToContibute.md) page was uploaded
  * The [Getting Started guide](GettingStartedPage.md) was uploaded
  * The [citations](citing.md) page was uploaded
  * See a presentation about [cufflink](http://dl.dropbox.com/u/35794017/daniel_combest_slides.pdf) given during the 6<sup>th</sup> OpenFOAM Workshop in June 2011.


---


# What is cufflink? #
Cuda For FOAM Link (cufflink) is an opensource library for linking numerical methods based on Nvidia's Compute Unified Device Architecture [(CUDA™)](http://developer.nvidia.com/cuda-toolkit-40) C/C++ programming language and OpenFOAM®.  Currently, the library utilizes the sparse linear solvers of [Cusp](http://code.google.com/p/cusp-library/) and methods from [Thrust](http://code.google.com/p/thrust/) to solve the linear Ax = b system derived from OpenFOAM's lduMatrix class and return the solution vector.  Cufflink is designed to utilize the course-grained parallelism of [OpenFOAM®](http://openfoam.org/) (via domain decomposition) to allow multi-GPU parallelism at the level of the linear system solver.  [Get Started](GettingStartedPage.md) with Cufflink.

# Current Features #
  * Currently only supports the [OpenFOAM-extend](http://www.extend-project.de) fork of the OpenFOAM code
  * A conjugate gradient solver based on [CUSP](http://code.google.com/p/cusp-library/) for symmetric matrices (e.g. pressure), with a choice of
    1. Diagonal preconditioner
    1. Sparse Approximate Inverse preconditioner
    1. Algebraic Multigrid (AMG) based on Smoothed Aggregation Precondtioner

  * A bi-conjugate gradient stabilized solver based on solver based on [CUSP](http://code.google.com/p/cusp-library/) for asymmetric matrices (e.g. velocity, epsilon, k), with a choice of
    1. Diagonal preconditioner
    1. Sparse Approximate Inverse preconditioner

  * Single GPU support (see [here](http://www.openfoamworkshop.org/6th_OpenFOAM_Workshop_2011/Program/Presentations/daniel_combest_slides.pdf) for some preliminary benchmarks and results)
  * Multi-GPU support (multi-gpu needs some more testing) via OpenFOAM's course grained parallelism achieved through [domain decomposition](http://www.openfoam.com/docs/user/running-applications-parallel.php#x12-820003.4).
  * Single Precision (sm\_10), Double precision (sm\_13), and Fermi Architecture (sm\_20) supported.  The double precision solvers are recommended over single precision due to known errors encountered in the Smoothed Aggregation Preconditioner in Single precision.

# System Requirements #
> Due to subtle differences in the forks of [OpenFOAM®](http://openfoam.org/), this library is only compatible with [OpenFOAM-extend](http://www.extend-project.de).  cufflink was developed and tested on [Linux Ubuntu 10.04 amd64 LTS](http://www.ubuntu.com/) with [CUDA 4.2](http://developer.nvidia.com/cuda-downloads), Thrust v1.5 (included in CUDA 4.2), and [Cusp 0.3.0](http://code.google.com/p/cusp-library/). To check your version of Thrust, Cusp, and cufflink, please compile the version.cu file in the cufflink install directory with
```
     nvcc version.cu -o version
```

and run in the `./version` command in the terminal to give the output

```
        Thrust v1.5
        Cusp   v0.3.0
        Cufflink   v0.1.0
```

> Check your version of nvcc with the command nvcc --version and the output should be similar to:

```
        nvcc: NVIDIA (R) Cuda compiler driver
        Copyright (c) 2005-2012 NVIDIA Corporation
        Built on Thu_Apr__5_00:24:31_PDT_2012
        Cuda compilation tools, release 4.2, V0.2.1221
```

> Also, ensure that your open mpi installation is installed correctly so that the runParallel test script will work and the multi-GPU solvers can be run.  To see whether you OpenFOAM install will work in parallel, test the Allrun script located in

```
   tutorials/multiphase/interFoam/laminar
```

folder by executing

```
   ./Allrun
```

Without openMPI working correctly, the multi-gpu solvers will not work.

# Suggested Hardware #
The amount of speedup is highly dependent on the GPU being used and the CPU the results are compared with.  As a basic suggestion, one should select at least a second generation Nvidia Tesla ([fermi](http://www.nvidia.com/object/fermi_architecture.html)) GPU with the appropriate power source and motherboard for optimal results.  Some information can be found [here](http://www.nvidia.com/object/tesla_computing_solutions.html).

# Disclaimer #
**Please note that cufflink is not approved or endorsed by Silicon Graphics International Corp., the owner of the OpenFOAM® trademark and producer of OpenFOAM® software.**

More information can be found at openfoam.org 
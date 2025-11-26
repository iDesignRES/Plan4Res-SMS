# SMSppThroughJulia.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com//SMSppThroughJulia.jl.svg?branch=master)](https://travis-ci.com//SMSppThroughJulia.jl)
[![codecov.io](http://codecov.io/github//SMSppThroughJulia.jl/coverage.svg?branch=master)](http://codecov.io/github//SMSppThroughJulia.jl?branch=master)

# How to get started

Prior to getting started one needs to install SMS++ on which this package is dependent. We will not require all subpackages of SMS++, but only:
- BundleSolver
- LagrangianDualSolver
- MILPSolver
- SMS++ (the core)
- Tools
- UCBlock

Here one can rely on the readily available installation instructions for SMS++ related to your platform:

https://gitlab.com/smspp/smspp-project#requirements

We will rely on internal use of the HiGHs solver, so be sure to at least add support for this solver. Support for the others can be deactivated, i.e., you can specify the additional options 

--without-cplex
--without-gurobi
--without-scip
--without-stopt

As these will not be needed.

<details>

<summary>Manual, i.e., - Die Hard - installation on Windows machines using MSBuild and VCPKG </summary>

### Manual installation on Windows machines

Download https://visualstudio.microsoft.com/fr/downloads/ but install only the command line tools

Check on/off the last version of the MS build tools for the command line ; Install ; You will need admin rights + internet access

Make some directory and

- git clone https://github.com/Microsoft/vcpkg.git

Open the Developer command prompt for VS

- In the just cloned vcpkg directory run : bootstrap-vcpkg.bat

- Next run : ./vcpkg.exe integrate install

- The subsequent commands have to be executed in the directory where the vcpkg has been (hand-)installed.

- Type : vcpkg.exe install zlib bzip2 blas lapack eigen3 glpk netcdf-cxx4 pthreads getopt boost --triplet x64-windows

- Wait some hours ( ~ 4h or 5h ) ...

- Go to the vcpkg directory then ports/coin-or-osi and edit portfile.cmake so as to allow for cplex:
        --with-cplex
        --with-cplex-lib=C:\\/IBM\\/ILOG\\/CPLEX_Studio221\\/cplex\\/lib\\/x64_windows_msvc14\\/stat_mda\\/cplex2210.lib
        --with-cplex-incdir=C:\\/IBM\\/ILOG\\/CPLEX_Studio221\\/cplex\\/include\\/ilcplex
		--with-cplex-cflags=-IC:\\/IBM\\/ILOG\\/CPLEX_Studio221\\/cplex\\/include\\/ilcplex
		--with-cplex-lflags=C:\\/IBM\\/ILOG\\/CPLEX_Studio221\\/cplex\\/lib\\/x64_windows_msvc14\\/stat_mda\\/cplex2210.lib
  (~ line 24)
  Alternatively type vcpkg.exe edit coin-or-osi (and the proper file should open by itself in your favourite editor)
  /!\ : it may be needed to not have any white spaces in the paths for this to work properly ; Also the \\/ for windows path specification are not optional !

Support for Gurobi can be enabled likewise by adding for instance

–with-gurobi
–with-gurobi-lib=C:\/gurobi1101\/win64/lib\/gurobi110.lib
–with-gurobi-incdir=C:\/gurobi1101\/win64\/include
–with-gurobi-cflags=-IC:\/gurobi1101\/win64\/include
–with-gurobi-lflags=C:\/gurobi1101\/win64\/lib\/gurobi110.lib

- Type vcpkg.exe install coin-or-osi --triplet x64-windows 

- Wait some 30 mins or so ; Check if your vcpkg/installed/x64-windows/lib directory has OsiCpx.lib

- Tyope vcpkg.exe install coin-or-clp coinutils pybind11 --triplet x64-windows

- Wait some more (1h30)

- Type vcpkg.exe install boost-mpi --triplet x64-windows

- In case of failure execute vcpkg\downloads\msmpisetup-10.1.12498.exe (or in fact the executable vcpkg tells you to install)

- Now try again : vcpkg.exe install boost-mpi --triplet x64-windows

- git clone https://gitlab.com/stochastic-control/vcpkg-registry

- Type vcpkg.exe install stopt --overlay-ports="YOURPATH"\vcpkg\vcpkg-registry\ports\stopt --triplet x64-windows

- Now we can get SMSpp ; Make an SMSpp directory and in it: 

	git clone --recurse-submodules https://gitlab.com/smspp/smspp-project -b develop

- Go to the SMSpp directory, make a buildVS directory

- Edit the extlib/ makefile-default paths to set the paths all properly

- Go in the buildVS directory ;

- cmake ../. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_TOOLCHAIN_FILE="YOURPATH"/vcpkg/scripts/buildsystems/vcpkg.cmake -Wno-dev -DMILPSolver_USE_GUROBI=OFF

Alternatively you can use the graphical interface of Cmake to toggle off any unwanted parts of SMS++ ;

- msbuild \"The SMS++ Project.sln"

> Yippee Ki‐Yay you are good to go

</details>

# How to use
Set up a directory cpp - where the dynamically compiled library emx_smspp_lib must go and some dependencies must go

Then follow up with 

using SMSppThroughJulia

SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_ramp.nc4")

Expected output:

[855.0, 829.9999999999993, 804.9999999999987, 779.9999999999983, 754.9999999999978, 729.9999999999974, 704.9999999999969, 679.9999999999965, 654.999999999996, 629.9999999999956, 604.9999999999951, 579.9999999999947, 554.9999999999942, 529.9999999999937, 504.9999999999933, 479.9999999999929, 454.9999999999925, 429.9999999999921, 404.9999999999917, 379.99999999999136, 354.999999999991, 329.9999999999907, 304.99999999999034, 279.99999999999]

SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_ramp2.nc4")
 
(In fact the same data but only the nature of the unit has changed which should change nothing to the output)

SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_mindown.nc4")

[854.9999999999994, 829.9999999999994, 804.9999999999994, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 
SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
, "test/data/NBlock_mindown2.nc4")

[855.0, 829.9999999999994, 804.999999999999, 800.0, 0.0, 0.0, 0.0, 0.0, 800.0, 830.0, 860.0, 890.0, 920.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0]
 
 SMSppThroughJulia.value_nuclear_on_price( [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0]
, "test/data/Nblock_mindown.nc4")

[855.0, 829.9999999999994, 804.999999999999, 800.0, 800.0, 800.0, 830.0, 860.0, 890.0, 920.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0, 925.0]
 
 SMSppThroughJulia.value_nuclear_on_price( [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_mod.nc4")

[887.0, 886.0, 885.0, 884.0, 883.0, 882.0, 881.0, 856.0, 855.0, 854.0, 853.0, 828.0, 827.0, 826.0, 825.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

 SMSppThroughJulia.value_nuclear_on_price( [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/NBlock_mod2.nc4")

[905.0, 910.0, 909.0, 884.0, 883.0, 882.0, 881.0, 856.0, 855.0, 854.0, 853.0, 828.0, 827.0, 826.0, 825.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 
 SMSppThroughJulia.value_nuclear_on_price( [-10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, -10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
, "test/data/Nblock_mod3.nc4")
 
[910.0000000000005, 911.0000000000005, 912.000000000001, 913.000000000001, 925.0, 925.0, 920.000000000033, 910.0, 900.0, 890.0, 880.0, 855.0, 845.0, 835.0, 825.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Alternative - automated test

Alternatively one can run

using SMSppThroughJulia

followed up by 

SMSppThroughJulia.test()

then one should see displayed

Testing : 8 / 8 - passed => indicating succes



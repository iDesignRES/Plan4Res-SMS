# Plan4Res - SMS++

SMS++ is a set of C++ classes intended to provide a system for modeling complex, block-structured mathematical models (in particular, but not exclusively, single-real-objective optimization problems), and solving them via sophisticated, structure-exploiting algorithms (in particular, but not exclusively, decomposition approaches and structured Interior-Point methods). SMS++ also supports and interfaces with classic solvers such as CPLEX, GUROBI, HiGHS, a bundle method and StOpt - an open source SDDP solver. The SMS++ system can be used to put together any block-structured mathematical programs and exploit the thus available structure.

Plan4Res is an instantation of the SMS++ framework for energy modelling. It comprises a three-layer set of "tools" that are embedded. These are:
  - Investment Layer
  - Seasonal Storage Valuation Layer
  - Unit-Commitment Layer (= Economic dispatch if no dynamic constraints are given)

Each layer relies on lower layers and can thus be consistently set up with more detailled smaller computations. Support for various solvers, whenever possible / implemented, can be harnessed through simple configuration files. As a result, the same model can be solved using various solution approaches and results compared. The Plan4Res toolset takes as input "netcdf" files that describe the system at hand. The system itself is agnostic to units and geography and as a result computations can be made in various units (for instance hydro storage in volume m3 or energy MWh) and at various geographical sizes. The Plan4Res system has been used in various EU sponsored projets : 
  - Plan4Res
  - OpenEntrance (for instance in the case study 1 : https://doi.org/10.5281/zenodo.7871106 , https://doi.org/10.5281/zenodo.7997102 , https://doi.org/10.5281/zenodo.8086793 -  but also 4,5 7: https://doi.org/10.5281/zenodo.8289102 )
  - OpenMod4Africa
  - CETPartnership "RESILIENT - Resilient Energy System Infrastructure Layouts for Industry, E-Fuels and Network Transition"
  - CETPartnership "Manoeuvre"

Moreover the mathematics and solution methods are state-of-the-art in mathematical programming.

For the purpose of iDesignRES two subfeatures are of particular interest:
  - The possibility to specify the network as Flow based Network, a DCNetwork or an AC Network ; The specific mode is triggered by simply changing the netcdf file given as input to the unit-commitment solver
  - The addition of a NuclearUnitBlock : a specific version of ThermalUnitBlock with specific constraints regarding the modulation of power output. Indeed Nuclear units can provide significant flexibility to the system, but with additional specific constraints.

# OS support

The framework is available on any specific platform and can be installed / compiled in two specific ways:
  - By pulling through the use of git all relevant source code from the link given below. Then by compiling using CMAKE, possibly in combination with VCPKG or by using classic makefiles ;
  - By using the singularity container.
The prefered mode should be the former, as the latter may get deprecated. 

# Where to find SMS++

[SMS++ Main Repository](https://gitlab.com/smspp/smspp-project)

[Plan4Res through Singularity](https://gitlab.com/cerl/plan4res/p4r-env)




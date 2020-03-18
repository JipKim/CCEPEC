# CCEPEC

This code was developed under Julia v1.2.0/JuMP0.20 by Jip Kim in 2019.
The following packages must be installed:

  - Gurobi
  - JuMP
  - Distributions
  - LinearAlgebra
 
To run the code, execute PH.jl, or include() it from a Julia prompt.

"ISONE8busTN" testsystem is given as follows in data folder:
  - Node.csv: Transmission network data
  - Line.csv: Transmission line data
  - Generator.csv: Generation data
  - We also specify additional input data in the beginning of PHiniKKT.jl/PHiterKKT.jl

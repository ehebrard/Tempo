# Tempo

## Build
Requires: cmake

1/ Create a build folder 
> mkdir build
> cd build

2/ Configure
> cmake -DCMAKE_CXX_COMPILER=g++-11 -DCMAKE_BUILD_TYPE=release ..

3/ Compile the executables
> make 

## Examples
There are a few examples of constraint programs in src/examples
1/ tempo_scheduler.cpp: solves disjunctive resource scheduling problems in various format 
2/ bibd: solves the "balanced incomplete block design" problem (pb 28 in the CSPLIB https://www.csplib.org/Problems/prob028/)
3/ tempo_sat.cpp: solves CNF formulas in dimacs format

To compile example xxx:
> make xxx 

## Running
Executables usually require a positional argument standing for the instance file, all other arguments have a flag. Run <exec> --help



## DEBUG CLAUSE LEARNING
1/ define the macro DBG_CL X (X is a number standing for the maximum number of generated clauses to debug)
2/ run tempo_scheduler with the option --dbg <filename> (it outputs clauses in <filename>)
3/ run the same cmd line as above, but with 'clausecheck' instead of tempo_scheduler


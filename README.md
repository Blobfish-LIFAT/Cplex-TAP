# Cplex-TAP
Finding exact solutions for the Travelling Analyst Problem.


## Build
This has only been tested on linux (ubuntu 20.04 and fedora 33) , you'll need c++ build tools for your distribution namely g++ and cmake

First clone the repo and swap the main for the production file
```shell
git clone https://github.com/Blobfish-LIFAT/Cplex-TAP
cd Cplex-TAP
mv exemple_production.cpp exemple_production.bak
cp prod_ready.bak exemple_production.cpp
```

Note : By default the solver times out after 1h this can be changed in the solver.ccp file line 67.


## Acknwoledgements
The Combo Knapsack algorithm and code (combo.h, combo.cpp) is the work of S.Martello, D.Pisinger and P.Toth. Used for academic purpose in accordence with it's licence. 


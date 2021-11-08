# Cplex-TAP
Finding exact solutions for the Travelling Analyst Problem.


## Build
### Requirements
This has only been tested on linux (ubuntu 20.04 and fedora 33) , you'll need c++ build tools for your distribution namely g++ and cmake

You'll have to install CPLEX on your system, the build scripts will only look for it in the default location if you installed it elsewhere for some reason you can edit the script FindCplex or create some symbolic links.
### Instructions
First clone the repo and swap the main for the production file
```shell
git clone https://github.com/Blobfish-LIFAT/Cplex-TAP
cd Cplex-TAP
mv exemple_production.cpp exemple_production.bak
cp prod_ready.bak exemple_production.cpp
```

Note : By default the solver times out after 1h this can be changed in the solver.ccp file line 67.

In order to build the binary 
```shell
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_DEPENDS_USE_COMPILER=FALSE -G "CodeBlocks - Unix Makefiles" ./
cmake --build ./ --target tap -- -j 6
```
The artifact should end up in ~/build/bin/tap

If the linker doesn't behave you'll have to add this line `set(CMAKE_CXX_LINK_EXECUTABLE  "${CMAKE_CXX_LINK_EXECUTABLE} -ldl")
` at the end of the file CmakeLists.txt and rerun the two previous commands.
## Acknwoledgements
The Combo Knapsack algorithm and code (combo.h, combo.cpp) is the work of S.Martello, D.Pisinger and P.Toth. Used for academic purpose in accordence with it's licence. 


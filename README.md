# 1. INSTALLING aScan

aScan depends on the following libraries that should be installed in your system in order to be able to compile aScan:

1) GSL - GNU Scientific Library https://www.gnu.org/software/gsl/
2) zlib - https://zlib.net/
3) bamtools https://github.com/pezmaster31/bamtools

For convenience, the bamtools library source is included in the aScan package, build it as follows.

# 1.1 COMPILING bamtools

The cmake utility (https://cmake.org/) is required to compile the bamtools library. 

From the aScan folder type:

cd bamtools

mkdir build

cd build

cmake ..

make

For further support on installing the bamtools refer to https://github.com/pezmaster31/bamtools/wiki/Building-and-installing

# 1.2 COMPILING aScan

To compile aScan run the following command from within the aScan folder

g++ aScan_dev.cpp -o aScan -I bamtools/src/ -lbamtools -L bamtools/build/src/api -lgsl -lgslcblas -lz -std=c++11 -Wall -pthread

This command assumes you compiled the bamtools library as described in 1.1. If you prefer to point to another bamtools installation, just modify the command to match your bamtools path.

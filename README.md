# aScan
g++ aScan_dev.cpp -o aScan -I bamtools/src/ -lbamtools -L bamtools/build/src/api -lgsl -lgslcblas -lz -std=c++11 -Wall -pthread

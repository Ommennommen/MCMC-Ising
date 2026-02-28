#pragma once
#include <memory>

struct Lattice {
    int N;
    std::unique_ptr<int[]> spins;

    Lattice(int n);

    int& at(int i, int j);
    const int& at(int i, int j) const;
};
void initializeLattice(Lattice& lattice);
void mcmc(Lattice& config, double beta);
int calcEnergy(Lattice& lattice);
int calcMag(Lattice& lattice);
void saveCSV(const std::string& filename, double* arr, int size);
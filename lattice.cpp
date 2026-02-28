#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <memory>
#include "lattice.h"


Lattice::Lattice(int n) : N(n) {
    spins = std::make_unique<int[]>(N * N);
}


int& Lattice::at(int i, int j) {
    return spins[i * N + j];
}

const int& Lattice::at(int i, int j) const {
    return spins[i * N + j];
}

void initializeLattice(Lattice& lattice) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> chance(0, 1);

    for (int i = 0; i < lattice.N; ++i) {
        for (int j = 0; j < lattice.N; ++j) {
            lattice.at(i, j) = chance(gen) == 0 ? -1 : 1;
        }
    }
}



void mcmc(Lattice& config, double beta){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> chance(0, config.N - 1);
    std::uniform_real_distribution<> rand(0, 1);

    int N = config.N;

    for (int i = 0; i < N; ++i){
        for (int j=0; j < N; ++j){
            int a = chance(gen);
            int b = chance(gen);
            int s = config.at(a, b);
            int nb = config.at((a+1)%N, b) + config.at(a, (b+1)%N) + config.at((a-1)%N, b) + config.at(a, (b-1)%N);

            int cost = 2*s*nb;

            if (cost < 0 || rand(gen) < exp(-cost*beta)){
                s *= -1;
                config.at(a, b) = s;
            }

        }
    }

}


int calcEnergy(Lattice& lattice) {
    int energy = 0;
    int N = lattice.N;
    for (int i =0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            int S = lattice.at(i, j);
            int nb = lattice.at((i+1)%N,j) + lattice.at(i,(j+1)%N) + lattice.at((i-1)%N, j) + lattice.at(i,(j-1)%N);
            energy += -nb*S;
        }
    }
    return energy/2;
}

int calcMag(Lattice& lattice) {
    int N = lattice.N;
    int mag = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j){
            mag += lattice.at(i, j);
        }
    }
    return mag;
}

void saveCSV(const std::string& filename, double* arr, int size) {
    std::ofstream file(filename);
    for (int i =0; i < size; ++i){
        if (i > 0) file << ",";
        file << arr[i];
    }
    file.close();
}


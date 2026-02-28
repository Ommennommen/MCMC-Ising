#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include "lattice.h"


int main() {

    // Parameters
    int N = 16; // Lattice size
    const int nT = 21;
    double step = (4.0 - 1.0) / (nT - 1);
    double temp[nT];
    for (int i = 0; i < nT; ++i){
        temp[i] = 1.0 + i * step;
    }
    int NC = 400;  //Number of configurations per temperature point
    int eqSteps = 1 << 8;  //Number of MC sweeps for equilibration
    int mcSteps = 1 << 9;  //Number of MC sweeps for calculation

    double E[nT];
    double M[nT];
    double C[nT];
    double X[nT];
    for (int i=0; i<nT; ++i){
        E[i] = 0;
        M[i] = 0;
        C[i] = 0;
        X[i] = 0;
    }
    double n1 = 1.0/(mcSteps*N*N); 
    double n2 = 1.0/(mcSteps*mcSteps*N*N);

    std::vector<int> data(nT*NC*N*N);

    for (int i = 0; i < nT; ++i){
        double iT = 1/temp[i];
        for (int j = 0; j < NC; ++j){
            Lattice lattice(N);
            initializeLattice(lattice);
            for (int k = 0; k < eqSteps; ++k){
                mcmc(lattice, iT);
            };
        
            for (int l = 0; l < lattice.N; ++l) {
                for (int m = 0; m < lattice.N; ++m) {
                    data[(i*NC + j) * N*N + l*N + m] = lattice.at(l, m);
                };
            }
        }


        double E1 = 0;
        double M1 = 0;
        double E2 = 0;
        double M2 = 0;
        double iT2=iT*iT;
        Lattice lattice(N);
        initializeLattice(lattice);
        for (int tt = 0; tt < mcSteps; ++tt) {
            mcmc(lattice, iT);       
            int Ene = calcEnergy(lattice);   //calculate the energy
            int Mag = calcMag(lattice);      //calculate the magnetisation

            E1 = E1 + Ene;
            M1 = M1 + Mag;
            M2 = M2 + Mag*Mag;
            E2 = E2 + Ene*Ene;
        }
        //divide by number of sites and iteractions to obtain intensive values    
        E[i] = n1*E1;
        M[i] = n1*M1;
        C[i] = (n1*E2 - n2*E1*E1)*iT2;
        X[i] = (n1*M2 - n2*M1*M1)*iT;

        std::cout << "Temperature " << i+1 << "/" << nT << std::endl;

    }
saveCSV("ising_temps.csv", temp, nT);
saveCSV("ising_energy.csv", E, nT);
saveCSV("ising_mag.csv", M, nT);
saveCSV("ising_specHeat.csv", C, nT);
saveCSV("ising_susc.csv", X, nT);

}
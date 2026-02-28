#include <iostream>
#include <random>
#include <cmath>
#include <fstream>


struct Lattice {
    int N;
    int* spins;

    Lattice(int n);
    ~Lattice();

    int& at(int i, int j);
    const int& at(int i, int j) const;
};



Lattice::Lattice(int n) : N(n) {
    spins = new int[N * N];
}

Lattice::~Lattice() {
    delete[] spins;
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

    int* data = new int[nT*NC*N*N];

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
    std::ofstream file("ising_configs.csv");
    for (int row = 0; row < nT*NC; ++row){
        for (int col = 0; col < N*N; ++col){
            file << data[row * N*N +col] << ",";
        }
        file << "\n";
    }
    file.close();

    std::ofstream file2("ising_temps.csv");
    for (int T =0; T < nT; ++T){
        if (T > 0) file2 << ",";
        file2 << temp[T];
    }
    file2.close();

    std::ofstream file3("ising_energy.csv");
    for (int e =0; e < nT; ++e){
        if (e > 0) file3 << ",";
        file3 << E[e];
    }
    file3.close();

    std::ofstream file4("ising_mag.csv");
    for (int m =0; m < nT; ++m){
        if (m > 0) file4 << ",";
        file4 << M[m];
    }
    file4.close();

    std::ofstream file5("ising_specHeat.csv");
    for (int c =0; c < nT; ++c){
        if (c > 0) file5 << ",";
        file5 << C[c];
    }
    file5.close();

    std::ofstream file6("ising_susc.csv");
    for (int s =0; s < nT; ++s){
        if (s > 0) file6 << ",";
        file6 << X[s];
    }
    file6.close();

    delete[] data;
}
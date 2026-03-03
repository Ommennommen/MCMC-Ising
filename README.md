# 2D Ising Model — MCMC Simulation
A C++ implementation of the 2D Ising model using Markov Chain Monte Carlo (MCMC) with the Metropolis algorithm. Written as a C++ learning project, translated from a prior Python implementation used in a machine learning research project.

## Physics
The Ising model describes a 2D lattice of magnetic spins, each either up (+1) or down (-1). Spins prefer to align with their neighbours, but thermal energy can flip them. The competition between these two effects produces a phase transition at the critical temperature $T_c \approx 2.27$ (in units where $J = k_B = 1$):
- **Low T**: Spins align, forming ordered ferromagnetic domains
- **High T**: Spins are random, producing a disordered paramagnetic state
- **Near $T_c$**: Large fluctuations, diverging susceptibility and specific heat

The Hamiltonian is:

$$H = -J \sum_{\langle i,j \rangle} s_i \cdot s_j$$

summed over nearest-neighbour pairs, with $J = 1$.

## Metropolis Algorithm
At each Monte Carlo step:
1. Randomly select a spin
2. Compute the energy change $\Delta E = 2s \cdot \sum_{\text{neighbours}}$ if flipped
3. Accept the flip if $\Delta E < 0$, or with probability $e^{-\beta \Delta E}$ otherwise

This satisfies detailed balance and samples the Boltzmann distribution. Periodic boundary conditions (torus topology) eliminate edge effects.

## Thermodynamic Observables
After equilibration, the following intensive quantities are measured via Monte Carlo averaging:

**Energy per spin:**
$$\langle E \rangle = \frac{1}{N^2} \sum E$$

**Magnetisation per spin:**
$$\langle M \rangle = \frac{1}{N^2} \sum |M|$$

**Specific heat:**
$$C = \beta^2 \left( \langle E^2 \rangle - \langle E \rangle^2 \right)$$

**Magnetic susceptibility:**
$$\chi = \beta \left( \langle M^2 \rangle - \langle M \rangle^2 \right)$$

All quantities are per spin.

## Performance
This C++ implementation completes a full simulation run (21 temperature points, 2000 configurations per temperature, 16×16 lattice) in approximately 13 minutes. The equivalent Python implementation required ~8 hours for the same parameters — a ~37x speedup.

## File Structure
```
MCMC/
├── main.cpp        # Simulation loop, data collection, CSV output
├── lattice.cpp     # Lattice implementation and physics functions
└── lattice.h       # Struct definition and function declarations
```

## Building
```bash
g++ main.cpp lattice.cpp -o MCMC
./MCMC
```

## Parameters
All parameters are set in `main.cpp`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N` | 16 | Lattice size ($N \times N$ grid) |
| `nT` | 21 | Number of temperature points |
| `NC` | 2000 | Configurations generated per temperature |
| `eqSteps` | 256 | MC sweeps for equilibration |
| `mcSteps` | 512 | MC sweeps for observable measurement |

Temperature range is $1.0$ to $4.0$, spanning the critical temperature $T_c \approx 2.27$.

## Output
The simulation produces six CSV files:

| File | Contents |
|------|----------|
| `ising_configs.csv` | All generated spin configurations ($n_T \times N_C$ rows, $N^2$ columns) |
| `ising_temps.csv` | Temperature values |
| `ising_energy.csv` | Average energy per spin vs temperature |
| `ising_mag.csv` | Average magnetisation per spin vs temperature |
| `ising_specHeat.csv` | Specific heat vs temperature |
| `ising_susc.csv` | Magnetic susceptibility vs temperature |

The configurations dataset can be used as input for machine learning models trained to predict temperature from spin configuration.

## Plotting
```bash
python plot.py
```
Requires `numpy` and `matplotlib`. Plots energy, magnetisation, specific heat, and susceptibility against temperature, showing the phase transition near $T_c \approx 2.27$.

## Background
This simulation was originally developed as part of an undergraduate machine learning research project investigating how phase transitions in the 2D Ising model are encoded in the weight matrices of convolutional neural networks. The original Python implementation has since been adopted for teaching purposes. This C++ rewrite was undertaken as a learning project to develop systems programming skills in preparation for robotics research.

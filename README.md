# 2D Ising Model — MCMC Simulation

A C++ implementation of the 2D Ising model using Markov Chain Monte Carlo (MCMC) with the Metropolis algorithm. Written as a C++ learning project, translated from a prior Python implementation used in a machine learning research project.

## Physics

The Ising model describes a 2D lattice of magnetic spins, each either up (+1) or down (-1). Spins prefer to align with their neighbours, but thermal energy can flip them. The competition between these two effects produces a phase transition at the critical temperature Tc ≈ 2.27 (in units where J = kB = 1):

- **Low T**: Spins align, forming ordered ferromagnetic domains
- **High T**: Spins are random, producing a disordered paramagnetic state
- **Near Tc**: Large fluctuations, diverging susceptibility and specific heat

The Hamiltonian is:
```
H = -J Σ s_i · s_j
```
summed over nearest-neighbour pairs, with J = 1.

## Metropolis Algorithm

At each Monte Carlo step:
1. Randomly select a spin
2. Compute the energy change ΔE = 2s·(sum of neighbours) if flipped
3. Accept the flip if ΔE < 0, or with probability exp(-βΔE) otherwise

This satisfies detailed balance and samples the Boltzmann distribution. Periodic boundary conditions (torus topology) eliminate edge effects.

## Thermodynamic Observables

After equilibration, the following intensive quantities are measured via Monte Carlo averaging:

| Observable | Formula |
|---|---|
| Energy | `<E> = n1 · ΣE` |
| Magnetisation | `<M> = n1 · Σ\|M\|` |
| Specific heat | `C = β² · (<E²> - <E>²)` |
| Susceptibility | `χ = β · (<M²> - <M>²)` |

All quantities are per spin.

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
| `N` | 16 | Lattice size (N×N grid) |
| `nT` | 21 | Number of temperature points |
| `NC` | 400 | Configurations generated per temperature |
| `eqSteps` | 256 | MC sweeps for equilibration |
| `mcSteps` | 512 | MC sweeps for observable measurement |

Temperature range is 1.0 to 4.0, spanning the critical temperature Tc ≈ 2.27.

## Output

The simulation produces six CSV files:

| File | Contents |
|------|----------|
| `ising_configs.csv` | All generated spin configurations (nT×NC rows, N² columns) |
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

Requires `numpy` and `matplotlib`. Plots energy, magnetisation, specific heat, and susceptibility against temperature, showing the phase transition near Tc ≈ 2.27.
# Kinetic Monte Carlo Simulation of hBN Growth

## Overview

This repository contains a lattice-based Kinetic Monte Carlo (KMC) simulator for the
epitaxial growth of monolayer hexagonal boron nitride (hBN) on a graphene substrate.

The simulator models the nonequilibrium growth process occurring during Molecular
Beam Epitaxy (MBE), where Boron (B) and Nitrogen (N) atoms are deposited onto a
graphene surface and subsequently undergo diffusion, desorption, and bond-preserving
"swing" movements.

Surface evolution is simulated using the rejection-free Gillespie Kinetic Monte Carlo
algorithm. Transition rates are computed from Arrhenius kinetics using activation
energy barriers obtained from Density Functional Theory (DFT) and Nudged Elastic
Band (NEB) calculations available at: 

Ren Zhong,
*Kinetic Monte Carlo Simulation of MBE Growth of Layered Hexagonal Boron Nitride*,
M.S. Thesis, Cornell University, 2019.

## Repository Structure

```
.
├── hBN.cpp
├── config.h
├── restart.h
├── README.md
├── LICENSE
├── CITATION.cff
└── images/

```
## Physical Model

The simulated system consists of

• Graphene substrate
• Boron adatoms
• Nitrogen adatoms
• Monolayer hexagonal Boron Nitride (hBN)

The substrate is represented as a periodic hexagonal lattice where each lattice site
corresponds to one possible adsorption position.

Each lattice site stores

0 → Empty site

1 → Boron atom

2 → Nitrogen atom

------------------------------------------------------------

## Lattice Geometry

The simulation domain is a periodic

    ROW × COL

hexagonal lattice.

Periodic boundary conditions are applied in both directions to eliminate edge effects.

Two neighbor definitions are used throughout the simulation:

Nearest Neighbors (NN)

These are the three neighboring sites capable of forming B-N chemical bonds.

Next Nearest Neighbors (NNN)

These are the six possible diffusion destinations for an atom.

Neighbor lookup tables are defined in

config.h

using

NN_x
NN_y
NNN_x
NNN_y

------------------------------------------------------------

## Physical Processes

Four microscopic processes are included in the simulation.

1. Deposition

Boron and Nitrogen atoms arrive alternately at a constant deposition rate.

Each arriving atom is placed on a randomly selected allowed lattice site.

------------------------------------------------------------

2. Surface Diffusion

Atoms thermally hop to one of their empty next-nearest-neighbor sites.

The hopping rate follows Arrhenius kinetics

r = ν exp(-Eb / kBT)

where

ν   = Attempt frequency

Eb  = Activation energy barrier

kB  = Boltzmann constant

T   = Substrate temperature

------------------------------------------------------------

3. Desorption

An isolated atom (zero B-N bonds) may leave the surface.

This process competes with surface diffusion.

------------------------------------------------------------

4. Swing Motion

When an atom is connected by exactly one bond, it may rotate around the bonded
neighbor without breaking the bond.

Swing events are assigned separate activation barriers.

------------------------------------------------------------

## Energy Barrier Model

Diffusion barriers depend on

• Atom type (B or N)

• Number of bonds before diffusion

• Number of bonds after diffusion

Separate barrier matrices are stored for

Energy_Barriers_Boron

Energy_Barriers_Nitrogen

These values originate from first-principles DFT/NEB calculations and determine
the Arrhenius transition rates used throughout the simulation.

------------------------------------------------------------

## Event Generation

Rather than rebuilding every possible event after each KMC step, the simulator
generates only events in the local neighborhood surrounding modified lattice sites.

The update procedure is

1. Remove outdated events around affected sites.

2. Scan nearby atoms.

3. Recalculate all valid diffusion, desorption and swing events.

4. Insert the new events into the global event list.

This localized update strategy significantly improves computational efficiency while
preserving the exact KMC dynamics.

------------------------------------------------------------

## Kinetic Monte Carlo Algorithm

Each simulation step performs

1. Generate all possible events.

2. Compute each event rate.

3. Calculate

       R = Σ ri

4. Draw a random number.

5. Select one event with probability

       Pi = ri / R

6. Advance physical time

       Δt = -ln(U)/R

7. Execute the selected event.

8. Update only the affected local region.

9. Repeat.

This is the standard rejection-free Gillespie algorithm.

------------------------------------------------------------

## Simulation Parameters

Simulation constants are defined in

config.h

Important parameters include

Temperature

Attempt frequency

Boltzmann constant

Lattice dimensions

Deposition interval

Maximum number of KMC iterations

Output frequency

Neighbor lookup tables

Energy barrier matrices

Changing these parameters allows different growth conditions to be studied without
modifying the simulation code.

------------------------------------------------------------

## Program Structure

config.h

Contains

• Physical constants

• Simulation parameters

• Neighbor lookup tables

• Energy barrier matrices

------------------------------------------------------------

restart.h

Implements

• Restart file loading

• Restart file saving

• Simulation continuation

------------------------------------------------------------

hBN.cpp

Contains

• Event generation

• Gillespie event selection

• Arrhenius rate calculation

• Deposition

• Diffusion

• Desorption

• Swing events

• Local event updates

• Simulation loop

------------------------------------------------------------

## Output Files

### outhex.txt

Contains atomic coordinates in OVITO format.

This file can be directly visualized using OVITO to observe surface morphology
evolution.

------------------------------------------------------------

### classhex.csv

Stores the complete simulation state including

• Lattice configuration

• Event counters

• Physical simulation time

• Restart information

This file enables interrupted simulations to continue from the previous state.

------------------------------------------------------------

## Terminal Output

During execution the program reports

Current iteration

Coverage

Number of deposited atoms

Number of desorbed atoms

Number of diffusion events

Number of swing events

Number of active events

Estimated physical simulation time

Wall-clock execution time

------------------------------------------------------------

## Compilation

Compile using

g++ -O3 -march=native hBN.cpp -o hBN

Run using

./hBN


------------------------------------------------------------


Acknowledgement

The physical model and energy barrier tables implemented in this code are based on the work of Ren Zhong, Kinetic Monte Carlo Simulation of MBE Growth of Layered Hexagonal Boron Nitride, M.S. Thesis, Cornell University, 2019. This repository provides an independent C++ implementation of the simulation algorithm with improvements in code organization, restart capability, and local event updates.

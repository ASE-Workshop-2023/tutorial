---
title: "Calculations in parallel"
teaching: 20
exercises: 20
questions:
objectives:
keypoints:
---

## Calculators using a parallel interpreter

GPAW is a bit special! It is an electronic structure code implemented as a Python library with C backend, and GPAW development is closely related to ASE development.

The difference between GPAW and other external Python calculators is that it includes a special Python interpreter for parallel calculations.

TODO:
- GPAW run in notebook using xc and kpts arguments
- time the run
- create a script file with equivalent settings and 4-way kpt parallelism
- compare runtime running on 4 cores

To begin with, we run a single-point energy calculation using Kohn-Sham density-functional theory (DFT).
- `xc` sets the exchange-correlation functional
- `kpts` sets the Brillouin-zone sampling
- `mode` sets the basis set; in this case a 400 eV cutoff plane-wave basis is used.

For more information about valid parameters, [see the GPAW docs](https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#parameters).

~~~
from gpaw import GPAW, PW

atoms = ase.build.bulk('Cu')
atoms.calc = GPAW(xc='PBE', kpts=(3, 3, 3), mode=PW(400))
energy = atoms.get_potential_energy()
print(energy)
~~~
{: .python}

~~~
from gpaw import GPAW, PW

atoms = ase.build.bulk('Cu')
atoms.calc = GPAW(xc='PBE', kpts=(3, 3, 3), mode=PW(400))
energy = atoms.get_potential_energy()
print(energy)
1
from gpaw import GPAW, PW
2
​
3
atoms = ase.build.bulk('Cu')
4
atoms.calc = GPAW(xc='PBE', kpts=(3, 3, 3), mode=PW(400))
5
energy = atoms.get_potential_energy()
6
print(energy)

  ___ ___ ___ _ _ _  
 |   |   |_  | | | | 
 | | | | | . | | | | 
 |__ |  _|___|_____|  22.8.0
 |___|_|             

User:   adam@Arctopus
Date:   Tue Apr  4 11:30:25 2023
Arch:   x86_64
Pid:    97252
CWD:    /home/adam/src/ase-tutorial-2023/development
Python: 3.10.0
gpaw:   /home/adam/.conda/envs/user-base/envs/ase-tutorials/lib/python3.10/site-packages/gpaw
_gpaw:  /home/adam/.conda/envs/user-base/envs/ase-tutorials/lib/python3.10/site-packages/
        _gpaw.cpython-310-x86_64-linux-gnu.so
ase:    /home/adam/src/ase/ase (version 3.23.0b1-70eab133b6)
numpy:  /home/adam/.conda/envs/user-base/envs/ase-tutorials/lib/python3.10/site-packages/numpy (version 1.23.5)
scipy:  /home/adam/.conda/envs/user-base/envs/ase-tutorials/lib/python3.10/site-packages/scipy (version 1.10.1)
libxc:  5.2.3
units:  Angstrom and eV
cores: 1
OpenMP: True
OMP_NUM_THREADS: 1

Input parameters:
  kpts: [3 3 3]
  mode: {ecut: 400.0,
         name: pw}
  xc: PBE

System changes: positions, numbers, cell, pbc, initial_charges, initial_magmoms 

Initialize ...
~~~
{: .output}

We should check the convergence of the energy with respect to k-point sampling. As well as accepting a simple mesh as e.g. `[3, 3, 3]`, we can pass a dictionary specification that generates a more specific mesh. In this case, we will shift the meshes to be off the Gamma-point.

We will also collect some timing information to understand how this impacts calculation cost. The output is directed to a text file to avoid cluttering our notebook.

This will take a few minutes to run: to see live output, open a terminal and use `tail -f kpts_serial.txt` to see this file grow. When it is finished, you can exit `tail` with ctrl-c.

~~~
from ase.utils.timing import Timer
timer = Timer()
energies, times, nkpts = [], [], []


for k in range(3,9):
    atoms.calc = GPAW(mode=PW(400), xc='PBE',
                      kpts={'size': [k, k, k],
                            'gamma': False},
                      txt='kpts_serial.txt')
    timer.start(str(k))
    energies.append(atoms.get_potential_energy())
    timer.stop(str(k))
    times.append(timer.get_time(str(k)))
    nkpts.append(len(atoms.calc.get_ibz_k_points()))
~~~
{: .python}

~~~
Timing:                              incl.     excl.
-----------------------------------------------------------
Hamiltonian:                         0.096     0.000   0.0% |
 Atomic:                             0.088     0.001   0.0% |
  XC Correction:                     0.087     0.087   0.8% |
 Calculate atomic Hamiltonians:      0.001     0.001   0.0% |
 Communicate:                        0.000     0.000   0.0% |
 Initialize Hamiltonian:             0.000     0.000   0.0% |
 Poisson:                            0.000     0.000   0.0% |
 XC 3D grid:                         0.007     0.007   0.1% |
LCAO initialization:                 0.344     0.082   0.7% |
 LCAO eigensolver:                   0.147     0.000   0.0% |
  Calculate projections:             0.000     0.000   0.0% |
  DenseAtomicCorrection:             0.000     0.000   0.0% |
  Distribute overlap matrix:         0.000     0.000   0.0% |
  Orbital Layouts:                   0.001     0.001   0.0% |
  Potential matrix:                  0.144     0.144   1.3% ||
  Sum over cells:                    0.001     0.001   0.0% |
 LCAO to grid:                       0.018     0.018   0.2% |
 Set positions (LCAO WFS):           0.097     0.012   0.1% |
  Basic WFS set positions:           0.002     0.002   0.0% |
  Basis functions set positions:     0.000     0.000   0.0% |
  P tci:                             0.014     0.014   0.1% |
  ST tci:                            0.039     0.039   0.3% |
  mktci:                             0.030     0.030   0.3% |
PWDescriptor:                        0.017     0.017   0.2% |
SCF-cycle:                           1.907     0.101   0.9% |
 Davidson:                           0.459     0.152   1.3% ||
  Apply H:                           0.033     0.027   0.2% |
   HMM T:                            0.006     0.006   0.1% |
  Subspace diag:                     0.076     0.004   0.0% |
   calc_h_matrix:                    0.056     0.023   0.2% |
    Apply H:                         0.034     0.027   0.2% |
     HMM T:                          0.006     0.006   0.1% |
   diagonalize:                      0.009     0.009   0.1% |
   rotate_psi:                       0.006     0.006   0.1% |
  calc. matrices:                    0.160     0.091   0.8% |
   Apply H:                          0.069     0.056   0.5% |
    HMM T:                           0.012     0.012   0.1% |
  diagonalize:                       0.023     0.023   0.2% |
  rotate_psi:                        0.014     0.014   0.1% |
 Density:                            0.191     0.000   0.0% |
  Atomic density matrices:           0.019     0.019   0.2% |
  Mix:                               0.064     0.064   0.6% |
  Multipole moments:                 0.001     0.001   0.0% |
  Pseudo density:                    0.106     0.015   0.1% |
   Symmetrize density:               0.091     0.091   0.8% |
 Hamiltonian:                        1.152     0.004   0.0% |
  Atomic:                            0.952     0.015   0.1% |
   XC Correction:                    0.937     0.937   8.2% |--|
  Calculate atomic Hamiltonians:     0.006     0.006   0.1% |
  Communicate:                       0.000     0.000   0.0% |
  Poisson:                           0.002     0.002   0.0% |
  XC 3D grid:                        0.189     0.189   1.7% ||
 Orthonormalize:                     0.003     0.000   0.0% |
  calc_s_matrix:                     0.001     0.001   0.0% |
  inverse-cholesky:                  0.000     0.000   0.0% |
  projections:                       0.001     0.001   0.0% |
  rotate_psi_s:                      0.000     0.000   0.0% |
Set symmetry:                        0.025     0.025   0.2% |
Other:                               8.970     8.970  79.0% |-------------------------------|
-----------------------------------------------------------
Total:                                        11.359 100.0%

Memory usage: 390.56 MiB
Date: Tue Apr  4 11:30:36 2023
~~~
{: .output}

~~~
fig, axes = plt.subplots(nrows=2, sharex=True)
axes[0].plot(nkpts, energies, 'o-')
axes[0].set_ylabel('energy / eV') 
axes[1].plot(nkpts, times, 'o-')
axes[1].set_ylabel('Calculation time / s')
axes[1].set_ylim([0, None])
axes[1].set_xlabel('number of k-points')
~~~
{: .python}

![](../fig/convergence_kpoint.png)

For this calculation, the computational cost per k-point is roughly linear, but the energy convergence is slow. We have multiple cores available, but how do we use them?

Create a script file with the following code, named "kpts_parallel.py"

~~~
from gpaw import GPAW, PW

import ase.build
from ase.parallel import parprint, world, paropen
from ase.utils.timing import Timer
import json

print("Hello from every process!")
parprint("Hello from one process!")

atoms = ase.build.bulk('Cu')

timer = Timer()
energies, times, nkpts = [], [], []

for k in range(3, 9):
    atoms.calc = GPAW(mode=PW(400), xc='PBE',
                      kpts={'size': [k, k, k],
                            'gamma': False},
                      txt='kpts_parallel.txt',
                      parallel={'kpt': True})
    timer.start(str(k))
  
    energy = atoms.get_potential_energy()
    
    timer.stop(str(k))
    energies.append(energy)
    times.append(timer.get_time(str(k)))
    nkpts.append(len(atoms.calc.get_ibz_k_points()))

with paropen('parallel_results.json', 'w') as fd:
    json.dump({'energies': energies,
               'times': times,
               'nkpts': nkpts},
              fd)
~~~
{: .python}

Run the script from the command-line with

~~~
gpaw -P 4 python kpts_parallel.py
~~~
{: .bash}

where `gpaw python` runs a special Python interpreter and the `-P 4` option runs this over 4 MPI tasks.

~~~
import json
with open('parallel_results.json', 'r') as fd:
    parallel_data = json.load(fd)

fig, axes = plt.subplots(nrows=3, sharex=True)
axes[0].plot(nkpts, energies, 'o-', label='serial')
axes[0].plot(parallel_data['nkpts'], parallel_data['energies'], 'o', label='parallel')
axes[0].set_ylabel('energy / eV')
axes[0].legend()

axes[1].plot(nkpts, times, 'o-', label='serial')
axes[1].plot(parallel_data['nkpts'], parallel_data['times'], 'o-', label='parallel')
axes[1].set_ylabel('Calculation time / s')
axes[1].legend()
axes[1].set_ylim([0, None])

axes[2].plot(nkpts, np.asarray(times) / parallel_data['times'], label='4 cores')
axes[2].set_ylabel('Speed-up factor')
axes[2].set_xlabel('number of k-points')
~~~
{: .python}

![](../fig/convergence_kpoint_parallel.png)

Hopefully, running in parallel did not change the results! We see that the parallel calculation is faster, but not 4 times faster.

The parallel Python process essentially runs 4 copies of the script. The GPAW part is smart enough to coordinate, but the copied processes can become obvious when we do things like `print`, and can cause data corruption when multiple processes try to write to the same file. `parprint` and `paropen` are provided in ASE to handle these common cases: for other scenarios it may be necessary to use logic with `ase.parallel.world`. For more information see the relevant [GPAW docs](https://wiki.fysik.dtu.dk/gpaw/documentation/parallel_runs/parallel_runs.html) and [ASE docs](https://wiki.fysik.dtu.dk/ase/ase/parallel.html)

## File-based calculators: parallel calculations

### Quantum Espresso

Quantum Espresso is a suite of programs for electronic structure calculations with a permissive Free Software license. Here we use the `pw.x` program is used for DFT calculations with a plane-wave basis set and pseudopotentials.

We need to find our pseudopotentials library and pass this information to the Calculator. In the virtual environment for these tutorials it can be found here:

~~~
from pathlib import Path
# Customise if necessary to espresso pseudopotential path
pseudo_dir = str(
    Path.home() / 'opt/espresso_sssp/SSSP_1.1.2_PBE_efficiency')
~~~
{: .python}

Unlike GPAW, we are going to call the program with MPI from within our regular Python interpreter. We define the mpi command when instantiating the calculator: the command might need to be tweaked for different machines with different parallel environments. This information is captured in a "Profile" object.

> ## Note  
> Profiles are a fairly new ASE feature and not yet used by all such Calculators. The existing way to manage these commands is by setting environment variables, e.g. ASE_ESPRESSO_COMMAND. Check the docs for each calculator to see how it is currently used.
{: .callout}

~~~
calculators.espresso import Espresso, EspressoProfile

profile = EspressoProfile(['mpirun', 'pw.x'])

calc = Espresso(profile=profile,
                pseudo_dir=pseudo_dir,
                kpts=(3, 3, 3),
                input_data={'control':  {'tprnfor': True,
                                         'tstress': True},
                            'system': {'ecutwfc': 50.}},
                pseudopotentials={'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'})
~~~
{: .python}

~~~
atoms = ase.build.bulk('Si')
atoms.calc = calc
atoms.get_potential_energy()
~~~
{: .python}

~~~
-310.1328387367529
~~~
{: .output}

Each Calculator has its own way of mapping keywords to the input syntax of the corresponding software code. The content of `input_data` may look a bit cryptic, but if you check `espresso.in` you can see there is a straightforward relationship between this data and the espresso input. In addition the structure information has been set, along with a few more automatically-generated keys.

> ## Exercise: Convergence testing
> Each Calculator has its own way of mapping keywords to the input syntax of the corresponding software code. The content of `input_data` may look a bit cryptic, but if you check `espresso.in` you can see there is a straightforward relationship between this data and the espresso input. In addition the structure information has been set, along with a few more automatically-generated keys.
>
> *Hint: to calculate the atomisation energy, you will need to compare the energy of the solid to a single atom in a large cell.*
>
> *Hint: you will need to calculate this property at several energies. Use Python functions and iteration constructs to avoid too much repetition.*
{: .challenge}

## Socket calculators

When computing the energy/forces for a series of related structures (e.g. during geometry optimisation), DFT codes typically re-use previously calculated wavefunctions to obtain a good starting estimate. When using a file-based ASE calculator, each calculation is essentially independent; not only do we miss out on wavefunction re-use, but we incur overhead as the code is "rebooted" (and memory re-allocated, pseudopotentials loaded, etc.)

This problem is addressed by codes implementing a socket interface: the calculation is initialised once and the code waits to receive a set of atomic positions. It calculates energies and forces, then awaits the next set of positions. Typically this is done by the "i-Pi protocol" developed to support calculation of nuclear quantum effects.

For the Espresso and Aims calculators in ASE, this is implemented as wrapper around the main calculator; for VASP and CP2K it works differently. Repeating our random-walk energy calculation from the Quippy/GAP example (in a smaller unit cell!):

~~~
from ase.calculators.socketio import SocketIOCalculator

calc = Espresso(profile=profile,
                pseudo_dir=pseudo_dir,
                kpts=(2, 2, 2),
                input_data={'control':  {'tprnfor': True,
                                         'tstress': True},
                            'system': {'ecutwfc': 40.}},
                pseudopotentials={'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'})

calc = SocketIOCalculator(calc=calc, unixsocket='new-random-walk-socket') 
~~~
{: .python}

This will produce a lot of output so we have not printed here!


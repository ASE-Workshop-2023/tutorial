---
title: "External calculators"
teaching: 15
exercises: 20
questions:
    - "How can I calculate standard properties using an external calculator?"
    - "How can I calculate standard properties using a machine learnt potential?"
    - "Why are socket calculators more computationally efficient?"
objectives:
    - "Calculate properties using an external calculator and machine learnt potential"
    - "Calculate properties using a socket interface"
keypoints:
    - "The `quippy` package provides a Python interface to a range of interatomic and tight-binding potentials"
    - "The workflow for external, file-based and built-in calculators is the same"
    - "Quantum Espresso is a suite of programs for electronic structure calculations "
    - "Socket interfaces allow more efficient communication and data re-use"
---

> ## Code connection
> In this episode we explore two external calculators: [Quippy](http://libatoms.github.io/QUIP/), which provides an interface to a range of interatomic and tight-binding potentials, including [Gaussian Approximation Potentials](https://libatoms.github.io/GAP/), and Socket calculations using the electronic structure code [Quantum Espresso](https://www.quantum-espresso.org/).
{: .callout}

### The `quippy` package provides a Python interface to a range of interatomic and tight-binding potentials

- Some calculators have interfaces which are not packaged with ASE, but available elsewhere.
- For example, the `quippy` package provides a Python interface to a range of interatomic and tight-binding potentials.
- In this episode we apply a machine-learning-based potential for Si.

> ## Getting the model and training data
> The documentation for Gaussian Approximation Potentials (GAP) links to a [few published GAP models](https://libatoms.github.io/GAP/data.html).
> To download and extract the data in a Jupyter notebook you can use bash.
> ~~~
> %%bash
>
> wget -q https://www.repository.cam.ac.uk/bitstream/handle/1810/317974/Si_PRX_GAP.zip
> unzip Si_PRX_GAP.zip
> ~~~
> {: .bash}
>
{: .callout}

> ## Note
> This requires a version of quippy that includes GAP. At the moment, `pip install` seems to work better than installing from conda-forge.
{: .callout}

### The workflow for external, file-based and built-in calculators is the same

- First, we import libraries and create an `Atoms` object; in this case, a silicon supercell.

~~~
from quippy.potential import Potential
si = ase.build.bulk('Si') * 4
~~~
{: .python}

- Second we attach a calculator, in this case a `Potential` object imported from `quippy`.

~~~
si.calc = Potential(param_filename='path/to/Si_PRX_GAP/gp_iter6_sparse9k.xml')
~~~
{: .python}

- Third, we calculate an energy; in this case we place this in a `for` loop and apply a random walk to the positions.

~~~
energies = []
for _ in range(10):
    si.rattle(stdev=0.01)
    energies.append(si.get_potential_energy())

fig, ax = plt.subplots()
ax.plot(energies, 'o-')
ax.set_xlabel('Random walk step')
ax.set_ylabel('Energy / eV')
~~~
{: .python}

<img src="../fig/energy_random_walk_plot.png" alt="Si random walk energy" width="400">

- This is a bit more expensive than EMT but still a lot cheaper than density-functional theory!
- A lot of work goes into developing a new potential, but with tools like quippy and ASE it is fairly easy for researchers to pick up the resulting model and apply it.

### Quantum Espresso is a suite of programs for electronic structure calculations

- Here we use the `pw.x` program for DFT calculations with a plane-wave basis set and pseudopotentials.
- First we need to find our pseudopotentials library and pass this information to the Calculator.
- The "Standard Solid-State Pseudopotentials" library is a general-purpose set available in "precision" and "efficiency" versions; for this tutorial we suggest to download the "efficiency" set [from this page](https://www.materialscloud.org/discover/sssp/table/efficiency).

~~~
from pathlib import Path
# Customise if necessary to espresso pseudopotential path
pseudo_dir = str(
    Path.home() / 'opt/espresso_sssp/SSSP_1.1.2_PBE_efficiency')
~~~
{: .python}

### Socket interfaces allow more efficient communication and data re-use

- Computing the energy/forces for a series of related structures (e.g. during geometry optimisation), DFT codes typically re-use previously calculated wavefunctions to obtain a good starting estimate.
- Using a file-based ASE calculator, each calculation is essentially independent; not only do we miss out on wavefunction re-use, but we incur overhead as the code is "rebooted" (and memory re-allocated, pseudopotentials loaded, etc.).
- This problem is addressed by codes implementing a socket interface.
- Here the calculation is initialised once and the code waits to receive a set of atomic positions. It calculates energies and forces, then awaits the next set of positions.
- Typically this is done by the "i-Pi protocol" developed to support calculation of nuclear quantum effects.
- For the Espresso and FHI-Aims calculators in ASE, this is implemented as wrapper around the main calculator; for VASP and CP2K it works differently.

~~~
from ase.calculators.qe import Espresso, EspressoProfile
from ase.calculators.socketio import SocketIOCalculator

profile = EspressoProfile(['pw.x'])

calc = Espresso(profile=profile,
                pseudo_dir=pseudo_dir,
                kpts=(2, 2, 2),
                input_data={'control':  {'tprnfor': True,
                                         'tstress': True},
                            'system': {'ecutwfc': 40.}},
                pseudopotentials={'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'})

calc = SocketIOCalculator(calc=calc, unixsocket='random-walk-socket')
~~~
{: .python}

> ## Exercise: Socket random walk calculation
> Repeat our random-walk energy calculation from the Quippy example using DFT implemented in Quantum Espresso.
>
> Hint: as this is a DFT calculation and we have only four compute cores, you will need to use a smaller unit cell.
{: .challenge}

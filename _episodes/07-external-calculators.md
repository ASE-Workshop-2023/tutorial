---
title: "Dynamics: geometry optimisation"
teaching: 20
exercises: 20
questions:
objectives:
keypoints:
    - "The `quippy` package provides a Python interface to a range of interatomic and tight-binding potentials"
    - "Socket calculators allow more efficient communication and data re-use"
---

> ## Code connection
> In this chapter we explore [Quippy](http://libatoms.github.io/QUIP/), which provides an interface to a range of interatomic and tight-binding potentials, including [Gaussian Approximation Potentials](https://libatoms.github.io/GAP/). We also briefly look at Socket calculations using [Quantum Espresso](https://www.quantum-espresso.org/).
{: .callout}

## The `quippy` package provides a Python interface to a range of interatomic and tight-binding potentials

- Some calculators have interfaces which are not packaged with ASE, but available elsewhere. 
- For example, the `quippy` package provides a Python interface to a range of interatomic and tight-binding potentials. 
- In this chapter we apply a machine-learning-based potential for Si.

> ## Getting the model and training data
> The documentation for Gaussian Approximation Potentials (GAP) links to a [few published GAP models](https://libatoms.github.io/GAP/data.html).
> Follow the "Si" link and download *Si_PRX_GAP.zip*, which contains the model and training data. 
> Move the .zip file somewhere you can find and extract the files.
{: .callout}

> ## Note
> This requires a version of quippy that includes GAP. At the moment, `pip install` seems to work better than installing from conda-forge.
{: .callout}

### We use the standard procedure for calculating properties

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

- Third, we calculate an energy; in this case we apply a random walk to the positions:

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

![](../fig/energy_random_walk_plot.png)

This is a bit more expensive than EMT but still a lot cheaper than density-functional theory! A lot of work goes into developing a new potential, but with tools like quippy and ASE it is fairly easy for researchers to pick up the resulting model and apply it.

## TODO

- geometry relaxation of silicon vacancy (see https://libatoms.github.io/QUIP/Tutorials/quippy-ase-interoperability.html)
- MD amorphous silicon? (https://pubs.acs.org/doi/10.1021/acs.jpclett.8b00902)
- exercise?

## Socket calculators allow more efficient communication and data re-use

- Computing the energy/forces for a series of related structures (e.g. during geometry optimisation), DFT codes typically re-use previously calculated wavefunctions to obtain a good starting estimate. 
- Using a file-based ASE calculator, each calculation is essentially independent; not only do we miss out on wavefunction re-use, but we incur overhead as the code is "rebooted" (and memory re-allocated, pseudopotentials loaded, etc.).
- This problem is addressed by codes implementing a socket interface.
- Here the calculation is initialised once and the code waits to receive a set of atomic positions. It calculates energies and forces, then awaits the next set of positions. 
- Typically this is done by the "i-Pi protocol" developed to support calculation of nuclear quantum effects.
- For the Espresso and Aims calculators in ASE, this is implemented as wrapper around the main calculator; for VASP and CP2K it works differently. Repeating our random-walk energy calculation from the Quippy/GAP example (in a smaller unit cell!):

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
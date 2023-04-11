---
title: "Dynamics: geometry optimisation"
teaching: 20
exercises: 20
questions:
objectives:
keypoints:
---

## Dynamics

## External library calculators

Some calculators are not packaged with ASE, but included in some other package. For example, the `quippy` package provides a Python interface to a range of interatomic and tight-binding potentials. Here we apply a machine-learning-based potential for Si.

The documentation for Gaussian Approximation Potentials (GAP) links to a few published GAP models: https://libatoms.github.io/GAP/data.html
Follow the "Si" link and download the *Si_PRX_GAP.zip* file containing the model and training data. Move the .zip file somewhere you can find and extract the files.

> ## Note
> This requires a version of quippy that includes GAP. At the moment, `pip install` seems to work better than installing from conda-forge.
{: .callout}

~~~
from pathlib import Path
from quippy.potential import Potential

gap = Potential(param_filename=str(Path.cwd() / 'Si_PRX_GAP/gp_iter6_sparse9k.xml'))
~~~
{: .python}

Let's build a Si supercell and apply a random walk to the positions:

~~~
si = ase.build.bulk('Si') * 4
si.calc = gap

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
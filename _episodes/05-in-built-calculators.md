---
title: "Built-in calculators"
teaching: 20
exercises: 20
questions:
objectives:
keypoints:
---

### Atoms objects can calculate  properties using an attached "Calculator"

Many software packages can be used to calculate properties from a set of atomic positions. In ASE, Atoms objects can try to calculate or fetch properties using an attached "Calculator".

Here we make tour of a variety of Calculators, highlighting some of the differences between them. A master list of the available Calculators can be found [here](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html).

### In-built calculators: simple potential models (EMT, LJ)

Much of the ASE documentation and tutorials makes use of the [inbuilt "EMT" calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/emt.html#pure-python-emt-calculator) - because it is convenient and "fast enough"!

This implements the Effective Medium Theory potential for Ni, Cu, Pd, Ag, Pt and Au. Some other elements are included "for fun", but really this is a method for alloys of those metals. For more information see this documentation and references for the (more efficient) implementation in the ASAP code: https://wiki.fysik.dtu.dk/asap/EMT

Considering the case of an infinite gold wire:

~~~
from ase import Atoms
from ase.calculators.emt import EMT

def make_wire(spacing: float = 2.5,
              box_size: float = 10.0) -> Atoms:

    wire = Atoms('Au',
                 positions=[[0., box_size / 2, box_size / 2]],
                 cell=[spacing, box_size, box_size],
                 pbc=[True, False, False])
    return wire

atoms = make_wire()
atoms.calc = EMT()
energy = atoms.get_potential_energy()
print(f"Energy: {energy} eV")
~~~
{: .python}

~~~
Energy: 0.9910548478768826 eV
~~~
{: .output}

~~~
import nglview

def show(atoms: Atoms) -> None:    
    view = nglview.show_ase(atoms)
    if any(atoms.pbc):
        view.add_unitcell()
    return view

show(atoms)
~~~
{: .python}

> ## Discussion
> Why did we need the parentheses () in the line `atoms.calc = EMT()`?
{: .discussion}

What can we do with an energy? We could look at how it varies with the atom spacing and fit a model.

~~~
import numpy as np
distances = np.linspace(2., 3., 21)

def get_energy(spacing):
    atoms = make_wire(spacing=spacing)
    atoms.calc = EMT()
    return atoms.get_potential_energy()

energies = list(map(get_energy, distances))
~~~
{: .python}

> ## Python tip
> if you need to apply a function to each element of some data, `map` can provide an elegant alternative to for-loops and list comprehensions!
{: .callout}

~~~
from numpy.polynomial import Polynomial
fit = Polynomial.fit(distances, energies, 3)
x = np.linspace(2., 3., 500)
~~~
{: .python}

~~~
%matplotlib inline
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
_ = ax.plot(distances, energies, 'o', label='calculated')
_ = ax.plot(x, fit(x), '-', label='cubic')
_ = ax.legend()
_ = ax.set_xlabel('Spacing / Å')
_ = ax.set_ylabel('Energy / eV')
~~~
{: .python}

![](./fig/energy_spacing_plot.png)

This somewhat resembles the Equation-of-State (EOS) curve for a solid. To see how you can fit a standard EOS model to a 3-D system and obtain an equilibrium volume, see the relevant [ASE tutorial](https://wiki.fysik.dtu.dk/ase/tutorials/eos/eos.html).

The EMT Calculator can also be used to obtain forces and unit cell stress:

~~~
print("Forces: ")
print(atoms.get_forces())

print("Stress: ")
print(atoms.get_stress())
~~~
{: .python}

~~~
Forces: 
[[0. 0. 0.]]
Stress: 
[ 0.00396458 -0.         -0.         -0.         -0.         -0.        ]
~~~
{: .output}

> ## Discussion
> Why are the forces exactly zero for this system?
{: .discussion}

Energy, forces and stress are standard "properties" in ASE, and we can check which properties are implemented by a particular calculator by inspecting the `implemented properties` attribute:

~~~
print(EMT.implemented_properties)
~~~
{: .python}

~~~
['energy', 'free_energy', 'energies', 'forces', 'stress', 'magmom', 'magmoms']
~~~
{: .output}

> ## Discussion
> Why do we *not* need to include parenthesis () here? Do we expect `EMT().implemented_properties` to work as well as `EMT.implemented_properties`?
{: .discussion}

It can be convenient to have e.g. forces attached to a particular Atoms object in this way, and will be used heavily by dynamics and optimizer routines in the next tutorial.

However, for Open Science purposes it is easier to store and share data that is not connected to a Calculator. (As we shall see later, Calculators might depend on a particular machine environment, memory state or software license.) We can request a standalone set of property data with e.g.:

~~~
properties = atoms.get_properties(['energy', 'forces', 'stress'])
print(properties)
~~~
{: .python}

~~~
(Properties({'energy': 0.9910548478768826, 'natoms': 1, 'energies': array([0.99105485]), 'free_energy': 0.9910548478768826, 'forces': array([[0., 0., 0.]]), 'stress': array([ 0.00396458, -0.        , -0.        , -0.        , -0.        ,
       -0.        ])})
~~~
{: .output}

This will not change even if the Atoms object is modified and properties are recalculated.

> ## Warning
> This is a new feature and does not yet work well for all calculators.
{: .callout}

#### Lennard-Jones

The classic [Lennard-Jones potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential) is implemented in `ase.calculators.lj`. You can set the $\epsilon$ and $\sigma$ parameters in the Calculator constructor:

~~~
from ase.calculators.lj import LennardJones

l = 4.1
atoms = Atoms('Xe2',
              positions=[[0., 0., -l / 2],
                         [0., 0., l / 2]],
              pbc=False)
atoms.calc = LennardJones(sigma=(4.1 / 2**(1/6)))

atoms.get_forces()
~~~
{: .python}

~~~

array([[ 0.00000000e+00,  0.00000000e+00, -6.49886649e-16],
       [ 0.00000000e+00,  0.00000000e+00,  6.49886649e-16]])
~~~
{: .output}

> ## Discussion
> Why are the forces so low at this geometry?
{: .discussion}

> ## Exercise: Lennard-Jones binding curve
> Try choosing a value of sigma and varying the distance between the atoms. Can you reproduce the classic plot of a Lennard-Jones binding curve?
{: .challenge}






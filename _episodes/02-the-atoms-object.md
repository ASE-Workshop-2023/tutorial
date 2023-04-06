---
title: "The Atoms Object"
teaching: 20
exercises: 20
questions:
- "How can I describe a molecule or crystal with code?"
- "How can I access simple structural properties?"
objectives:
- "List the key features of ASE"
- "Identify and describe the core ASE classes"
- "Understand how ASE fits within the larger atomistic modelling software ecosystem"
keypoints:
- "Molecules and materials are represented by the `Atoms` class"
- "There are various ways to visualise a structure"
- "To describe crystals we use the `pbc` keyword"
- "To get structural properties we use \"getter\" methods on the `Atoms` instance"

---

### Molecules and materials are represented by the `Atoms` class

We can define a molecule with lists of symbols and positions:

~~~
from ase import Atoms

d = 1.10
molecule = Atoms(['N', 'N'], positions=[(0., 0., 0.), (0., 0., d)])
~~~
{: .python}

For convenience we can compress the list of symbols to a chemical formula:

~~~
molecule = Atoms('N2', positions=[(0., 0., 0.), (0., 0., d)])
~~~
{: .python}

### There are various ways to visualise a structure

It can be useful to visualise our structure to make sure it is reasonable. 
`ase.visualize.view` provides a simple structure viewer in a floating window; this is quite useful when working on a Python script, but can be a little bit annoying when using a Jupyter notebook.

Try visualising the structure by running the code below; you can spin the molecule around with right-click-and-drag, and zoom with mouse wheel.

~~~
from ase.visualize import view
view(molecule)
~~~
{: .python}

Some alternative viewers are available for Jupyter notebooks; here we will use `nglview`. It should be pre-installed on the virtual machines for the workshop. In this viewer left-click-and-drag is used for rotation.

~~~
import nglview
nglview.show_ase(molecule)
~~~
{: .python}

### To describe crystals we use the `pbc` keyword

Many interesting systems are crystals, described by atomic positions in a periodic unit cell. There are two relevant keyword settings for an `Atoms` object: the unit cell itself and the periodic boundary conditions (PBC).

~~~
a = 5.387
crystal = Atoms('Zn4S4',
                scaled_positions=[[0., 0., 0.],
                                  [0., 0.5, 0.5],
                                  [0.5, 0., 0.5],
                                  [0.5, 0.5, 0.],
                                  [0.25, 0.75, 0.75],
                                  [0.25, 0.25, 0.25],
                                  [0.75, 0.75, 0.25],
                                  [0.75, 0.25, 0.75]],
               cell=[a, a, a],
               pbc=True)

def show(atoms: Atoms) -> None:    
    view = nglview.show_ase(atoms)
    if any(atoms.pbc):
        view.add_unitcell()
    return view

show(crystal)
~~~
{: .python}

Note that we used a few tricks to make writing in the structure a bit easier:

- The symbols were compressed to 'Zn4S4'
- Instead of working out positions in Angstrom, scaled_positions relative to lattice vectors were used
- The cell with just 3 values, so it is assumed to be cubic. In other cases we might use the full 3x3 matrix, e.g. ```cell=[[a, 0, 0], [0, a, 0], [0, 0, a]]```
We set ```pbc=True``` to indicate periodic boundary conditions in all directions. These can also be specified along each direction, e.g. ```pbc=[True, True, False]``` for a "slab" calculation with exposed surfaces.
- We also defined a Python function ```show()``` which will show Atoms including a unit cell with `nglview`. This will save us from writing three lines of code every time!

### To get structural properties we use "getter" methods on the `Atoms` instance

Now that we have some Atoms objects we can see what information is available from them. We call some "getter" methods on the `molecule` object, which is an instance of `Atoms`.

~~~
print("N2 positions")
print(molecule.get_positions(), end="\n\n")

print("N2 symbols")
print(molecule.get_chemical_symbols(), end="\n\n")

print("N2 masses")
print(molecule.get_masses(), end="\n\n")

print("N2 center of mass")
print(molecule.get_center_of_mass())
~~~
{: .python}

~~~
N2 positions
[[0.  0.  0. ]
 [0.  0.  1.1]]

N2 symbols
['N', 'N']

N2 masses
[14.007 14.007]

N2 center of mass
[0.   0.   0.55]
~~~
{: .output}

The first two attributes here are not surprising; they are the information we provided. The masses were not provided: when the `Atoms` was created, ASE found some standard values and included them in the data set.

If we like, we can override this default and include some isotopic effects:

~~~
d = 1.10
isotope = Atoms('N2',
                positions=[(0., 0., 0.), (0., 0., d)],
                masses=[13.006, 14.003])

print("13N-14N masses")
print(isotope.get_masses(), end="\n\n")

print("13N-14N center of mass")
print(isotope.get_center_of_mass())
~~~
{: .python}

~~~
13N-14N masses
[13.006 14.003]

13N-14N center of mass
[0.         0.         0.57030249]
~~~
{: .output}

The center of mass was not defined when the `Atoms` are created; it is derived from the other properties. So if we modify the masses using a setter, it should be recalculated correctly.

~~~
isotope = molecule.copy()
print("Center of mass before modifying masses:")
print(isotope.get_center_of_mass(), end='\n\n')

isotope.set_masses([13.006, 14.003])
print("Center of mass after modifying masses:")
print(isotope.get_center_of_mass())
~~~
{: .python}

~~~
Center of mass before modifying masses:
[0.   0.   0.55]

Center of mass after modifying masses:
[0.         0.         0.57030249]
~~~
{: .output}

### In a jupyter notebook we can get the "docstring" of a method or function by adding ? to the name

~~~
isotope.get_center_of_mass?
~~~
{: .python}

~~~
Signature: isotope.get_center_of_mass(scaled=False)
Docstring:
Get the center of mass.

If scaled=True the center of mass in scaled coordinates
is returned.
File:      ~/src/ase/ase/atoms.py
Type:      method
~~~
{: .output}

And we can get access to the available methods and properties with tab-completion. In a Jupyter notebook or IPython terminal, try:

~~~
crystal.[TAB]
~~~
{ :.python}

where `[TAB]` means "hit the TAB key". You should see that the `Atoms` object has a lot of features available!

Not all of them ll work until we start using the `Calculator` class, but the rest of this tutorial should include some useful ones.

> ## Exercise
> Use tab-completion and docstrings to explore the features of `Atoms`. In the 
> ZnS structure, find the distance between the first Zn atom and the four S
> atoms. Are some of them nearer than others?
{: .challenge}
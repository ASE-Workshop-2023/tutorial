---
title: "Building and Manipulating Atoms"
teaching: 20
exercises: 20
questions:
    - "How can I build structures without specifying all atomic positions?"
objectives:
    - "Build a molecule using the built-in database"
    - "Build a crystal using built-in crystal structure types"
    - "Build (optimal) supercell expansions"
    - "Use NumPy array operations to remove, add or swap atom(s)"
keypoints:
---

In this chapter we explore the [`ase.build` module](https://wiki.fysik.dtu.dk/ase/ase/build/build.html), which contains tools for building structures using parameters rather than detailed lists of positions.

### A set of simple molecules are 
Definitions for a set of simple molecules (the "G2" set, plus a few extra) are included with ASE. So in fact the easiest way to get an N2 molecule is

~~~
import ase.build
g2_n2 = ase.build.molecule('N2')
show(g2_n2)
~~~
{ .python}

TODO: look at the type to show it is `Atoms`

And it's just as easy to get a buckyball! This feature is very useful for testing and trying things out.

~~~
show(ase.build.molecule('C60'))
~~~
{: .python}

### Crystals

The equivalent tool for crystals is `ase.build.bulk`. This includes lattice parameters for some elemental reference states (the list is [in the code here](https://gitlab.com/ase/ase/-/blob/6ac638d0c699f7bc80c10a8dccb7d42eda011be2/ase/data/__init__.py#L578)),
but we can also use known lattice parameters to build structures. So we get copper "for free":

~~~
show(ase.build.bulk('Cu', cubic=True))
~~~
{: .python}

but for ZnS we have to provide a bit more information:

~~~
show(
    ase.build.bulk('ZnS',
                   crystalstructure='zincblende',
                   a=5.387,
                   cubic=True)
)
~~~
{: .python}

### Supercells

A compact notation can be used to create a repeated unit cell. Starting with a cubic unit cell of Si:

~~~
si = ase.build.bulk('Si', cubic=True)
show(si)
~~~
{: .python}

equal repetition in each direction is possible with an integer multiplication:

~~~
show(si * 4)
~~~
{: .python}

A 3-list or 3-tuple gives the expansion in each direction:

~~~
show(si * [2, 4, 1])
~~~
{: .python}

For more complex transformations, we need to pass a 3x3 matrix to `ase.build.make_supercell`. I f we start with the non-cubic _primitive_ cell

~~~
si_prim = ase.build.bulk('Si')
show(si_prim)
~~~
{: .python}

cubic supercells can be formed with the correct transformation matrix.

~~~
si_prim.cell
~~~
{: .python}

~~~
show(
    ase.build.make_supercell(si_prim, [[1, 1, -1], [1, -1, 1], [-1, 1, 1]])
)
~~~
{: .python}

How do we find such matrices? This case is a known "textbook" example, but we can also perform a numerical search to find the matrix giving the most cubic result. This can be useful when setting up supercell calculations to model dilute defects and maximise the distance between periodic images.

~~~
ase.build.find_optimal_cell_shape(si_prim.cell, 4, 'sc', verbose=True)
~~~
{: .python}

~~~
target metric (h_target):
[[1. 0. 0.]
 [0. 1. 0.]
 [0. 0. 1.]]
normalization factor (Q): 0.184162
idealized transformation matrix:
[[-1.  1.  1.]
 [ 1. -1.  1.]
 [ 1.  1. -1.]]
closest integer transformation matrix (P_0):
[[-1  1  1]
 [ 1 -1  1]
 [ 1  1 -1]]
smallest score (|Q P h_p - h_target|_2): 0.000000
optimal transformation matrix (P_opt):
[[-1  1  1]
 [ 1 -1  1]
 [ 1  1 -1]]
supercell metric:
[[5.43 0.   0.  ]
 [0.   5.43 0.  ]
 [0.   0.   5.43]]
determinant of optimal transformation matrix: 4
array([[-1,  1,  1],
       [ 1, -1,  1],
       [ 1,  1, -1]])
~~~
{: .output}

### Slicing and splicing

An Atoms object can be treated as a list of Atom objects, with numpy-like array slicing available.

~~~
for atom in crystal:
    print(atom.symbol, atom.position, atom.mass)
~~~
{: .python}

~~~
Zn [0. 0. 0.] 65.38
Zn [0.     2.6935 2.6935] 65.38
Zn [2.6935 0.     2.6935] 65.38
Zn [2.6935 2.6935 0.    ] 65.38
S [1.34675 4.04025 4.04025] 32.06
S [1.34675 1.34675 1.34675] 32.06
S [4.04025 4.04025 1.34675] 32.06
S [4.04025 1.34675 4.04025] 32.06
~~~
{: .output}

If we index multiple atoms, an Atoms object is returned:

~~~
zinc_indices = [i for i, atom in enumerate(crystal) if atom.symbol == 'Zn']
zinc_sublattice = crystal[zinc_indices]
show(zinc_sublattice)
~~~
{: .python}

Individual atoms can be appended:

~~~
from ase import Atom
composite = zinc_sublattice.copy()
composite.append(Atom('He', position=(1.34675, 4.04025, 4.04025)))
show(composite)
~~~
{: .python}

Or entire Atoms can be combined with `+`. (The first Atoms takes precedence to determine the cell etc.) So to obtain a distorted sphalerite cell, moving the S sublattice along the x coordinate relative to the Zn sublattice:

~~~
sulfur_sublattice = crystal[4:]
sulfur_sublattice.translate([.3, 0., 0.])
show(zinc_sublattice + sulfur_sublattice)
~~~
{: .python}

While Atoms is not exactly a regular Python object, it plays nicely with the delete operation. So to create a zinc-vacancy defect:

~~~
zinc_vacancy = crystal.copy()
del zinc_vacancy[0]
show(zinc_vacancy)
~~~
{: .python}

and to create antisite disorder, we can swap two positions from the positions array.

~~~
antisite = crystal.copy()
antisite.positions[[0, 4]] = antisite.positions[[4, 0]]
show(antisite)
~~~
{: .python}



> ## Exercise: Water animation
> Create an animation of a 
> water molecule being wrapped in 
> a C60 cage (or something even 
> cooler!)
> 
> Hints:
> - The GIF animation will need to be generated with a list of Atoms objects
> - Molecules can be combined with +
> - To get the wrapping effect we need to keep adding atoms that are near to the atoms already added
> - To avoid writing too much repetitive code, use Python's looping tools
> 
> > ## Solution
> > 
> > ~~~
> > water = 
> > ase.build.molecule('H2O')
> > water.center()
> > 
> > bucky = ase.build.molecule('C60')
> > bucky.center()
> > 
> > start_atom = 36
> > distances = bucky.get_all_distances()[start_atom]
> > sorted_bucky_indices = sorted(enumerate(distances),
> >                               key = (lambda x: x[1]))
> > sorted_bucky_indices = [i for i, _ in sorted_bucky_indices]
> > sorted_bucky = bucky[sorted_bucky_indices]
> > 
> > frames = [water.copy()]
> > 
> > for i in range(len(sorted_bucky)):
> >     frames.append(water + sorted_bucky[:i + 1])
> > 
> > from ase.io.animation import write_gif
> > _ = write_gif('wrapped_molecule.gif', frames)
> > ~~~
> > {: .python}
> {: .solution}
{: .challenge}


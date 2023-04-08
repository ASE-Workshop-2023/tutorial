---
title: "Reading and Writing Atoms"
teaching: 10
exercises: 15
questions:
- "How do I read structures from a file?"
- "How do I write structures to a file?"
- "Which file formats does ASE support?"
objectives:
- "Read in a structure from a file"
- "Write your structure to a file"
keypoints:
---

### ASE can read and write a variety of file formats

Even with the shortcuts above, writing out structures in Python syntax can be a bit cumbersome. There are many established file formats for this data, and ASE contains read/write functions for some of them.

For example, we might obtain a structure from an online database. Here is the [Materials Project entry for sphalerite](https://materialsproject.org/materials/mp-10695); use the "export as" button on the structure visualiser to obtain a .cif file.

~~~
import ase.io
from pathlib import Path

# You may need to change the path to match the location/file you downloaded
imported_crystal = ase.io.read(Path.home() / "Downloads/ZnS.cif", format='cif')

show(imported_crystal)
~~~
{: .python}

Using the virtual desktop, have a look at this .cif file in a text editor. CIF is quite a complicated format because it is designed to hold a lot of data relevant to crystallography. By comparison, the `FHI-aims` quantum chemistry code has a very simple structure format. Let's write the ZnS structure in this format instead:

~~~
ase.io.write('geometry.in', imported_crystal, scaled=True)
~~~
{: .python}

Taking a peek at this file:

~~~
%%bash
cat geometry.in
~~~
{: .sh}

~~~
#===============================================================================
# Created using the Atomic Simulation Environment (ASE)

# Mon Apr  3 11:11:50 2023

#=======================================================
lattice_vector 5.3873657499999998 0.0000000000000000 0.0000000000000000 
lattice_vector 0.0000000000000000 5.3873657499999998 0.0000000000000000 
lattice_vector 0.0000000000000000 0.0000000000000000 5.3873657499999998 
atom_frac 0.0000000000000000 0.0000000000000000 0.0000000000000000 Zn
atom_frac 0.5000000000000000 0.5000000000000000 0.0000000000000000 Zn
atom_frac 0.5000000000000000 0.0000000000000000 0.5000000000000000 Zn
atom_frac 0.0000000000000000 0.5000000000000000 0.5000000000000000 Zn
atom_frac 0.2500000000000000 0.2500000000000000 0.2500000000000000 S
atom_frac 0.2500000000000000 0.7500000000000000 0.7500000000000000 S
atom_frac 0.7500000000000000 0.7500000000000000 0.2500000000000000 S
atom_frac 0.7500000000000000 0.2500000000000000 0.7500000000000000 S
~~~
{: .output}

Note that
- we didn't specify the FHI-aims format; it was correctly guessed from the filename. (This also works for ase.io.read.)
- we added a `scaled=True` option to write in fractional coordinates; by default this writer uses Cartesian coordinates.

How does this work? There are lots of relevant functions in modules under `ase.io`: `ase.io.read()` and `ase.io.write()` will automatically dispatch to these functions. If we look at `ase.io.write?` we get the signature

~~~
ase.io.write(
    filename: Union[str, pathlib.PurePath, IO],
    images: Union[ase.atoms.Atoms, Sequence[ase.atoms.Atoms]],
    format: str = None,
    parallel: bool = True,
    append: bool = False,
    **kwargs: Any,
) -> None
~~~
{: .output}

`**kwargs` means "all remaining keyword arguments"; this allows extra options like `scaled=True` to be passed to writers that understand them.

More information about the supported formats can be found at https://wiki.fysik.dtu.dk/ase/ase/io/io.html and a summary can be produced from the command-line 

~~~
%%bash
ase info --formats
~~~
{: .bash}

To find the extra supported `**kwargs`, look at the documentation of the lower-level functions. For example we can find the `scaled=True` option for FHI-aims at https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#ase.io.aims.write_aims

> ## Exercise: Converting structures
> Import a structure file relevant to your own research, 
> and write it to a different format. See what keywords 
> are available for your favourite formats; for example, 
> VASP users are likely to be interested in using  
> `vasp5=True`.
>
{: .challenge}
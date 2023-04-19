---
title: "What is ASE?"
teaching: 10
exercises: 0
questions:
- "How can ASE help me in my research?"
- "In what way is the ASE code structured?"
objectives:
- "List the key features of ASE"
- "Identify and describe the core ASE classes"
- "Describe how ASE fits within the larger atomistic modelling software ecosystem"
- "Understand the difference between class, instance and method in object oriented programming"
keypoints:
- "ASE is a Python library for atomistic modelling"
- "ASE uses several design patterns from object oriented programming (OOP)"
- "The `Atoms` class is used to represent molecules and materials"
- "The `Calculator` class calculates basic properties of an `Atoms` object"
- "ASE plays nicely with a variety of atomistic modelling tools"
---

### ASE is a Python library for atomistic modelling

- The atomic simulation environment (ASE) is Python library for _atomistic modelling_; it allows you to set up, run and analyse atomistic simulations using Python scripts. 
- Systems are described by a set of atomic positions, and calculations are derived from these.
- There are two core classes in ASE: `Atoms` and `Calculator`.

> ## NumPy
> ASE makes strong use of [NumPy arrays](). NumPy arrays have a host of useful features such as
> array slicing, and ensure that there is efficient performance even for thousands or millions
> of atoms. If you are unfamiliar with NumPy arrays there is a 
> brief overview on the [ASE website]().
{: .callout}

### ASE uses several design patterns from object oriented programming (OOP). 

- You do not need to be familiar with the intricacies of OOP to start coding with
ASE, however some basic terminology is useful.
	- An *object* is an combination of data alongside operations for manipulating and returning information on that data.
	- A *class* is a category of objects - for example, `Composer()`.
	- An *instance* is an object that belongs to a class - for example, the instance `beethoven` belongs to the class `Composer`.
	- A *method* is a fuction that belongs to an object - for example, the class `Composer`  contains the function `count_symphonies()`.

- To use ASE you will need to create an instance of a pre-defined class. 
- To take an example most will be familiar with, lets assume you have a class called `Composer()`.
- We can create an instance of that class called `beethoven`.

~~~
beethoven = Composer(birth_year=1770)
~~~
{: .python}

- This then allows us to access the methods and data associated with `beethoven` and the `Composer` class more generally:

> ~~~
> print(beethoven.count_symphonies())
> print(beethoven.birth_year)
> ~~~
> {: .python}
> ~~~
> 9
> 1770
> ~~~
> {: .output}
{: .callout}

### The `Atoms` class is used to represent molecules and materials

- When working with ASE, molecules and materials are represented by an `Atoms` object. 
- `Atoms` represents (at minimum) a collection of atoms of any chemical species and associated positions

### The `Calculator` class calculates basic properties of an `Atoms` object

- A `Calculator` object can be attached to an `Atoms` object. 
- The `Calculator` takes the atomic numbers and positions from `Atoms` and calculates basic properties such as energy or forces. 
- `Calculator` objects can themselves be split into three types: 
	- in-built calculators that run the simulation within the same Python interpreter process. - file-based calculators that run the simulation as a sub-process, with communication mediated through input and output files.
	- calculators that run the simulation as a sub-process, with communication via pipes. Each approach has its own advantages and limitations which will be outlined later in the course.

> ## Further information
> For further information about code structure, implementation and features we recommend
> reading [this topical review](https://dx.doi.org/10.1088/1361-648X/aa680e) from Larsen et al.
{: .callout}

### ASE plays nicely with a variety of atomistic modelling tools

- The modular design of ASE allows it to play nicely with ("interoperate") with many different atomistic simulation codes. 
- This includes codes implementing density functional theory, (semi-)empirical methods, tight-binding models and classical interatomic potentials. 
- For a full list of supported codes please see the [ASE documentation]().

> ## Course structure
> The structure of this course is designed to reflect the structure of ASE. 
> First we introduce the `Atoms` object (chapters 2-4) and `Calculators` object (chapters 5-6). 
> We then outline how to simulate more complex and computationally demanding systems by parallelising over multiple compute cores (chapter 7). 
> Finally we apply our understanding to three common tasks: molecular dynamics, geometry optimisation and electronic structure calculations (chapter 9).
> As the course progresses we move from smaller code snippets introducing the core ASE objects, to extended pieces of code for more complex tasks.
> Throughout the course you will see boxes like this. They will provide Python hints, signposting to other resources, discussion prompts or exercises for you to complete.
{: .callout}






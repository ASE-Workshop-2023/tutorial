---
title: "Extra topics"
questions:
    - "Why are socket calculators more computationally efficient?"
objectives:
    - "Calculate properties using a socket interface"
keypoints:
    - "Socket interfaces allow more efficient communication and data re-use"
---

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

> ## Quantum Espresso setup
> The Quantum Espresso calculator was introduced in [a previous tutorial](../08-calculations-in-parallel/index.html),
> if you haven't done this tutorial yet, look there for some setup information.
{: .callout}

> ## Exercise: Socket random walk calculation
> Repeat our random-walk energy calculation from the [Quippy example](../07-external-calculations/index.html) using DFT implemented in Quantum Espresso.
>
> Hint: as this is a DFT calculation and we have only four compute cores, you will need to use a smaller unit cell.
{: .challenge}

> ## Exercise: Geometry optimisation
> For a small unit cell of your choice, try performing geometry optimization with the QuasiNewton optimizer.
> How does the performance compare between using Espresso as a FileIOCalculator and as a SocketIOCalculator?
{: .challenge}

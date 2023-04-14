---
title: "Electronic structure"
teaching: 20
exercises: 20
questions:
objectives:
keypoints:
---

> ## Code connection
> In this chapter we will use [Quantum Espresso](https://www.quantum-espresso.org/) to calculate the electronic structure (bandstructure and density of states) of silicon.
{: .callout}

### We can use the `tocell()` method to build a FCC lattice

- The first step, as always, is to import libraries and build an `Atoms` object
- Here we create an `FCC` class instance and then convert this to a `Cell` object

~~~
from ase import Atoms
from ase.lattice import FCC
from ase.units import Bohr

si_fcc = Atoms(symbols='Si2',
               cell = FCC(a=10.20*Bohr).tocell(),
               scaled_positions=[(0.0, 0.0, 0.0),
                                 (0.25, 0.25, 0.25)],
               pbc=True)
~~~
{: .python}

## setup basic PW calculator for scf calculation

1. The main input variables: 
   * we specify in control that we want to do a self-consistent calculation, a prefix to identify the calculation, the directory where the program saves the data files, pseudopotentials' location;
   * in system we specify the two cutoff for the FFT grid, their values are read from the dataset provided with the SSSP pseudopotential library;
   * in electrons we select the iterative diagonalization algorithm and the convergence threshold. 

~~~
from pathlib import Path
import json
pseudo_dir = str(Path.home() / 'SSSP-1.2.1_PBE_efficiency')
outdir = Path('./si_fcc/out).absolute() 
with open(str(pseudo_dir / 'SSSP_1.2.1_PBE_efficiency.json')) as fj:
  dj = json.load(fj) 
ecutwfc = dj['Si']['cutoff_wfc']
ecutrho = dj['Si']['cutoff_rho'] 
control = {'calculation': "scf",
           'prefix': "si_fcc",
           'outdir': f"{outdir}",
           'pseudo_dir': pseudo_dir,
           'tprnfor': True,
           'tstress': True}
system = {'ecutwfc': ecutwfc,
          'ecutrho': ecutrho}
electrons = {'diagonalization':'davidson',
             'conv_thr': 1.e-8}
~~~
{: .python}

2. The pseudopotentials for each atomic species, taken from the SSSP library
3. The MP mesh for kpoints 
4. The offset of the MP mesh 
5. The command for starting the calculations 

~~~
atomic_species = {'Si':dj['Si']['filename']}
kpts = [8,8,8] #8X8X8 MP k-point mesh 
koffset = [1,1,1] 

~~~
{: .python}

6. Setup how pw.x is run and save the info in EspressoProfile instances:
   we setup two profiles:
   *  one where we run with 4 MPI ranks with the default parallelization
   *  one where we run with 4 MPI ranks using pools parallelization, distributing k-points in 4 pools 

~~~
from ase.calculator.espresso import EspressoProfile
qepath = str(Path.home() / 'qe-7.2')
plain_argv = ['mpirun', '-np', '4', f"{qepath}/bin/pw.x"]
pools_argv = ['-nk', '4']
profile = EspressoProfile(plain_argv)
profile_4pools = EspressoProfile(plain_argv + pools_argv)    
~~~
{: .python}

7. We create the Espresso calculator for the scf calculation. We initialize it with a dedicated directory where the calculator will save the input and the standard-output file of pw.x. As we are going to run few calculations one after the other we will use a specific directory for each of them so that we don't overwrite the input and output files. The binary files that have to be read by the successive calculations will instead be saved in the same directory specified by the outdir input variable defined above. 

~~~
from ase.calculators.espresso import Espresso 
scf_directory = Path('./si_fcc/scf/).absolute()

scf_calc = Espresso(directory=scf_directory,
                profile=profile,
                input_data={'control': control,
                            'system': system,
                            'electrons': electrons},
                pseudopotentials=atomic_species,
                kpts=[8,8,8],
                koffset=[1,1,1])
~~~
{: .python}

8. We now use the calculator to run an scf calculation for our system.
~~~
si_fcc.calc = scf_calc
props_scf = si_fcc.get_properties(['energy'])
round(props_scf['energy'],8) 
~~~
{: .python}

~~~
-215.69011555
~~~
{: .output}

## Bandstructure

Setup the calculator for computing the bands. 
1. First we generate the k-point path using ASE bandpath method:

~~~
kpointpath = si_fcc.cell.bandpath(path='LGXWKG', density=20)
kpointpath.special_points
~~~
{: .python}

~~~
{'G': array([0., 0., 0.]),
 'K': array([0.375, 0.375, 0.75 ]),
 'L': array([0.5, 0.5, 0.5]),
 'U': array([0.625, 0.25 , 0.625]),
 'W': array([0.5 , 0.25, 0.75]),
 'X': array([0.5, 0. , 0.5])}
~~~
{: .output}

7. Then we setup the new calculator for bands:
   * we set 'bands' as 'calculation in this calculator
   * we add some conduction bands nbnd 4 --> 12  
   * we use the pools parallized profile 
   * we raise e bit the convergence threshold 1.e-8 --> 1.e-6

~~~
control.update({'calculation':"bands"})
system.update({'nbnd':12})
electrons.update({'conv_thr': 1.e-6})
bands_directory = Path('./si_fcc/bands').absolute()
bands_calc = Espresso(directory=bands_directory,
                  input_data={'control': control,
                              'system': system,
                              'electrons':electrons},
                  pseudopotentials=atomic_species, 
                  profile=profile_4pools,
                  kpts=kpointpath) 
~~~
{: .python}

8. We run the bands calculation with the new calculator. 

~~~
si_fcc_bands = si_fcc.copy()
si_fcc_bands.calc = bands_calc
results_bands = si_fcc_bands.get_properties(['eigenvalues']) 
~~~
{: .python}

9. We now use ASE's BandsStructure class to plot the bands:
   * we extract the band's  energies from the results dictionary of the bands_calc calculator;
   * we use as reference level the Fermi level that has been computed in the scf calculation;
   * as path argument we use again the kpointspath generated above. 

~~~
from ase.spectrum.band_structure import BandStructure
si_fcc_band_structure = BandStructure(path=kpointpath, energies=results_bands['eigenvalues'], reference=si_fcc.calc.get_fermi_level())
bandplot = band_structure.plot(emin=-6, emax=15)
bandplot.figure.savefig('band_structure_Si.png')
~~~
{: .python}
![](../fig/band_structure_Si.png)

10. Now we setup the calculator for computing the DOS, in this case we will perform first an nscf calculation with a symmetrized, denser  MP grid. 
   * we need to change verbosity to 'high' because we are using more than 100 kpoint and the printout of energies in the standard output wound otherwise be suppressed. 
   * we further increase the number of bands to avoid artifacts at higher energies in our plot. 

~~~
control.update({'calculation':'nscf',
                'verbosity': 'high'})
system.update({'nbnd': 20})
nscf_directory = Path('./si_fcc/nscf').absolute() 
nscf_calc = QEpw(directory =nscf_directory,
                 profile = profile_4pools,
                 input_data={'control': control,
                             'system': system,
                             'electron': electrons}, 
                 pseudopotentials=atomic_species, 
                 kpts=[16,16,16],
                 koffset=[0,0,0]) 
si_fcc_nscf = si_fcc.copy()
si_fcc_nscf.calc = nscf_calc
nscf_results = nscf_calc.get_properties(['eigenvalues'])
~~~ 
{: .python}

11. We use ASE's DOS class to compute the Density of States, the initiatialization of the DOS class extracts all information needed directly from the nscf calculator that we pass as argument to the DOS init method: 

~~~
from ase.dft.dos import DOS 
from matplotlib  import pyplot as plt

dos = DOS(nscf_calc, npts = 200,width=0.2) 
energies=dos.energies
dosvals = dos.get_dos() 
plt.plot(energies, dosvals)
plt.axis([-13, 15, 0, 5])
plt.savefig('DOS_Si.png')
~~~
{: .python}

![](../fig/DOS_Si.png)

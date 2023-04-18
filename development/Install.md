# Install

- Install Anaconda

- Install Xcode Command Line Tools
>  xcode-select --install

- Install Homebrew
> /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
> echo 'export PATH=/usr/local/bin:$PATH' >> ~/.bash_profile

- Install ASE and GPAW dependencies

brew install libxc
brew install open-mpi
brew install fftw

- Conda:

```
conda create --name ase-tutorial
conda activate ase-tutorial
# https://www.vasp.at/forum/viewtopic.php?t=18626
conda install -c conda-forge mopac nglview qe ipywidgets=7
conda install numpy scipy matplotlib jupyter
conda install pip
pip install ase gpaw quippy-ase
```

run 

- extensions
```
jupyter-nbextension enable nglview --py --sys-prefix
```

- GPAW dataset

https://wiki.fysik.dtu.dk/gpaw/setups/setups.html#installation-of-paw-datasets

- navigate to where you want to install and: 
`gpaw install-data ./`. CLick y.

- QE dataset?

- Quippy model?



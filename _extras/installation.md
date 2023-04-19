---
title: Installation
---

Note! Not yet complete....

Note that installation will depend somewhat on the operating system and computer being used. 
There are notes below as guidance; follow them at your own risk! They have been tested on macOS Big Sur 11.5.2.

### 1) Install Anaconda or Miniconda

We will use conda to manage our computer environment. This comes in two forms; Anaconda or Miniconda. 
If you are new to programming, you may prefer the former. 
For Miniconda see [here](https://docs.conda.io/en/latest/miniconda.html), for Anaconda see [here](https://www.anaconda.com/download).

### 2) Install Xcode Command Line Tools (Mac only)

This is for Mac only.
There is a strong chance these may already be installed, but if not: 

```
xcode-select --install
```

### 3) Install Homebrew

We need Homebrew for the next step.

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
echo 'export PATH=/usr/local/bin:$PATH' >> ~/.bash_profile
```

### 4) Install ASE and GPAW dependencies

We use the Homebrew `brew` command to install dependencies for ASE and GPAW.

```
brew install libxc
brew install open-mpi
brew install fftw
```

### 5) Conda:

Use conda to install ase and more.

```
conda create --name ase-tutorial
conda activate ase-tutorial

conda install -c conda-forge mopac nglview qe gpaw jupyter quippy-ase ipywidgets==7.6.5
conda install pip
pip install git+https://gitlab.com/ase/ase.git@master
```

### 6) Enable extensions

This is required for visualising structures in a Jupyter Notebook using nglviewer.

```
jupyter nbextension install widgetsnbextension --py --sys-prefix 
jupyter nbextension enable --py --sys-prefix widgetsnbextension 
jupyter-nbextension enable nglview --py --sys-prefix
```

### 7) Files and environment variables

You will also need to download some data files and set environment variables. These are detailed in the relevant chapters.

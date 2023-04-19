---
title: "Setup"
teaching: 10
exercises: 0
questions:
    - "How do I setup my computing environment for the Daresbury workshop tutorials?"
objectives:
    - "Access virtual environment for the tutorials"
keypoints:
    - "We will use a STFC cloud environment"
    - "Use your username given during registration"
    - "It is possible to create a new instance" 
    - "There are several useful tools installed"
    
---

### We will use a STFC cloud environment 

- We will use STFC cloud environment https://training.jupyter.stfc.ac.uk. 
- It runs a custom baked docker image of Ubuntu Jammy Jellyfish


### Use your username given during registration

- Go to training.jupyter.stfc.ac.uk into your browser
  - Click on **Signup** 
  - Provide the username given during registration
  - Choose password 
  - Click **Create User**

<img src="../fig/step_0.png" alt="landing page screen" width="1200">

<img src="../fig/step_1.png" alt="Signup page" width="1200">

- Authorization happens behind the scenes; if successful you will see something like the image below.

<img src="../fig/step_2.png" alt="Signup page" width="1200">

- Login with the credentials from above

<img src="../fig/step_3.png" alt="Login Screen" width="1200">

- If all ok you shall see something like the image below.

<img src="../fig/step_4.png" alt="Login success" width="1200">

### It is possible to create a new instance 

- Sometimes you might need to create a new instance; for example, if something goes wrong or we need to use an updated image.
- In this case, follow the steps below.

1) Go to the hub settings: `File -> Hub Control Panel`

<img src="../fig/step_5.png" alt="hub control panel settings" width="1200">

2) Stop the instance: click on on the `Stop My Server` button 
3) Logout: Click `Logout`.

<img src="../fig/step_6.png" alt="stop server and logout" width="1200">

4) Start the instance

<img src="../fig/step_7.png" alt="start a new instance" width="1200">

5) Choose the instance: (Test) - ASE Image

<img src="../fig/step_8.png" alt="stop server and logout" width="1200">

### There are several useful tools installed

- In addition to the packages required for the tutorial, there are several other useful tools installed.

#### Browsers

**Mozilla Firefox** is installed on the machine.

#### Compilers

The **GNU** toolchain is used throughout the summer school and are available at the unix prompt.

* **gcc**: the C/C++ compiler
* **gfortran**: the fortran compiler (it assumes f77 and f95 for ``*``.f and ``*``.f90 respecively). Some of the codes may be in fixed format requiring the compiler flag -ffixed-form.
* **python3** is available on the machine, use python3, be aware that python will give you python2.

#### Plotting Packages

Two graphics packages are available for plotting graphs: **gnuplot** and **xmgrace**. You can also use matplotlib from python.

#### Molecular Graphics Packages

**Jmol, VESTA, AVOGADRO, VMD and xcrysden** is also available. In order to use Jmol type *jmol* on the command line.

#### Editors

There are several editors available. You should choose whichever you are confortable with.

* **vi** the venerable UNIX screen mode editor.
* **vim** the improved venerable UNIX screen mode editor.
* **emacs** probably the commonest full-screen UNIX editor.
* **gedit** gui editor

#### Terminals

When one refers to terminal, console or command line, usually means a shell window. Gnome Terminal, xterm and uxterm are available,
You can click on the terminal icon to get one in the desktop or in the jupyter hub.

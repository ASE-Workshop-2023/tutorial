---
title: "Dynamics"
teaching: 20
exercises: 20
questions:
objectives:
keypoints:
---

> ## Code connection
> In this chapter we will use the [EMT](), [MOPAC]() and [Quippy]() calculators to perform molecular dynamics and geometry optimisations.
{: .callout}

### Energies and forces can be used to update structures

- In the previous tutorials we created `Atoms` objects and used Calculators to get properties including energies and forces - but we didn't do very much with this information.
- A very useful thing to do with this information is to update our strucuture!
- ASE includes algorithms for molecular dynamics, geometry optimization and global optimization.
- It is a useful toolkit for experimentation and method-development in this area, or development of multi-step pipelines.
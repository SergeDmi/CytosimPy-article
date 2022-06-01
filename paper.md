---
title: 'PyCytosim : a python interface for Cytosim'
tags:
  - C++
  - Python
  - Langevin
  - Biology
  - Simulation
authors:
  - name: Serge Dmitrieff^[*]
    orcid: 0000-0002-1362-5670
    affiliation: "1"
  - name: François Nédélec
    orcid: 0000-0002-8141-5288
    affiliation: "2"
affiliations:
 - name: Institut Jacques Monod CNRS UMR7592 and Universite Paris Diderot, 75205 Paris Cedex 13 France
   index: 1
 - name: Sainsbury Laboratory, University of Cambridge, Great Britain
   index: 2
date: 1 May 2022
bibliography: paper.bib
---

# Summary

Cytosim is a state-of-the art simulator for semi-flexible polymer dynamics, and cytoskeletal networks in particular. It is a fast C++ implementation of the Langevin equations (for mechanics) and Gillespie algorithm (for biochemical reactions), that encompases a broad range of biological objects : cytoskeleton fibers, molecular motors, and other objects.  
We developped a python interface with cytosim, that allows real-time access to native cytosim objects. Thus is it possible to not only access the simulation state - e.g. to analyze it - but also to modify the simulation state using python commands.


# Statement of need

Since its 2007 release, cytosim has attracted a lot of use [CITE] given its capacity to rapidly simulate large networks. The input interface to the simulation lies in a configuration file, interpreted by the simulation engine, calling the classes defined in the C++ code. Beyond its graphical interface, the output of the simulation is done by reporting in text files, that are typically later analyzed by users through data analysis scripts, usually written in python or matlab. 

Through our experience with the simulation and its users, we realized that this workflow offered several caveats. To run simulations for a range of parameter values, a list of configuration files had to be generated. To do this automatically, a utility (preconfig.py) was provided [CITE] but it did not cover all use cases. The configuration files also do not scripting, such as dynamic events (e.g. finish the simulation if there are more than 100 filaments). The output to text files was limited in scope, and any new desired reporting had to be coded in C++, which offered several challenges.  Lastly, the user could not add or change features without going through the lengthy and technical process of creating a new C++ class, implementing the features in C++, and consequently altering the compilation process. Given the complexity of biological processes, it was inevitable that users would require novel features, and were not able to implement them in the simulation.

Therefore, we developed PyCytosim, a python interface for Cytosim. This allows user to directly access the native cytosim objects (C++ class instances) in python. Thus objects can be created, the simulation can be altered, and the results can be analyzed at runtime, in python. Thus, a user can seamlessly run and analyze simulations for ranges of parameters, while keeping a dynamic control on simulations. 

# Implementation

PyCytosim is implemented in C++17 and relies on the core of Cytosim core with hardly any code change. It adds C++17 python bindings using the header library pybind11. Overwhelmingly, PyCytosim just provides access to native Cytosim objects, and their member functions. 

## Cytosim

Cytosim solves the Langevin equation for a set of objects [CITE]. Since Cytosim uses an implicit scheme to solve this equation, the forces are usually computed through a local quadratic approximations :  
$$ f(x+\delta x) = f(x) + \delta x \partial f/\partial x \big|_x $$  
The vectors $f$ and matrices $\partial f/\partial x$ are stored as a large vector and matrix for all the objects of the simulation. Additionally, biochemical reactions can be defined for some objects. Those are implemented through the Gillespie algorithm.

The objects themselves are instances of C++ classes. Cytosim classes follow the object-oriented paradigm and uses extensively class inheritence. For example, a cytosim filament is an instance of class Fiber, that derives from class Meca (a set of points on which forces can be implemented), that itself derives from Object. One or several objects are associated to an instance of Property, that contains the properties of these objects, e.g. the rigidity of a fiber, the stiffness of a crosslinker, etc. 

Thus a Cytosim simulation essentially stores a (usually large) collection of objects and a (usually small) collection of object properties, as well as a force vector $f$ and a force derivative matrix $\partial f/\partial x$ stored in an instance of the class Meca. Object mechanics, and mechanical interactions between objects will add elements to $f$ and  $\partial f/\partial x$ through accessory functions implemented in the class Meca. 

## PyCytosim

PyCytosim is designed to provide a python interface to a running cytosim simulation - and thus can access cytosim objects in memory at runtime. Most cytosim objects, such as fibers (instances of class Fiber) are passed to python by reference, so that the destruction of the python variable does not cause the destruction of the cytosim object in memory. The python class is accessible in the PyCytosim module as *cytosim.Class* (e.g. *cytosim.Fiber*).  

C++ classes have member variable and functions, that can be public or private. Most public functions of Cytosim classes have been bound or translated to the python classes. For exemple, Fiber::property() is a C++ member function of class Fiber that returns (a pointer to) the fiber's property. In python it is available as Fiber.property() and returns a reference to the property itself.
 Because it is possible to directly access native Cytosim objects, an effort was put to bind as many cytosim objects, and their members, as possible. Thus these objects can be directly manipulated in python, as well as passed as arguments to the (python-bound) C++ member functions. For instance a instance of the python class *cytosim.Fiber* can be passed to a binding of any C++ function that takes a Fiber as input, e.g. *Simul::Fibers::delete(fiber)*. 

For some variable types, a conversion step Python->C++ or C++->Python is necessary.For example the cytosim C++ class *Vector* is converted to a numpy array (and other way round), since using a numpy array is more natural in python than Cytosim's vector. In these cases, the conversion is usually performed in the binding rather than by the user. Therefore, for the user, PyCytosim behaves as a regular python module where (nearly) no conversion has to be explicitely performed.  

A cytosim object (Fiber, Solid, etc.) contains an array points stored as double precision in memory. It is possible to directly access these arrays as numpy arrays without copy operations - at the risk of altering them. In this case an explicit conversion has to be perfomed :   
```python
fiber = simul.fibers[0]  
points = np.array(fiber.data(), copy = False)
```

While this may be useful for the sake of performance, or for efficiently interfacing cytosim with another simulation, it is usually to manipulate a copy of the points, in which case no explicit conversion needs be performed :  
```python
fiber = simul.fibers[0]  
points = fiber.points()
```

Here, Fiber.points() is one of the few accessory functions that are specific to PyCytosim to facilitate the interface to Cytosim objects. Most members and member functions are merely binding to the cytosim members ; for instance, *Simul.fibers* is a binding to the member *fibers* of the cytosim class Simul (e.g. *Simul::fibers*). Therefore the Cytosim doctumentation readily documents the vast majority of PyCytosim.



## Availability

Plyssim is available on [github](https://github.com/SergeDmi/Plyssim), and compiles on Mac and Linux, with Clang, GCC, and ICC (albeit slower with ICC) provided required libraries are installed. A docker image is provided on [dockerhub](https://hub.docker.com/r/sergedmi/plyssim).

# Plyconvert
PlyConvert is a simple command-line utility tool to perform simple manipulations on ply files. It is written in Python 3 and uses the plydata structure from the package plyfile as a backend. Typically it is used to re-scale, align, rename or convert ply files. This can be performed for a single file, or in batch, searching recursively in folders.

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

S.D. aknowledges Nicolas Minc for getting him interested in cell wall mechanics, the entire Minc lab past and present for interactions, and CNRS-Momentum for the funding.

# References 

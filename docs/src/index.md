# Crystals

This package provides a basic framework for dealing with crystalline structures in Julia,
from an atomistic point of view. The main object type, `Crystal`, defines a crystalline
structure in terms of its periodicity, as given by a cell, and atomic sites with different
properties.

Other than the straightforward manipulations of crystalline structures via standard linear
algebra and container API, this package defines a few more complex functions, such as
methods to create supercells and derived primitive cells, or to compute the point-group and
space-group of a crystal.

One of the package's main feature/bug/inconvenience is that it mandates users fully specify
the physical units (say, nanometers, or even kilowatts per seconds if such is your fancy) of
the crystal, thanks to [Unitful](ajkeller34/Unitful.jl). However, as a small give-away to
practical considerations, it does allow crystal sites to be defined and manipulated in terms
of `cartesian` or `fractional` coordinates quite transparently.

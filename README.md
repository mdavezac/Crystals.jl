[![Build Status](https://travis-ci.org/mdavezac/Crystals.jl.svg?branch=master)](https://travis-ci.org/mdavezac/Crystals.jl)

# Crystals

A `Crystal` declares an atomic crystalline structure, e.g an inifinite periodic
arrangement of atoms.

## Usage

The constructor takes at the very least an `n` by `n` array defining the periodicity of the
crystal, i.e. the crystal cell. The cell must have physical units attached to it.

    using Crystals
    crystal = Crystal(eye(3)u"nm")
    @assert crystal.cell === [1 0 0; 0 1 0; 0 0 1]

However, it can also accept atomic positions and any other array of atomic
properties:

    using Crystals
    crystal = Crystal(eye(2)u"km",
                      position=transpose([1 1; 2 3; 4 5])u"m",
                      species=["Al", "O", "O"],
                      label=[:+, :-, :-])
    @assert crystal[:position] == transpose([1 1; 2 3; 4 5])u"m"

Note that the positions form an `n` by `N` array where `N` is the number of atoms. This is
the logical mathematical format if when performing matrix operations on the positions.
However, from an input point of view, it is easier to think one atom at at time, rather than
one coordinate (and all atoms) at a time. Hence the transpose.

Several methods are available to manipulate the crystal cell (`supercell`, `primitive`,
`space_group`, etc..). Check the documentation for more information.

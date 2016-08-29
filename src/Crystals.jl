"""
Crystals
--------

A `Crystal` declares an atomic crystalline structure, e.g an inifinite periodic
arrangement of atoms.

Usage
-----

The constructor takes at the very least an `n` by `n` array  (`2 ≤ n ≤ 6`) and a scale:

    using Crystals
    crystal = Crystal(eye(3), 1)
    @assert crystal.cell === [1 0 0; 0 1 0; 0 0 1]
    @assert crystal.scale === 1

However, it can also accept atomic positions and any other array of atomic
properties:

    crystal = Crystal(eye(2),
                      position=transpose([1 1; 2 3; 4 5]),
                      species=["Al", "O", "O"],
                      label=[:+, :-, :-])
"""
module Crystals
export Position, PositionArray, PositionDataArray, Crystal, volume
export gruber
export hart_forcade, is_periodic, into_cell, origin_centered, into_voronoi,
       supercell
export smith_normal_form
export point_group_operations, inner_translations
export Lattices

module Constants
  const default_tolerance = 1e-8
end

include("Structure.jl")
using .Structure

include("Gruber.jl")
using .Gruber

include("SNF.jl")
using .SNF

include("utilities.jl")
using .Utilities

include("SpaceGroup.jl")
using .SpaceGroup

module Lattices
  using Crystals.Crystal
  include("Bravais.jl")
  include("Binary.jl")
end
end # module

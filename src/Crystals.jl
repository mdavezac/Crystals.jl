"""
Crystals
--------

A `Crystal` declares an atomic crystalline structure, e.g an inifinite periodic
arrangement of atoms.

Usage
-----

The constructor takes at the very least an `n` by `n` array  (`2 ≤ n ≤ 6`) and a scale:

    using Crystals
    crystal = Crystal(eye(3)u"nm", 1)
    @assert crystal.cell === [1 0 0; 0 1 0; 0 0 1]
    @assert crystal.scale === 1

However, it can also accept atomic positions and any other array of atomic
properties:

    crystal = Crystal(eye(2)u"nm",
                      position=transpose([1 1; 2 3; 4 5]),
                      species=["Al", "O", "O"],
                      label=[:+, :-, :-])
"""
module Crystals
using Unitful: @u_str
export @u_str

# export Position, PositionArray, PositionDataArray, is_fractional
export Crystal, is_fractional, volume, round!
export eachatom
export smith_normal_form
export gruber, niggly
export hart_forcade, is_periodic, to_fractional, to_cartesian, into_cell, origin_centered
export into_voronoi, supercell, cell_parameters
export point_group_operations, inner_translations, is_primitive
# , primitive,
#        space_group
export Lattices

include("Logging.jl")

module Constants
  const default_tolerance = 1e-8
end

# include("Positions.jl")
# using .Positions

include("Structures.jl")
using .Structures

include("CrystalAtoms.jl")
using .CrystalAtoms

include("SNF.jl")
using .SNF

include("utilities.jl")
using .Utilities

include("Gruber.jl")
using .Gruber

include("SpaceGroup.jl")
using .SpaceGroup

module Lattices
  using Crystals.Structures.Crystal
  using Unitful: @u_str
  include("Bravais.jl")
  include("Binary.jl")
  include("A2BX4.jl")
end
end # module

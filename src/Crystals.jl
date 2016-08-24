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
                      label=[:+, :-, :-]
"""
module Crystals

using DataFrames: DataFrame, DataArray, ColumnIndex, index, NA
import DataFrames: deleterows!, hcat!, nullable!, pool!, ourshowcompact, nrow
using FixedSizeArrays: FixedVectorNoTuple

include("Structure.jl")
# include("SpaceGroup.jl")
include("Gruber.jl")
import .Gruber: gruber
const default_tolerance = Gruber.default_tolerance

include("SNF.jl")
import .SNF: smith_normal_form

include("utilities.jl")

export Crystal, PositionArray, PositionDataArray, deleterows!,
      nullable!, Position, gruber, smith_normal_form, into_cell, into_voronoi,
      is_periodic, origin_centered, supercell

end # module

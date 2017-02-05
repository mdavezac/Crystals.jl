Atomic positions can be initialized, modified, and accessed in cartesian or fractional
coordinates. Cartesian coordinates refer to *real world* positions in the same physical
units as the crystal cell. Fractional coordinates however are in units *of* the crystal
cell.

## Creating and accessing fractional and cartesian coordinates

In the following, we create a crystal structure using fractional coordinates, through
the simple expedience of *not* specifying actual units.

```@meta
CurrentModule = Crystals
DocTestSetup = quote
    using Crystals
end
```
```jldocset
frac_crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u"nm",
                       position=[0, 0, 0],
                       position=[0.25, 0.25, 0.25])
@assert frac_crystal[:position] === frac_crystal[:fractional]
@assert frac_crystal[:cartesian] ≈ frac_crystal.cell * frac_crystal[:fractional]
display(frac_crystal[:cartesian])

# output
3×2 Array{Quantity{Float64, Dimensions:{𝐋}, Units:{nm}},2}:
 0.0 nm  1.05 nm
 0.0 nm  1.05 nm
 0.0 nm  1.05 nm
```

Note that querying `:position` returns fractional coordinates. If we create a structure with
cartesian coordinates instead -- by calling the constructor with positions that have units
-- then querying `:position` would return the cartesian coordinates.

```jldocset
cart_crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u"nm",
                       position=[0, 0, 0]u"nm",
                       position=[1.05, 1.05, 1.05]u"nm")
@assert cart_crystal[:position] === cart_crystal[:fractional]
@assert cart_crystal[:fractional] ≈ inv(cart_crystal.cell) * cart_crystal[:cartesian]
display(cart_crystal[:fractional])

# output
3×2 Array{Float64,2}:
  0.0  1.05
  0.0  1.05
  0.0  1.05
```

Of course, in either case, we can access either `:cartesian` or `:fractional` coordinates
through the relevant column name. However, depending on how the crystal was created, one of
these calls will be essentially a no-op, and the other will involve a matrix-matrix
multiplication (and possibly computing the inverse of a matrix). To obtain a specific kind
of crystal from any other crystal, one can simply use the bracket operator:


```jldocset
crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u"nm",
                  position=[0, 0, 0]u"nm",
                  position=[1.05, 1.05, 1.05]u"nm",
                  species=["Si", "Si"])
crystal[[:fractional, :species]]

# output
cell(nm):
  0.0 2.1 2.1
  2.1 0.0 2.1
  2.1 2.1 0.0
│ Atom │ fractional       │ species │
├──────┼──────────────────┼─────────┤
│ 1    │ (0.0,0.0,0.0)    │ "Si"    │
│ 2    │ (0.25,0.25,0.25) │ "Si"    │
```

Note that the column name explicitly specifies `fractional`, as opposed to `cartesian`.
This call (indeed, all bracket operator call) will always create a new instance of a
`Crystal`, whether it is strictly needed or not.

## Setting coordinates from fractional or cartesian inputs

When setting a position, it is not necessary to specify whether the input is fractional or
cartesian. If the input has physical units, then it is cartesian. If it doesn't, then it is
fractional. The appropriate transformation is applied before setting the position in the
crystal.

```jldocset
frac_crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u"nm",
                       position=[0, 0, 0],
                       position=[0.25, 0.25, 0.25])
frac_crystal[2, :position] = [1, 1, 1]u"nm"
frac_crystal

# output
cell(nm):
  0.0 2.1 2.1
  2.1 0.0 2.1
  2.1 2.1 0.0
│ Atom │ fractional                   │
├──────┼──────────────────────────────┤
│ 1    │ (0.0,0.0,0.0)                │
│ 2    │ (0.238095,0.238095,0.238095) │
```


## Misc

The function `is_fractional` can tell us programmatically whether a structure is fractional
or cartesian.

```@docs
is_fractional
```

Apart from `:cartesian` and `:fractional`, there are three other special column names that
ease access to specific atomic properties. `:x`, `:y`, `:z` will return an array
representing the corresponding coordinate, in the same system -- cartesian or fractional --
as the crystal (and as `:position`). `:x`, `:y`, `:z`  can be used to set a specific
coordinate as well. Only three such special names are provided. 

```jldoctest
crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u"nm",
                  position=[0, 0, 0],
                  position=[0.25, 0.27, 0.25])

println("Y coordinate: ", crystal[2, :y])
crystal[2, :y] = 0.25
crystal

# output
Y coordinate: 0.27
cell(nm):
  0.0 2.1 2.1
  2.1 0.0 2.1
  2.1 2.1 0.0
│ Atom │ fractional       │
├──────┼──────────────────┤
│ 1    │ (0.0,0.0,0.0)    │
│ 2    │ (0.25,0.25,0.25) │
```

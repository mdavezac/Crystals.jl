## Simple Construction

A `Crystal` declares an atomic crystalline structure, e.g. an infinite periodic arrangement
of atoms. The constructor takes at the very least an `n` by `n` array defining
the periodicity of the crystal, i.e. the crystal cell. The cell must have physical units
attached to it.

```jldoctest
using Crystals
crystal = Crystal(eye(3)u"nm")

# output
cell(nm):
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
```

However, it can also accept atomic positions and any other array of atomic
properties:

```@meta
DocTestSetup = quote
  using Crystals
end
```
```jldoctest
crystal = Crystal(eye(2)u"km",
                  position=transpose([1 1; 2 3; 4 5])u"m",
                  species=["Al", "O", "O"],
                  label=[:+, :-, :-])

# output
cell(m):
  1000.0 0.0
  0.0 1000.0
│ Atom │ Cartesian │ species │ label │
├──────┼───────────┼─────────┼───────┤
│ 1    │ (1.0,1.0) │ "Al"    │ +     │
│ 2    │ (2.0,3.0) │ "O"     │ -     │
│ 3    │ (4.0,5.0) │ "O"     │ -     │
```

Note that the positions form an `n` by `N` array where `N` is the number of atoms. This is
the logical mathematical format if when performing matrix operations on the positions.
However, from an input point of view, it is easier to think one atom at at time, rather than
one coordinate (and all atoms) at a time. Hence the transpose.



```@eval
using Base.Markdown
using Crystals
result = """
!!! note

    The following are reserved column names that have special meaning. They can be used in
    most, but not all circumstances. For instance, `:x`, `:y`, and `:x` cannot be used in
    the constructor.
"""
result *= "    - :" * join(map(string, Crystals.Structures.RESERVED_COLUMNS), "\n    - :")
Base.Markdown.parse(result)
```

## Accessing the cell and atomic sites

Access to the crystal cell happens via `.` call, `crystal.cell`. Atomic properties on the
other hand are accessed and modified through the square bracket operator. There are several
ways of doing this, more or less reflecting what can be done with a
[DataFrame](www.github.com/JuliaStags/DataFrames.jl):

```@meta
DocTestSetup = quote
  using Crystals
  crystal = Crystal(eye(2)u"km",
                  position=transpose([1 1; 2 3; 4 5])u"m",
                  species=["Al", "O", "O"],
                  label=[:+, :-, :-])
end
```

```jldoctest
println(crystal[:label])
println(crystal[[1, 3], [:species, :label]])

crystal[1, :position] = [0, 4]u"m"

# output
Symbol[:+,:-,:-]
2×2 DataFrames.DataFrame
│ Row │ species │ label │
├─────┼─────────┼───────┤
│ 1   │ "Al"    │ +     │
│ 2   │ "O"     │ -     │
2-element Array{Quantity{Int64, Dimensions:{𝐋}, Units:{m}},1}:
 0 m
 4 m
```

However, the main access pattern is atom/row-based. Using an integer or a sequence of
integers will create a new `Crystal` structure with only the selected atoms.

```jldoctest
crystal[1]

# output
cell(m):
  1000.0 0.0
  0.0 1000.0
│ Atom │ Cartesian │ species │ label │
├──────┼───────────┼─────────┼───────┤
│ 1    │ (1.0,1.0) │ "Al"    │ +     │
```

Note that the return is still a crystalline structure. In a way, we are selecting an atom
*and* it's periodic image, rather than single atom.


## Appending and removing atoms

To add atoms to a `Crystal` instance, it is recommended to use the `push!`, `append!`, and
`vcat` functions. Removing atoms and atomic properties from a structure can be done quite
easily with `delete!`. Or alternatively, one can select a few atoms and properties using the
bracket notation, e.g. `crystal[[1, 2], [:species, :label]]`.

```@docs
push!
append!
vcat
delete!
DataFrames.deleterows!
empty!
```

!!! warning

    Manipulating directly the fields `positions` and `properties` of a `Crystal` instance is
    discouraged. In practice, the only requirement is that the two represent the same number
    of atoms (with the exception of crystals with no atomic properties outside of the
    positions, for which the `properties` field is empty). However, no provision is made
    anywhere in the code for abnormal `Crystal` structure. This means the code may crash in
    creative ways.


## Iterating through a structure

It is a possible to iterate through all the atoms of a structure. The object yielded
references directly the parent structure. Atomic properties and positions can be accessed
and modified one site at a time.

```jldoctest
for atom in eachatom(crystal)
    println("Atom: ", atom[:species], " and ", atom[:label])
    atom[:species] *= string(atom[:label])
end
crystal

# output
Atom: Al and +
Atom: O and -
Atom: O and -
cell(m):
  1000.0 0.0
  0.0 1000.0
│ Atom │ Cartesian │ species │ label │
├──────┼───────────┼─────────┼───────┤
│ 1    │ (1.0,1.0) │ "Al+"   │ +     │
│ 2    │ (2.0,3.0) │ "O-"    │ -     │
│ 3    │ (4.0,5.0) │ "O-"    │ -     │
```

!!! warning

    Modifying the number of sites in the crystal structure during iteration will result in
    undefined behaviour.

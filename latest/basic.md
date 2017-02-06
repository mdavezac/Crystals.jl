
<a id='Simple-Construction-1'></a>

## Simple Construction


A `Crystal` declares an atomic crystalline structure, e.g. an infinite periodic arrangement of atoms. The constructor takes at the very least an `n` by `n` array defining the periodicity of the crystal, i.e. the crystal cell. The cell must have physical units attached to it.


```julia
using Crystals
crystal = Crystal(eye(3)u"nm")

# output
cell(nm):
  1.0 0.0 0.0
  0.0 1.0 0.0
  0.0 0.0 1.0
```


However, it can also accept atomic positions and any other array of atomic properties:




```julia
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


Note that the positions form an `n` by `N` array where `N` is the number of atoms. This is the logical mathematical format if when performing matrix operations on the positions. However, from an input point of view, it is easier to think one atom at at time, rather than one coordinate (and all atoms) at a time. Hence the transpose.


!!! note
    The following are reserved column names that have special meaning. They can be used in most, but not all circumstances. For instance, `:x`, `:y`, and `:x` cannot be used in the constructor.

      * :position
      * :fractional
      * :cartesian
      * :x
      * :y
      * :z



<a id='Accessing-the-cell-and-atomic-sites-1'></a>

## Accessing the cell and atomic sites


Access to the crystal cell happens via `.` call, `crystal.cell`. Atomic properties on the other hand are accessed and modified through the square bracket operator. There are several ways of doing this, more or less reflecting what can be done with a [DataFrame](www.github.com/JuliaStags/DataFrames.jl):




```julia
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


However, the main access pattern is atom/row-based. Using an integer or a sequence of integers will create a new `Crystal` structure with only the selected atoms.


```julia
crystal[1]

# output
cell(m):
  1000.0 0.0
  0.0 1000.0
│ Atom │ Cartesian │ species │ label │
├──────┼───────────┼─────────┼───────┤
│ 1    │ (1.0,1.0) │ "Al"    │ +     │
```


Note that the return is still a crystalline structure. In a way, we are selecting an atom *and* it's periodic image, rather than single atom.


<a id='Appending-and-removing-atoms-1'></a>

## Appending and removing atoms


To add atoms to a `Crystal` instance, it is recommended to use the `push!`, `append!`, and `vcat` functions. Removing atoms and atomic properties from a structure can be done quite easily with `delete!`. Or alternatively, one can select a few atoms and properties using the bracket notation, e.g. `crystal[[1, 2], [:species, :label]]`.

<a id='Base.push!' href='#Base.push!'>#</a>
**`Base.push!`** &mdash; *Function*.



```julia
push!(crystal, position; kwargs...)

```

Appends an atom to a crystal structure. The position of the atom is a necessary argument, whether in Cartesian or in fractional coordinates. If keyword arguments are present, then they represent atomic properties for the atom being added. Properties that are not explicitly given are set to `NA`. Similarly, new properties that were not present in the crystal structure previously are `NA` except for the newly added atom.

**Examples**

```jldocset
push!(crystal, [10, 20]u"nm", species="B", μ=1)
crystal

# output
cell(m):
1000.0 0.0
0.0 1000.0
│ Atom │ position        │ species │ label │ μ  │
├──────┼─────────────────┼─────────┼───────┼────┤
│ 1    │ (1.0,1.0)       │ "Al"    │ +     │ NA │
│ 2    │ (2.0,3.0)       │ "O"     │ -     │ NA │
│ 3    │ (4.0,5.0)       │ "O"     │ -     │ NA │
│ 4    │ (1.0e-8,2.0e-8) │ "B"     │ NA    │ 1  │
```


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L231-L257' class='documenter-source'>source</a><br>

<a id='Base.append!' href='#Base.append!'>#</a>
**`Base.append!`** &mdash; *Function*.



```julia
append!(crystal, other; check_periodicity)

```

Appends one or more crystal structure to the first structure. Unless `check_periodicity` is `false`, the structures must have the exact same periodicity. An error will result otherwise.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L617-L623' class='documenter-source'>source</a><br>

<a id='Base.vcat' href='#Base.vcat'>#</a>
**`Base.vcat`** &mdash; *Function*.



```julia
vcat(crystal, other; check_periodicity)

```

Concatenates crystals together. The lattices must be compatible.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L601-L605' class='documenter-source'>source</a><br>

<a id='Base.delete!' href='#Base.delete!'>#</a>
**`Base.delete!`** &mdash; *Function*.



```
Base.delete!(crystal::Crystal, col::Symbol)
Base.delete!(crystal::Crystal, col::AbstractVector{Symbol})
```

Deletes one or more atomic property. Positions cannot be deleted.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L537-L542' class='documenter-source'>source</a><br>


```
Base.delete!(crystal::Crystal, ::Colon)
```

Deletes all atomic properties except for positions.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L551-L555' class='documenter-source'>source</a><br>


```
delete!(crystal::Crystal, rows::Integer)
delete!(crystal::Crystal, rows::AbstractVector{Integer})
delete!{T <: Integer}(crystal::Crystal, rows::Range{T})
delete!(crystal::Crystal, rows::Colon)
```

Alias for `deleterows`[@ref].


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L591-L598' class='documenter-source'>source</a><br>

<a id='DataFrames.deleterows!' href='#DataFrames.deleterows!'>#</a>
**`DataFrames.deleterows!`** &mdash; *Function*.



```
deleterows!(crystal::Crystal, rows::Integer)
deleterows!(crystal::Crystal, rows::AbstractVector{Integer})
deleterows!{T <: Integer}(crystal::Crystal, rows::Range{T})
deleterows!(crystal::Crystal, rows::Colon)
```

Deletes one (single integer), a few (sequence or range), or all (colon) atoms in the structure.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L570-L578' class='documenter-source'>source</a><br>

<a id='Base.empty!' href='#Base.empty!'>#</a>
**`Base.empty!`** &mdash; *Function*.



```julia
empty!(crystal)

```

Deletes all atomic sites, both properties and positions.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L560-L564' class='documenter-source'>source</a><br>


!!! warning
    Manipulating directly the fields `positions` and `properties` of a `Crystal` instance is discouraged. In practice, the only requirement is that the two represent the same number of atoms (with the exception of crystals with no atomic properties outside of the positions, for which the `properties` field is empty). However, no provision is made anywhere in the code for abnormal `Crystal` structure. This means the code may crash in creative ways.



<a id='Iterating-through-a-structure-1'></a>

## Iterating through a structure


It is a possible to iterate through all the atoms of a structure. The object yielded references directly the parent structure. Atomic properties and positions can be accessed and modified one site at a time.


```julia
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
    Modifying the number of sites in the crystal structure during iteration will result in undefined behaviour.



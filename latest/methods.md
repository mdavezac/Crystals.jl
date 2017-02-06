
<a id='Querying-and-property-methods-1'></a>

## Querying and property methods



<a id='Crystals.Structures.is_fractional' href='#Crystals.Structures.is_fractional'>#</a>
**`Crystals.Structures.is_fractional`** &mdash; *Function*.



```julia
is_fractional(crystal)

```

True if the crystal structure is fractional.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L90-L94' class='documenter-source'>source</a><br>

<a id='Crystals.Structures.are_compatible_lattices' href='#Crystals.Structures.are_compatible_lattices'>#</a>
**`Crystals.Structures.are_compatible_lattices`** &mdash; *Function*.



```
are_compatible_lattices(lattices...)
```

True if the lattices are mathematically equivalent. Two lattices are equivalent if they represent the same periodicity. In practice, this means the two lattices have the same volume, and their cell vectors are integer linear combinations of one another.

**Parameters**

  * `lattices::Vararg{Union{Matrix, Crystal}}`: Any number of lattices
  * `digits::Integer`: when checking the cells are integer co-linear, the product $A^{-1}B$ is first rounded to this number of digits
  * `rtol::Real`: relative tolerance when checking the volumes correspond (no units). Default to 1.0e-8.
  * `atol::Real`: absolute tolerance when checking the volumes correspond (no units). Default to 1.0e-8.

**Examples**

```julia
using Crystals, Unitful
crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u"nm")
cells = Matrix{Int64}[eye(3)]
while length(cells) < 5
    cell = rand(-5:5, (3, 3))
    volume(cell) == 1 && push!(cells, cell)
end
lattices = [crystal.cell * c for c in cells]
are_compatible_lattices(crystal, lattices...)

# output
true
```


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L143-L176' class='documenter-source'>source</a><br>

<a id='Crystals.Structures.volume' href='#Crystals.Structures.volume'>#</a>
**`Crystals.Structures.volume`** &mdash; *Function*.



```julia
volume(crystal)

```

Returns the volume of a `Crystal` instance or of a cell. It comes down to computing $|det(A)|$ where $A$ is the crystal cell.

**Examples**

```jlcon
julia> using Crystals

julia> volume(Crystal([0 2 2; 2 0 2; 2 2 0]u"nm"))
16.0 nm^3

julia> volume([0 2 2; 2 0 2; 2 2 0]u"nm")
16.0 nm^3

julia> volume([0 2 2; 2 0 2; 2 2 0])
16.0
```


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L106-L126' class='documenter-source'>source</a><br>

<a id='Crystals.Utilities.cell_parameters' href='#Crystals.Utilities.cell_parameters'>#</a>
**`Crystals.Utilities.cell_parameters`** &mdash; *Function*.



```
cell_parameters(a::Quantity, b::Quantity, c::Quantity,
                α::Quantity=(π/2)u"rad", β::Quantity=(π/2)u"rad",
                γ::Quantity=(π/2)u"rad")
```

Computes the cell matrix from the cell parameters [a, b, c, α, β, γ].


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/utilities.jl#L301-L307' class='documenter-source'>source</a><br>


```
cell_parameters(cell::AbstractMatrix)
cell_parameters(cell::Crystal)
```

Parameters (a, b, c, α, β, γ) of the input cell returned in a named tuple.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/utilities.jl#L323-L328' class='documenter-source'>source</a><br>

<a id='Crystals.Gruber.gruber' href='#Crystals.Gruber.gruber'>#</a>
**`Crystals.Gruber.gruber`** &mdash; *Function*.



```
gruber(cell::Matrix;
       tolerance::Real=default_tolerance, itermax::Unsigned=50,
       max_no_change::Unsigned=10)
```

Determines Gruber cell of an input cell.

The Gruber cell is an optimal parameterization of a lattice, e.g. shortest cell-vectors and angles closest to 90 degrees. The new cell is in the same basis as the origin cell: no rotation has been incurred. The cell parameters are uniquely determined, even though the cell itself is not (certain symmetry operations leaving the structure unchanged may yield a more recognizable cell). If you want a unique Cartesian cell (in a different Cartesian basis), use the `niggly` algorithm.

**Arguments**

  * `cell::Matrix`: the input lattice cell-vectors. Cannot be singular.
  * `itermax::Integer`: maximum number of iterations before bailing out
  * `tolerance::Number`: tolerance parameter when comparing real numbers
  * `max_no_change::Integer`: Maximum number of times to go through algorithm without changes before bailing out


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Gruber.jl#L110-L131' class='documenter-source'>source</a><br>

<a id='Crystals.Gruber.niggly' href='#Crystals.Gruber.niggly'>#</a>
**`Crystals.Gruber.niggly`** &mdash; *Function*.



```
niggly(cell::Matrix; kwargs...)
```

Determines a unique Cartesian cell equivalent to the input, with the shortest possible vectors and squarest angles. For an explanation of the parameters, see `gruber`. In practice, this function computes the cell-parameters of a `gruber` cell and then reconstructs the cell matrix. Hence, the result should be quite unique for any lattice representation, including any rotation of the underlying Cartesian coordinates.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Gruber.jl#L204-L212' class='documenter-source'>source</a><br>

<a id='Crystals.Utilities.is_periodic' href='#Crystals.Utilities.is_periodic'>#</a>
**`Crystals.Utilities.is_periodic`** &mdash; *Function*.



```
is_periodic(a::AbstractVector, b::AbstractVector, cell::AbstractMatrix;
            tolerance::Real=1.0e-8)
```

True if the positions are one-to-one periodic with respect to the input cell.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/utilities.jl#L158-L163' class='documenter-source'>source</a><br>


```
is_periodic(a::AbstractMatrix, b::AbstractVector, cell::AbstractMatrix;
            tolerance::Real=1.0e-8)
```

Array of boolean describing whether positions in `a` are periodic with positions in `b`.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/utilities.jl#L172-L177' class='documenter-source'>source</a><br>

<a id='Crystals.SpaceGroup.is_primitive' href='#Crystals.SpaceGroup.is_primitive'>#</a>
**`Crystals.SpaceGroup.is_primitive`** &mdash; *Function*.



```
is_primitive(cartesian::AbstractMatrix, cell::AbstractMatrix, species::AbstractVector;
             tolerance::Real=1.0e-8)
is_primitive(crystal::Crystal, col::Union{Symbol, AbstractVector{Symbol}}; kwargs...)
is_primitive(crystal::Crystal; kwargs...)
```

True if the crystal structure is primitive, e.g. not a supercell, e.g. not reducible to an equivalent lattice with fewer atoms.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/SpaceGroup.jl#L190-L198' class='documenter-source'>source</a><br>


Typical container and DataFrame methods are also available, such as `names`, `size`, `ndims`, `nrow`, `ncol`, `endof`.


<a id='Looping-1'></a>

## Looping

<a id='Crystals.CrystalAtoms.eachatom' href='#Crystals.CrystalAtoms.eachatom'>#</a>
**`Crystals.CrystalAtoms.eachatom`** &mdash; *Function*.



Iterator over each atom in the crystal 


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/CrystalAtoms.jl#L30' class='documenter-source'>source</a><br>


`eachindex` is also available. It returns the range over atoms indices `1:size(crystal, 1)`.


<a id='Modifying,-building-up-and-building-down-1'></a>

## Modifying, building up and building down

<a id='Base.round' href='#Base.round'>#</a>
**`Base.round`** &mdash; *Function*.



```julia
round(crystal, args)

```

Rounds the cell and positions of a crystal. See `Base.round` for possible parameters.

**Examples**

```julia
using Crystals
crystal = Crystal([0 0.501 0.501; 0.496 0.001 0.497; 0.497 0.497 0]u"nm",
                  position=[0.001, -0.001, -0.001]u"nm",
                  position=[0.25, 0.251, -0.247]u"nm")
round(crystal, 2)

# output

cell(nm):
  0.0 0.5 0.5
  0.5 0.0 0.5
  0.5 0.5 0.0
│ Atom │ Cartesian         │
├──────┼───────────────────┤
│ 1    │ (0.0,-0.0,-0.0)   │
│ 2    │ (0.25,0.25,-0.25) │
```


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L653-L678' class='documenter-source'>source</a><br>

<a id='Crystals.Structures.round!' href='#Crystals.Structures.round!'>#</a>
**`Crystals.Structures.round!`** &mdash; *Function*.



```julia
round!(crystal, args)

```

Rounds the cell and positions of a crystal. See `round` for possible parameters.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Structures.jl#L634-L638' class='documenter-source'>source</a><br>

<a id='Crystals.Utilities.into_cell' href='#Crystals.Utilities.into_cell'>#</a>
**`Crystals.Utilities.into_cell`** &mdash; *Function*.



```julia
into_cell(pos, cell; tolerance)

```

Folds periodic positions into cell


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/utilities.jl#L203-L207' class='documenter-source'>source</a><br>

<a id='Crystals.Utilities.into_voronoi' href='#Crystals.Utilities.into_voronoi'>#</a>
**`Crystals.Utilities.into_voronoi`** &mdash; *Function*.



```
into_voronoi(positions::AbstractArray, cell::AbstractMatrix; extent::Integer=1)
```

Folds positions into first Brillouin zone of the input cell. Makes a well-meaning effort at returning the periodic image with the smallest possible norm. It recenter the atoms around the origin and then looks for the smallest periodic images within `-extent:extent` cells. If the cell is quite pathological, then the result will not be within the Voronoi cell.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/utilities.jl#L233-L240' class='documenter-source'>source</a><br>

<a id='Crystals.Utilities.origin_centered' href='#Crystals.Utilities.origin_centered'>#</a>
**`Crystals.Utilities.origin_centered`** &mdash; *Function*.



```
origin_centered(positions::AbstractArrays, cell::AbstractMatrix)
```

Folds positions back to origin, such that each fractional component $x_f$ is between $-0.5\leq x_f < 0.5$. If the input is in Cartesian (fractional) coordinates, then the result is also in Cartesian (fractional) coordinates.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/utilities.jl#L222-L228' class='documenter-source'>source</a><br>


[`is_primitive`](methods.md#Crystals.SpaceGroup.is_primitive) can be queried to figure out whether a crystal is a primitive cell or a supercell.


<a id='Point-group,-space-group,-and-geometry-1'></a>

## Point-group, space-group, and geometry

<a id='Crystals.SpaceGroup.point_group' href='#Crystals.SpaceGroup.point_group'>#</a>
**`Crystals.SpaceGroup.point_group`** &mdash; *Function*.



Finds and stores point group operations for a given lattice

A lattice is defined by a 3x3 matrix or cell.  Rotations are determined from G-vector triplets with the same norm as the unit-cell vectors.

Implementation taken from [ENUM](http://enum.sourceforge.net/).


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/SpaceGroup.jl#L47-L54' class='documenter-source'>source</a><br>

<a id='Crystals.SpaceGroup.space_group' href='#Crystals.SpaceGroup.space_group'>#</a>
**`Crystals.SpaceGroup.space_group`** &mdash; *Function*.



```julia
space_group(crystal; kwargs...)

```

Computes the space-group operations of a crystal. By default, all atomic properties are considered when determining whether atomic sites are related by symmetry. However, it is possible to specify a subset of atomic properties. An empty subset of atomic properties implies that all atoms sites are equivalent. It is also possible to specify an array of integers, serving as labels for each atomic site.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/SpaceGroup.jl#L353-L361' class='documenter-source'>source</a><br>


```julia
space_group(cell, positions, species; kwargs...)

```

```
Computes the space-group of a crystal specified using standard Julia types.
```


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/SpaceGroup.jl#L376-L380' class='documenter-source'>source</a><br>

<a id='Crystals.Utilities.hart_forcade' href='#Crystals.Utilities.hart_forcade'>#</a>
**`Crystals.Utilities.hart_forcade`** &mdash; *Function*.



```julia
hart_forcade(lattice, supercell; digits)

```

Computes the cyclic group of a supercell with respect to a lattice. It makes it possible to identify the class of periodically equivalent cell that a given position within the supercell belongs to. The function returns a named tuple with the transform and the quotient.

**Examples**

```jldocset
using Crystals, Unitful
fcc = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"
supercell = [0 2 2; 0 -4 2; -1 0 -2]
ht = hart_forcade(fcc, fcc * supercell)
display(ht)

println("Positions in supercell:")
for index in CartesianRange((ht.quotient...))
    position = inv(ht.transform) * [index[u] for u in eachindex(ht.quotient)]
    println("- ", ustrip(position), " (", unit(eltype(position)), ")")
end

# output

Hart-Forcade transform
- transform (nm^-1): [-1.0 -1.0 1.0; -1.0 1.0 1.0; -1.0 1.0 3.0]
- quotient: [1,2,6]

Positions in supercell:
- [-1.0,0.0,0.0] (nm)
- [-2.0,0.5,-0.5] (nm)
- [-0.5,0.0,0.5] (nm)
- [-1.5,0.5,0.0] (nm)
- [0.0,0.0,1.0] (nm)
- [-1.0,0.5,0.5] (nm)
- [0.5,0.0,1.5] (nm)
- [-0.5,0.5,1.0] (nm)
- [1.0,0.0,2.0] (nm)
- [0.0,0.5,1.5] (nm)
- [1.5,0.0,2.5] (nm)
- [0.5,0.5,2.0] (nm)
```


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/utilities.jl#L34' class='documenter-source'>source</a><br>

<a id='Crystals.SpaceGroup.primitive' href='#Crystals.SpaceGroup.primitive'>#</a>
**`Crystals.SpaceGroup.primitive`** &mdash; *Function*.



```
primitive(crystal::Crystal; tolerance::Real=1.0e-8)
```

Computes the primitive cell of the input crystal. If the crystal is primitive, it is returned as is, otherwise a new crystal is returned.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/SpaceGroup.jl#L210-L215' class='documenter-source'>source</a><br>

<a id='Crystals.Utilities.supercell' href='#Crystals.Utilities.supercell'>#</a>
**`Crystals.Utilities.supercell`** &mdash; *Function*.



```julia
supercell(lattice, supercell; site_id, tolerance)

```

Creates a supercell from an input lattice.

# Parameters

  * `lattice::Crystal`: the original lattice
  * `supercell::AbstractMatrix`: the cell of the supercell in Cartesian coordinates
  * `site_id::Bool`: Whether to add/modify an atomic property indicating the index of the site in the original lattice
  * `cell_id::Bool`: Whether to add/modify an atomic property indicating the index of the cell the site belongs to


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/utilities.jl#L259-L272' class='documenter-source'>source</a><br>


<a id='Predefined-lattices-1'></a>

## Predefined lattices


A small number of standard lattices are available for construction in the exported `Lattices` submodule.


For instance, the body-centered lattice:


```julia
using Crystals
Lattices.bcc()
```

```
cell(nm):
  -0.5 0.5 0.5
  0.5 -0.5 0.5
  0.5 0.5 -0.5
│ Atom │ Cartesian     │
├──────┼───────────────┤
│ 1    │ (0.0,0.0,0.0) │
```


Or more complex lattices, such as b5 spinels:


```julia
Lattices.b5()
```

```
cell(nm):
  0.0 0.5 0.5
  0.5 0.0 0.5
  0.5 0.5 0.0
│ Atom │ Cartesian           │ species │
├──────┼─────────────────────┼─────────┤
│ 1    │ (0.5,0.5,0.5)       │ 'A'     │
│ 2    │ (0.5,0.25,0.25)     │ 'A'     │
│ 3    │ (0.25,0.5,0.25)     │ 'A'     │
│ 4    │ (0.25,0.25,0.5)     │ 'A'     │
│ 5    │ (0.875,0.875,0.875) │ 'B'     │
│ 6    │ (0.125,0.125,0.125) │ 'B'     │
│ 7    │ (0.25,0.25,0.25)    │ 'X'     │
│ 8    │ (0.25,0.5,0.5)      │ 'X'     │
│ 9    │ (0.5,0.25,0.5)      │ 'X'     │
│ 10   │ (0.5,0.5,0.25)      │ 'X'     │
│ 11   │ (0.75,0.75,0.75)    │ 'X'     │
│ 12   │ (0.75,0.5,0.5)      │ 'X'     │
│ 13   │ (0.5,0.75,0.5)      │ 'X'     │
│ 14   │ (0.5,0.5,0.75)      │ 'X'     │
```

<a id='Crystals.Lattices.b5' href='#Crystals.Lattices.b5'>#</a>
**`Crystals.Lattices.b5`** &mdash; *Function*.



```julia
b5(T; unit)
b5()

```

b5 (spinel) lattice 


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/A2BX4.jl#L1-L5' class='documenter-source'>source</a><br>

<a id='Crystals.Lattices.bcc' href='#Crystals.Lattices.bcc'>#</a>
**`Crystals.Lattices.bcc`** &mdash; *Function*.



```julia
bcc()
bcc(T; unit)

```

Creates a body-centered lattice with a single site 


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Bravais.jl#L1-L5' class='documenter-source'>source</a><br>

<a id='Crystals.Lattices.diamond' href='#Crystals.Lattices.diamond'>#</a>
**`Crystals.Lattices.diamond`** &mdash; *Function*.



```julia
diamond()
diamond(T; unit)

```

Diamond lattice 


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Binary.jl#L16-L20' class='documenter-source'>source</a><br>

<a id='Crystals.Lattices.fcc' href='#Crystals.Lattices.fcc'>#</a>
**`Crystals.Lattices.fcc`** &mdash; *Function*.



```julia
fcc(T; unit)
fcc()

```

Creates an face centered lattice with a single site 


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Bravais.jl#L7-L11' class='documenter-source'>source</a><br>

<a id='Crystals.Lattices.rock_salt' href='#Crystals.Lattices.rock_salt'>#</a>
**`Crystals.Lattices.rock_salt`** &mdash; *Function*.



```julia
rock_salt()
rock_salt(T; unit)

```

Rock-salt lattice 


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Binary.jl#L1-L5' class='documenter-source'>source</a><br>

<a id='Crystals.Lattices.wurtzite' href='#Crystals.Lattices.wurtzite'>#</a>
**`Crystals.Lattices.wurtzite`** &mdash; *Function*.



```julia
wurtzite()
wurtzite(T; unit)

```

Wurtzite lattice 


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Binary.jl#L22-L26' class='documenter-source'>source</a><br>

<a id='Crystals.Lattices.zinc_blende' href='#Crystals.Lattices.zinc_blende'>#</a>
**`Crystals.Lattices.zinc_blende`** &mdash; *Function*.



```julia
zinc_blende()
zinc_blende(T; unit)

```

Zinc-blende lattice 


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Binary.jl#L9-L13' class='documenter-source'>source</a><br>


<a id='Output-and-Logging-1'></a>

## Output and Logging


Information and errors during calculations are displayed using an internal log provided by [Lumberjack](https://www.github.com/WestleyArgentum/Lumberjack.jl). By default, only critical errors result in output. The verbosity can be set manually.

<a id='Crystals.Log.set_log_level' href='#Crystals.Log.set_log_level'>#</a>
**`Crystals.Log.set_log_level`** &mdash; *Function*.



```julia
set_log_level()
set_log_level(level)

```

Modifies log-level of all "trucks" in Crystals logs. The input should be one of "debug", "info", "warn", "error", from least to most verbose.


<a target='_blank' href='https://github.com/mdavezac/Crystals.jl/tree/7f715d66462d7db5977704035cf74211abd39dac/src/Logging.jl#L32-L37' class='documenter-source'>source</a><br>


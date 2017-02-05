## Querying and property methods

```@meta
CurrentModule = Crystals
```
```@docs
is_fractional
are_compatible_lattices
volume
cell_parameters
gruber
niggly
is_periodic
is_primitive
```

Typical container and dataframe methods are also available, such as `names`, `size`,
`ndims`, `nrow`, `ncol`, `endof`.

## Looping

```@docs
eachatom
```

`eachindex` is also available. It returns the range over atoms indices `1:size(crystal, 1)`.

## Modifying, building up and building down

```@docs
round
round!
into_cell
into_voronoi
origin_centered
```

[`is_primitive`](@ref) can be queried to figure out whether a crystal is a primitive cell or
a supercell.

## Point-group, space-group, and geometry

```@docs
point_group
space_group
hart_forcade
primitive
supercell
```

## Predefined lattices

A small number of standard lattices are available for construction in the exported
`Lattices` submodule.

For instance, the body-centered lattice:

```@example 1
using Crystals
Lattices.bcc()
```

Or more complex lattices, such as b5 spinels:

```@example 1
Lattices.b5()
```


```@autodocs
Modules = [Crystals.Lattices]
```


## Output and Logging

Information and errors during calculations are displayed using an internal log provided by
[Lumberjack](https://www.github.com/WestleyArgentum/Lumberjack.jl). By default, only
critical errors result in output. The verbosity can be set manually.

```@docs
Crystals.Log.set_log_level
```

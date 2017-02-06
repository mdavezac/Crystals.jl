var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Crystals-1",
    "page": "Home",
    "title": "Crystals",
    "category": "section",
    "text": "This package provides a basic framework for dealing with crystalline structures in Julia, from an atomistic point of view. The main object type, Crystal, defines a crystalline structure in terms of its periodicity, as given by a cell, and atomic sites with different properties.Other than the straightforward manipulations of crystalline structures via standard linear algebra and container API, this package defines a few more complex functions, such as methods to create supercells and derived primitive cells, or to compute the point-group and space-group of a crystal.One of the package's main feature/bug/inconvenience is that it mandates users fully specify the physical units (say, nanometers, or even kilowatts per seconds if such is your fancy) of the crystal, thanks to Unitful. However, as a small give-away to practical considerations, it does allow crystal sites to be defined and manipulated in terms of cartesian or fractional coordinates quite transparently."
},

{
    "location": "basic.html#",
    "page": "Basic Usage",
    "title": "Basic Usage",
    "category": "page",
    "text": ""
},

{
    "location": "basic.html#Simple-Construction-1",
    "page": "Basic Usage",
    "title": "Simple Construction",
    "category": "section",
    "text": "A Crystal declares an atomic crystalline structure, e.g. an infinite periodic arrangement of atoms. The constructor takes at the very least an n by n array defining the periodicity of the crystal, i.e. the crystal cell. The cell must have physical units attached to it.using Crystals\ncrystal = Crystal(eye(3)u\"nm\")\n\n# output\ncell(nm):\n  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0However, it can also accept atomic positions and any other array of atomic properties:DocTestSetup = quote\n  using Crystals\nendcrystal = Crystal(eye(2)u\"km\",\n                  position=transpose([1 1; 2 3; 4 5])u\"m\",\n                  species=[\"Al\", \"O\", \"O\"],\n                  label=[:+, :-, :-])\n\n# output\ncell(m):\n  1000.0 0.0\n  0.0 1000.0\n│ Atom │ Cartesian │ species │ label │\n├──────┼───────────┼─────────┼───────┤\n│ 1    │ (1.0,1.0) │ \"Al\"    │ +     │\n│ 2    │ (2.0,3.0) │ \"O\"     │ -     │\n│ 3    │ (4.0,5.0) │ \"O\"     │ -     │Note that the positions form an n by N array where N is the number of atoms. This is the logical mathematical format if when performing matrix operations on the positions. However, from an input point of view, it is easier to think one atom at at time, rather than one coordinate (and all atoms) at a time. Hence the transpose.using Base.Markdown\nusing Crystals\nresult = \"\"\"\n!!! note\n\n    The following are reserved column names that have special meaning. They can be used in\n    most, but not all circumstances. For instance, `:x`, `:y`, and `:x` cannot be used in\n    the constructor.\n\"\"\"\nresult *= \"    - :\" * join(map(string, Crystals.Structures.RESERVED_COLUMNS), \"\\n    - :\")\nBase.Markdown.parse(result)"
},

{
    "location": "basic.html#Accessing-the-cell-and-atomic-sites-1",
    "page": "Basic Usage",
    "title": "Accessing the cell and atomic sites",
    "category": "section",
    "text": "Access to the crystal cell happens via . call, crystal.cell. Atomic properties on the other hand are accessed and modified through the square bracket operator. There are several ways of doing this, more or less reflecting what can be done with a DataFrame:DocTestSetup = quote\n  using Crystals\n  crystal = Crystal(eye(2)u\"km\",\n                  position=transpose([1 1; 2 3; 4 5])u\"m\",\n                  species=[\"Al\", \"O\", \"O\"],\n                  label=[:+, :-, :-])\nendprintln(crystal[:label])\nprintln(crystal[[1, 3], [:species, :label]])\n\ncrystal[1, :position] = [0, 4]u\"m\"\n\n# output\nSymbol[:+,:-,:-]\n2×2 DataFrames.DataFrame\n│ Row │ species │ label │\n├─────┼─────────┼───────┤\n│ 1   │ \"Al\"    │ +     │\n│ 2   │ \"O\"     │ -     │\n2-element Array{Quantity{Int64, Dimensions:{𝐋}, Units:{m}},1}:\n 0 m\n 4 mHowever, the main access pattern is atom/row-based. Using an integer or a sequence of integers will create a new Crystal structure with only the selected atoms.crystal[1]\n\n# output\ncell(m):\n  1000.0 0.0\n  0.0 1000.0\n│ Atom │ Cartesian │ species │ label │\n├──────┼───────────┼─────────┼───────┤\n│ 1    │ (1.0,1.0) │ \"Al\"    │ +     │Note that the return is still a crystalline structure. In a way, we are selecting an atom and it's periodic image, rather than single atom."
},

{
    "location": "basic.html#Base.push!",
    "page": "Basic Usage",
    "title": "Base.push!",
    "category": "Function",
    "text": "push!(crystal, position; kwargs...)\n\n\nAppends an atom to a crystal structure. The position of the atom is a necessary argument, whether in Cartesian or in fractional coordinates. If keyword arguments are present, then they represent atomic properties for the atom being added. Properties that are not explicitly given are set to NA. Similarly, new properties that were not present in the crystal structure previously are NA except for the newly added atom.\n\nExamples\n\nusing Crystals\ncrystal = Crystal(eye(2)u\"km\",\n                  tpositions=[1 1; 2 3; 4 5]u\"m\",\n                  species=[\"Al\", \"O\", \"O\"],\n                  label=[:+, :-, :-])\npush!(crystal, [10, 20]u\"nm\", species=\"B\", μ=1)\ncrystal\n\n# output\n\ncell(m):\n  1000.0 0.0\n  0.0 1000.0\n│ Atom │ Cartesian       │ species │ label │ μ  │\n├──────┼─────────────────┼─────────┼───────┼────┤\n│ 1    │ (1.0,1.0)       │ \"Al\"    │ +     │ NA │\n│ 2    │ (2.0,3.0)       │ \"O\"     │ -     │ NA │\n│ 3    │ (4.0,5.0)       │ \"O\"     │ -     │ NA │\n│ 4    │ (1.0e-8,2.0e-8) │ \"B\"     │ NA    │ 1  │\n\n\n\n"
},

{
    "location": "basic.html#Base.append!",
    "page": "Basic Usage",
    "title": "Base.append!",
    "category": "Function",
    "text": "append!(crystal, other; check_periodicity)\n\n\nAppends one or more crystal structure to the first structure. Unless check_periodicity is false, the structures must have the exact same periodicity. An error will result otherwise.\n\n\n\n"
},

{
    "location": "basic.html#Base.vcat",
    "page": "Basic Usage",
    "title": "Base.vcat",
    "category": "Function",
    "text": "vcat(crystal, other; check_periodicity)\n\n\nConcatenates crystals together. The lattices must be compatible.\n\n\n\n"
},

{
    "location": "basic.html#Base.delete!",
    "page": "Basic Usage",
    "title": "Base.delete!",
    "category": "Function",
    "text": "Base.delete!(crystal::Crystal, col::Symbol)\nBase.delete!(crystal::Crystal, col::AbstractVector{Symbol})\n\nDeletes one or more atomic property. Positions cannot be deleted.\n\n\n\nBase.delete!(crystal::Crystal, ::Colon)\n\nDeletes all atomic properties except for positions.\n\n\n\ndelete!(crystal::Crystal, rows::Integer)\ndelete!(crystal::Crystal, rows::AbstractVector{Integer})\ndelete!{T <: Integer}(crystal::Crystal, rows::Range{T})\ndelete!(crystal::Crystal, rows::Colon)\n\nAlias for deleterows[@ref].\n\n\n\n"
},

{
    "location": "basic.html#DataFrames.deleterows!",
    "page": "Basic Usage",
    "title": "DataFrames.deleterows!",
    "category": "Function",
    "text": "deleterows!(crystal::Crystal, rows::Integer)\ndeleterows!(crystal::Crystal, rows::AbstractVector{Integer})\ndeleterows!{T <: Integer}(crystal::Crystal, rows::Range{T})\ndeleterows!(crystal::Crystal, rows::Colon)\n\nDeletes one (single integer), a few (sequence or range), or all (colon) atoms in the structure.\n\n\n\n"
},

{
    "location": "basic.html#Base.empty!",
    "page": "Basic Usage",
    "title": "Base.empty!",
    "category": "Function",
    "text": "empty!(crystal)\n\n\nDeletes all atomic sites, both properties and positions.\n\n\n\n"
},

{
    "location": "basic.html#Appending-and-removing-atoms-1",
    "page": "Basic Usage",
    "title": "Appending and removing atoms",
    "category": "section",
    "text": "To add atoms to a Crystal instance, it is recommended to use the push!, append!, and vcat functions. Removing atoms and atomic properties from a structure can be done quite easily with delete!. Or alternatively, one can select a few atoms and properties using the bracket notation, e.g. crystal[[1, 2], [:species, :label]].push!\nappend!\nvcat\ndelete!\nDataFrames.deleterows!\nempty!warning: Warning\nManipulating directly the fields positions and properties of a Crystal instance is discouraged. In practice, the only requirement is that the two represent the same number of atoms (with the exception of crystals with no atomic properties outside of the positions, for which the properties field is empty). However, no provision is made anywhere in the code for abnormal Crystal structure. This means the code may crash in creative ways."
},

{
    "location": "basic.html#Iterating-through-a-structure-1",
    "page": "Basic Usage",
    "title": "Iterating through a structure",
    "category": "section",
    "text": "It is a possible to iterate through all the atoms of a structure. The object yielded references directly the parent structure. Atomic properties and positions can be accessed and modified one site at a time.for atom in eachatom(crystal)\n    println(\"Atom: \", atom[:species], \" and \", atom[:label])\n    atom[:species] *= string(atom[:label])\nend\ncrystal\n\n# output\nAtom: Al and +\nAtom: O and -\nAtom: O and -\ncell(m):\n  1000.0 0.0\n  0.0 1000.0\n│ Atom │ Cartesian │ species │ label │\n├──────┼───────────┼─────────┼───────┤\n│ 1    │ (1.0,1.0) │ \"Al+\"   │ +     │\n│ 2    │ (2.0,3.0) │ \"O-\"    │ -     │\n│ 3    │ (4.0,5.0) │ \"O-\"    │ -     │warning: Warning\nModifying the number of sites in the crystal structure during iteration will result in undefined behaviour."
},

{
    "location": "cartesian.html#",
    "page": "Cartesian and Fractional Coordinates",
    "title": "Cartesian and Fractional Coordinates",
    "category": "page",
    "text": "Atomic positions can be initialized, modified, and accessed in Cartesian or fractional coordinates. Cartesian coordinates refer to real world positions in the same physical units as the crystal cell. Fractional coordinates however are in units of the crystal cell."
},

{
    "location": "cartesian.html#Creating-and-accessing-fractional-and-Cartesian-coordinates-1",
    "page": "Cartesian and Fractional Coordinates",
    "title": "Creating and accessing fractional and Cartesian coordinates",
    "category": "section",
    "text": "In the following, we create a crystal structure using fractional coordinates, through the simple expedience of not specifying actual units.CurrentModule = Crystals\nDocTestSetup = quote\n    using Crystals\n    using Unitful\nendfrac_crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u\"nm\",\n                       position=[0, 0, 0],\n                       position=[0.25, 0.25, 0.25])\n@assert frac_crystal[:position] === frac_crystal[:fractional]\n@assert frac_crystal[:cartesian] ≈ frac_crystal.cell * frac_crystal[:fractional]\nunits = unit(eltype(frac_crystal[:cartesian]))\nprintln(ustrip(frac_crystal[:cartesian]), \" (\", units, \")\")\n\n# output\n[0.0 1.05; 0.0 1.05; 0.0 1.05] (nm)Note that querying :position returns fractional coordinates. If we create a structure with Cartesian coordinates instead – by calling the constructor with positions that have units – then querying :position would return the Cartesian coordinates.cart_crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u\"nm\",\n                       position=[0, 0, 0]u\"nm\",\n                       position=[1.05, 1.05, 1.05]u\"nm\")\n@assert cart_crystal[:position] === cart_crystal[:cartesian]\n@assert cart_crystal[:fractional] ≈ inv(cart_crystal.cell) * cart_crystal[:cartesian]\nprintln(cart_crystal[:fractional])\n\n# output\n[0.0 0.25; 0.0 0.25; 0.0 0.25]Of course, in either case, we can access either :cartesian or :fractional coordinates through the relevant column name. However, depending on how the crystal was created, one of these calls will be essentially a no-op, and the other will involve a matrix-matrix multiplication (and possibly computing the inverse of a matrix). To obtain a specific kind of crystal from any other crystal, one can simply use the bracket operator:crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u\"nm\",\n                  position=[0, 0, 0]u\"nm\",\n                  position=[1.05, 1.05, 1.05]u\"nm\",\n                  species=[\"Si\", \"Si\"])\ncrystal[[:fractional, :species]]\n\n# output\ncell(nm):\n  0.0 2.1 2.1\n  2.1 0.0 2.1\n  2.1 2.1 0.0\n│ Atom │ fractional       │ species │\n├──────┼──────────────────┼─────────┤\n│ 1    │ (0.0,0.0,0.0)    │ \"Si\"    │\n│ 2    │ (0.25,0.25,0.25) │ \"Si\"    │Note that the column name explicitly specifies fractional, as opposed to cartesian. This call (indeed, all bracket operator call) will always create a new instance of a Crystal, whether it is strictly needed or not."
},

{
    "location": "cartesian.html#Setting-coordinates-from-fractional-or-Cartesian-inputs-1",
    "page": "Cartesian and Fractional Coordinates",
    "title": "Setting coordinates from fractional or Cartesian inputs",
    "category": "section",
    "text": "When setting a position, it is not necessary to specify whether the input is fractional or Cartesian. If the input has physical units, then it is Cartesian. If it doesn't, then it is fractional. The appropriate transformation is applied before setting the position in the crystal.frac_crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u\"nm\",\n                       position=[0, 0, 0],\n                       position=[0.25, 0.25, 0.25])\nfrac_crystal[2, :position] = [1, 1, 1]u\"nm\"\nfrac_crystal\n\n# output\ncell(nm):\n  0.0 2.1 2.1\n  2.1 0.0 2.1\n  2.1 2.1 0.0\n│ Atom │ fractional                   │\n├──────┼──────────────────────────────┤\n│ 1    │ (0.0,0.0,0.0)                │\n│ 2    │ (0.238095,0.238095,0.238095) │"
},

{
    "location": "cartesian.html#Accessing-and-Modifying-specific-components-1",
    "page": "Cartesian and Fractional Coordinates",
    "title": "Accessing and Modifying specific components",
    "category": "section",
    "text": "Apart from :cartesian and :fractional, there are three other special column names that ease access to specific atomic properties. :x, :y, :z will return an array representing the corresponding coordinate, in the same system – Cartesian or fractional – as the crystal (and as :position). :x, :y, :z  can be used to set a specific coordinate as well. Only three such special names are provided. crystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u\"nm\",\n                  position=[0, 0, 0],\n                  position=[0.25, 0.27, 0.25])\n\nprintln(\"Y coordinate: \", crystal[2, :y])\ncrystal[2, :y] = 0.25\ncrystal\n\n# output\nY coordinate: 0.27\ncell(nm):\n  0.0 2.1 2.1\n  2.1 0.0 2.1\n  2.1 2.1 0.0\n│ Atom │ fractional       │\n├──────┼──────────────────┤\n│ 1    │ (0.0,0.0,0.0)    │\n│ 2    │ (0.25,0.25,0.25) │"
},

{
    "location": "methods.html#",
    "page": "API Catalogue",
    "title": "API Catalogue",
    "category": "page",
    "text": ""
},

{
    "location": "methods.html#Crystals.Structures.is_fractional",
    "page": "API Catalogue",
    "title": "Crystals.Structures.is_fractional",
    "category": "Function",
    "text": "is_fractional(crystal)\n\n\nTrue if the crystal structure is fractional.\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Structures.are_compatible_lattices",
    "page": "API Catalogue",
    "title": "Crystals.Structures.are_compatible_lattices",
    "category": "Function",
    "text": "are_compatible_lattices(lattices...)\n\nTrue if the lattices are mathematically equivalent. Two lattices are equivalent if they represent the same periodicity. In practice, this means the two lattices have the same volume, and their cell vectors are integer linear combinations of one another.\n\nParameters\n\nlattices::Vararg{Union{Matrix, Crystal}}: Any number of lattices\ndigits::Integer: when checking the cells are integer co-linear, the product A^-1B is first rounded to this number of digits\nrtol::Real: relative tolerance when checking the volumes correspond (no units). Default to 1.0e-8.\natol::Real: absolute tolerance when checking the volumes correspond (no units). Default to 1.0e-8.\n\nExamples\n\nusing Crystals, Unitful\ncrystal = Crystal([0 2.1 2.1; 2.1 0 2.1; 2.1 2.1 0]u\"nm\")\ncells = Matrix{Int64}[eye(3)]\nwhile length(cells) < 5\n    cell = rand(-5:5, (3, 3))\n    volume(cell) == 1 && push!(cells, cell)\nend\nlattices = [crystal.cell * c for c in cells]\nare_compatible_lattices(crystal, lattices...)\n\n# output\ntrue\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Structures.volume",
    "page": "API Catalogue",
    "title": "Crystals.Structures.volume",
    "category": "Function",
    "text": "volume(crystal)\n\n\nReturns the volume of a Crystal instance or of a cell. It comes down to computing det(A) where A is the crystal cell.\n\nExamples\n\njulia> using Crystals\n\njulia> volume(Crystal([0 2 2; 2 0 2; 2 2 0]u\"nm\"))\n16.0 nm^3\n\njulia> volume([0 2 2; 2 0 2; 2 2 0]u\"nm\")\n16.0 nm^3\n\njulia> volume([0 2 2; 2 0 2; 2 2 0])\n16.0\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Utilities.cell_parameters",
    "page": "API Catalogue",
    "title": "Crystals.Utilities.cell_parameters",
    "category": "Function",
    "text": "cell_parameters(a::Quantity, b::Quantity, c::Quantity,\n                α::Quantity=(π/2)u\"rad\", β::Quantity=(π/2)u\"rad\",\n                γ::Quantity=(π/2)u\"rad\")\n\nComputes the cell matrix from the cell parameters [a, b, c, α, β, γ].\n\n\n\ncell_parameters(cell::AbstractMatrix)\ncell_parameters(cell::Crystal)\n\nParameters (a, b, c, α, β, γ) of the input cell returned in a named tuple.\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Gruber.gruber",
    "page": "API Catalogue",
    "title": "Crystals.Gruber.gruber",
    "category": "Function",
    "text": "gruber(cell::Matrix;\n       tolerance::Real=default_tolerance, itermax::Unsigned=50,\n       max_no_change::Unsigned=10)\n\nDetermines Gruber cell of an input cell.\n\nThe Gruber cell is an optimal parameterization of a lattice, e.g. shortest cell-vectors and angles closest to 90 degrees. The new cell is in the same basis as the origin cell: no rotation has been incurred. The cell parameters are uniquely determined, even though the cell itself is not (certain symmetry operations leaving the structure unchanged may yield a more recognizable cell). If you want a unique Cartesian cell (in a different Cartesian basis), use the niggly algorithm.\n\nArguments\n\ncell::Matrix: the input lattice cell-vectors. Cannot be singular.\nitermax::Integer: maximum number of iterations before bailing out\ntolerance::Number: tolerance parameter when comparing real numbers\nmax_no_change::Integer: Maximum number of times to go through algorithm without changes before bailing out\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Gruber.niggly",
    "page": "API Catalogue",
    "title": "Crystals.Gruber.niggly",
    "category": "Function",
    "text": "niggly(cell::Matrix; kwargs...)\n\nDetermines a unique Cartesian cell equivalent to the input, with the shortest possible vectors and squarest angles. For an explanation of the parameters, see gruber. In practice, this function computes the cell-parameters of a gruber cell and then reconstructs the cell matrix. Hence, the result should be quite unique for any lattice representation, including any rotation of the underlying Cartesian coordinates.\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Utilities.is_periodic",
    "page": "API Catalogue",
    "title": "Crystals.Utilities.is_periodic",
    "category": "Function",
    "text": "is_periodic(a::AbstractVector, b::AbstractVector, cell::AbstractMatrix;\n            tolerance::Real=1.0e-8)\n\nTrue if the positions are one-to-one periodic with respect to the input cell.\n\n\n\nis_periodic(a::AbstractMatrix, b::AbstractVector, cell::AbstractMatrix;\n            tolerance::Real=1.0e-8)\n\nArray of boolean describing whether positions in a are periodic with positions in b.\n\n\n\n"
},

{
    "location": "methods.html#Crystals.SpaceGroup.is_primitive",
    "page": "API Catalogue",
    "title": "Crystals.SpaceGroup.is_primitive",
    "category": "Function",
    "text": "is_primitive(cartesian::AbstractMatrix, cell::AbstractMatrix, species::AbstractVector;\n             tolerance::Real=1.0e-8)\nis_primitive(crystal::Crystal, col::Union{Symbol, AbstractVector{Symbol}}; kwargs...)\nis_primitive(crystal::Crystal; kwargs...)\n\nTrue if the crystal structure is primitive, e.g. not a supercell, e.g. not reducible to an equivalent lattice with fewer atoms.\n\n\n\n"
},

{
    "location": "methods.html#Querying-and-property-methods-1",
    "page": "API Catalogue",
    "title": "Querying and property methods",
    "category": "section",
    "text": "CurrentModule = Crystalsis_fractional\nare_compatible_lattices\nvolume\ncell_parameters\ngruber\nniggly\nis_periodic\nis_primitiveTypical container and DataFrame methods are also available, such as names, size, ndims, nrow, ncol, endof."
},

{
    "location": "methods.html#Crystals.CrystalAtoms.eachatom",
    "page": "API Catalogue",
    "title": "Crystals.CrystalAtoms.eachatom",
    "category": "Function",
    "text": "Iterator over each atom in the crystal \n\n\n\n"
},

{
    "location": "methods.html#Looping-1",
    "page": "API Catalogue",
    "title": "Looping",
    "category": "section",
    "text": "eachatomeachindex is also available. It returns the range over atoms indices 1:size(crystal, 1)."
},

{
    "location": "methods.html#Base.round",
    "page": "API Catalogue",
    "title": "Base.round",
    "category": "Function",
    "text": "round(crystal, args)\n\n\nRounds the cell and positions of a crystal. See Base.round for possible parameters.\n\nExamples\n\nusing Crystals\ncrystal = Crystal([0 0.501 0.501; 0.496 0.001 0.497; 0.497 0.497 0]u\"nm\",\n                  position=[0.001, -0.001, -0.001]u\"nm\",\n                  position=[0.25, 0.251, -0.247]u\"nm\")\nround(crystal, 2)\n\n# output\n\ncell(nm):\n  0.0 0.5 0.5\n  0.5 0.0 0.5\n  0.5 0.5 0.0\n│ Atom │ Cartesian         │\n├──────┼───────────────────┤\n│ 1    │ (0.0,-0.0,-0.0)   │\n│ 2    │ (0.25,0.25,-0.25) │\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Structures.round!",
    "page": "API Catalogue",
    "title": "Crystals.Structures.round!",
    "category": "Function",
    "text": "round!(crystal, args)\n\n\nRounds the cell and positions of a crystal. See round for possible parameters.\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Utilities.into_cell",
    "page": "API Catalogue",
    "title": "Crystals.Utilities.into_cell",
    "category": "Function",
    "text": "into_cell(pos, cell; tolerance)\n\n\nFolds periodic positions into cell\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Utilities.into_voronoi",
    "page": "API Catalogue",
    "title": "Crystals.Utilities.into_voronoi",
    "category": "Function",
    "text": "into_voronoi(positions::AbstractArray, cell::AbstractMatrix; extent::Integer=1)\n\nFolds positions into first Brillouin zone of the input cell. Makes a well-meaning effort at returning the periodic image with the smallest possible norm. It recenter the atoms around the origin and then looks for the smallest periodic images within -extent:extent cells. If the cell is quite pathological, then the result will not be within the Voronoi cell.\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Utilities.origin_centered",
    "page": "API Catalogue",
    "title": "Crystals.Utilities.origin_centered",
    "category": "Function",
    "text": "origin_centered(positions::AbstractArrays, cell::AbstractMatrix)\n\nFolds positions back to origin, such that each fractional component x_f is between -05leq x_f  05. If the input is in Cartesian (fractional) coordinates, then the result is also in Cartesian (fractional) coordinates.\n\n\n\n"
},

{
    "location": "methods.html#Modifying,-building-up-and-building-down-1",
    "page": "API Catalogue",
    "title": "Modifying, building up and building down",
    "category": "section",
    "text": "round\nround!\ninto_cell\ninto_voronoi\norigin_centeredis_primitive can be queried to figure out whether a crystal is a primitive cell or a supercell."
},

{
    "location": "methods.html#Crystals.SpaceGroup.point_group",
    "page": "API Catalogue",
    "title": "Crystals.SpaceGroup.point_group",
    "category": "Function",
    "text": "Finds and stores point group operations for a given lattice\n\nA lattice is defined by a 3x3 matrix or cell.  Rotations are determined from G-vector triplets with the same norm as the unit-cell vectors.\n\nImplementation taken from ENUM.\n\n\n\n"
},

{
    "location": "methods.html#Crystals.SpaceGroup.space_group",
    "page": "API Catalogue",
    "title": "Crystals.SpaceGroup.space_group",
    "category": "Function",
    "text": "space_group(crystal; kwargs...)\n\n\nComputes the space-group operations of a crystal. By default, all atomic properties are considered when determining whether atomic sites are related by symmetry. However, it is possible to specify a subset of atomic properties. An empty subset of atomic properties implies that all atoms sites are equivalent. It is also possible to specify an array of integers, serving as labels for each atomic site.\n\n\n\nspace_group(cell, positions, species; kwargs...)\n\n\nComputes the space-group of a crystal specified using standard Julia types.\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Utilities.hart_forcade",
    "page": "API Catalogue",
    "title": "Crystals.Utilities.hart_forcade",
    "category": "Function",
    "text": "hart_forcade(lattice, supercell; digits)\n\n\nComputes the cyclic group of a supercell with respect to a lattice. It makes it possible to identify the class of periodically equivalent cell that a given position within the supercell belongs to. The function returns a named tuple with the transform and the quotient.\n\nExamples\n\nusing Crystals, Unitful\nfcc = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u\"nm\"\nsupercell = [0 2 2; 0 -4 2; -1 0 -2]\nht = hart_forcade(fcc, fcc * supercell)\n\nprintln(ht)\nprintln(\"Positions in supercell:\")\nfor index in CartesianRange((ht.quotient...))\n    position = inv(ht.transform) * [index[u] for u in eachindex(ht.quotient)]\n    println(\"- \", ustrip(position), \" (\", unit(eltype(position)), \")\")\nend\n\n# output\n\nHart-Forcade transform\n- transform (nm^-1): [-1.0 -1.0 1.0; -1.0 1.0 1.0; -1.0 1.0 3.0]\n- quotient: [1,2,6]\n\nPositions in supercell:\n- [-1.0,0.0,0.0] (nm)\n- [-2.0,0.5,-0.5] (nm)\n- [-0.5,0.0,0.5] (nm)\n- [-1.5,0.5,0.0] (nm)\n- [0.0,0.0,1.0] (nm)\n- [-1.0,0.5,0.5] (nm)\n- [0.5,0.0,1.5] (nm)\n- [-0.5,0.5,1.0] (nm)\n- [1.0,0.0,2.0] (nm)\n- [0.0,0.5,1.5] (nm)\n- [1.5,0.0,2.5] (nm)\n- [0.5,0.5,2.0] (nm)\n\n\n\n"
},

{
    "location": "methods.html#Crystals.SpaceGroup.primitive",
    "page": "API Catalogue",
    "title": "Crystals.SpaceGroup.primitive",
    "category": "Function",
    "text": "primitive(crystal::Crystal; tolerance::Real=1.0e-8)\n\nComputes the primitive cell of the input crystal. If the crystal is primitive, it is returned as is, otherwise a new crystal is returned.\n\n\n\n"
},

{
    "location": "methods.html#Crystals.Utilities.supercell",
    "page": "API Catalogue",
    "title": "Crystals.Utilities.supercell",
    "category": "Function",
    "text": "supercell(lattice, supercell; site_id, tolerance)\n\n\nCreates a supercell from an input lattice.\n\n# Parameters\n\nlattice::Crystal: the original lattice\nsupercell::AbstractMatrix: the cell of the supercell in Cartesian coordinates\nsite_id::Bool: Whether to add/modify an atomic property indicating the index of the site in the original lattice\ncell_id::Bool: Whether to add/modify an atomic property indicating the index of the cell the site belongs to\n\n\n\n"
},

{
    "location": "methods.html#Point-group,-space-group,-and-geometry-1",
    "page": "API Catalogue",
    "title": "Point-group, space-group, and geometry",
    "category": "section",
    "text": "point_group\nspace_group\nhart_forcade\nprimitive\nsupercell"
},

{
    "location": "methods.html#Crystals.Lattices.b5",
    "page": "API Catalogue",
    "title": "Crystals.Lattices.b5",
    "category": "Function",
    "text": "b5(T; unit)\nb5()\n\n\nb5 (spinel) lattice \n\n\n\n"
},

{
    "location": "methods.html#Crystals.Lattices.bcc",
    "page": "API Catalogue",
    "title": "Crystals.Lattices.bcc",
    "category": "Function",
    "text": "bcc()\nbcc(T; unit)\n\n\nCreates a body-centered lattice with a single site \n\n\n\n"
},

{
    "location": "methods.html#Crystals.Lattices.diamond",
    "page": "API Catalogue",
    "title": "Crystals.Lattices.diamond",
    "category": "Function",
    "text": "diamond()\ndiamond(T; unit)\n\n\nDiamond lattice \n\n\n\n"
},

{
    "location": "methods.html#Crystals.Lattices.fcc",
    "page": "API Catalogue",
    "title": "Crystals.Lattices.fcc",
    "category": "Function",
    "text": "fcc()\nfcc(T; unit)\n\n\nCreates an face centered lattice with a single site \n\n\n\n"
},

{
    "location": "methods.html#Crystals.Lattices.rock_salt",
    "page": "API Catalogue",
    "title": "Crystals.Lattices.rock_salt",
    "category": "Function",
    "text": "rock_salt(T; unit)\nrock_salt()\n\n\nRock-salt lattice \n\n\n\n"
},

{
    "location": "methods.html#Crystals.Lattices.wurtzite",
    "page": "API Catalogue",
    "title": "Crystals.Lattices.wurtzite",
    "category": "Function",
    "text": "wurtzite()\nwurtzite(T; unit)\n\n\nWurtzite lattice \n\n\n\n"
},

{
    "location": "methods.html#Crystals.Lattices.zinc_blende",
    "page": "API Catalogue",
    "title": "Crystals.Lattices.zinc_blende",
    "category": "Function",
    "text": "zinc_blende()\nzinc_blende(T; unit)\n\n\nZinc-blende lattice \n\n\n\n"
},

{
    "location": "methods.html#Predefined-lattices-1",
    "page": "API Catalogue",
    "title": "Predefined lattices",
    "category": "section",
    "text": "A small number of standard lattices are available for construction in the exported Lattices submodule.For instance, the body-centered lattice:using Crystals\nLattices.bcc()Or more complex lattices, such as b5 spinels:Lattices.b5()Modules = [Crystals.Lattices]"
},

{
    "location": "methods.html#Crystals.Log.set_log_level",
    "page": "API Catalogue",
    "title": "Crystals.Log.set_log_level",
    "category": "Function",
    "text": "set_log_level(level)\nset_log_level()\n\n\nModifies log-level of all \"trucks\" in Crystals logs. The input should be one of \"debug\", \"info\", \"warn\", \"error\", from least to most verbose.\n\n\n\n"
},

{
    "location": "methods.html#Output-and-Logging-1",
    "page": "API Catalogue",
    "title": "Output and Logging",
    "category": "section",
    "text": "Information and errors during calculations are displayed using an internal log provided by Lumberjack. By default, only critical errors result in output. The verbosity can be set manually.Crystals.Log.set_log_level"
},

]}

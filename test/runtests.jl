module CrystalTest
using Crystal
using FactCheck: @fact, facts, context, exitstatus, roughly
using DataFrames: nrow, ncol, NA, DataArray

contains(x) = y -> x ∈ y

facts("Construction") do

  for ndim in [2, 3]
      structure = Structure(eye(ndim))
      @fact abs(structure.cell - eye(ndim)) .< 1e-8 --> all
      @fact length(structure.positions) --> 0
      @fact length(structure) --> 2
      @fact ndims(structure) --> ndim
      @fact nrow(structure.properties) --> 0
      @fact ncol(structure.properties) --> 1
      @fact names(structure.properties) --> contains(:specie)
      @fact structure.scale --> roughly(1e0)
  end
end

for ndim in [2, 3]
  # Does not make use of indexing interface, for better test separation
  # This is not a good example on how to use Structure
  facts("Pushing values") do


    context("Push atom without attributes") do
      structure = Structure(eye(ndim))
      push!(structure, ones(ndim), "Al")
      @fact size(structure.positions, 2) --> 1
      @fact size(structure.positions, 1) --> ndim
      @fact nrow(structure.properties) --> 1
      @fact ncol(structure.properties) --> 1
      @fact ndims(structure) --> ndim
      @fact length(structure) --> 2
      @fact structure.positions[:, 1] --> roughly(ones(ndim))
      @fact structure.properties[1, :specie] --> "Al"
    end


    context("Push atom with attributes") do
      structure = Structure(eye(ndim))
      push!(structure, ones(ndim), "Al"; fruit="abricot")
      @fact size(structure.positions, 2) --> 1
      @fact size(structure.positions, 1) --> ndim
      @fact structure.positions[:, 1] --> roughly(ones(ndim))
      @fact ndims(structure) --> ndim
      @fact length(structure) --> 3

      @fact nrow(structure.properties) --> 1
      @fact ncol(structure.properties) --> 2
      @fact structure.positions[:, 1] --> roughly(ones(ndim))
      @fact structure.properties[1, :specie] --> "Al"
      @fact structure.properties[1, :fruit] --> "abricot"
    end
  end

  facts("Indexing") do
    structure = Structure(5eye(ndim), 10)
    push!(structure, [0, 0, 0][1:ndim], "C"; charge=4)
    push!(structure, [0.25, 0.25, 0.25][1:ndim], "O"; charge=-2, label=1)
    push!(structure, -[0.25, 0.25, 0.25][1:ndim], "O"; charge=-2, label=2)

    context("Column") do
      @fact structure[:specie] --> ["C", "O", "O"]
      @fact structure[1] --> structure[:specie]
      @fact structure[:charge] --> [4, -2, -2]
      @fact structure[2] --> structure[:charge]
      atoms = transpose([0 0 0; 0.25 0.25 0.25; -0.25 -0.25 -0.25])
      @fact structure[:position] --> roughly(atoms[1:ndim, :])
      @fact structure[end] --> roughly(atoms[1:ndim, :])
      @fact structure[4] --> roughly(atoms[1:ndim, :])
    end

    context("Columns") do
      partial = structure[[:charge, :label]]
      @fact names(partial) --> [:charge, :label, :position]
      @fact partial.cell --> roughly(structure.cell)
      @fact partial.scale --> roughly(structure.scale)
      @fact partial[:charge] --> structure[:charge]
      @fact partial[:position] --> structure[:position]

      full = structure[:]
      @fact names(full) --> names(structure)
      @fact full.cell --> roughly(structure.cell)
      @fact full.scale --> roughly(structure.scale)
      @fact full[:specie] --> structure[:specie]
      @fact full[:charge] --> structure[:charge]
      @fact full[:position] --> structure[:position]
      @fact full --> x -> x ≢ structure
    end
  end
end

exitstatus()
end

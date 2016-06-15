module CrystalTest
using Crystal
using FactCheck: @fact, facts, context, exitstatus, roughly
using DataFrames: nrow, ncol

contains(x) = y -> x ∈ y

facts("Construction") do

  for ndim in [2, 3]
      structure = Structure(eye(ndim))
      @fact abs(structure.cell - eye(ndim)) .< 1e-8 --> all
      @fact length(structure.positions) --> 0
      @fact nrow(structure.properties) --> 0
      @fact ncol(structure.properties) --> 1
      @fact names(structure.properties) --> contains(:specie)
      @fact length(structure.scale) --> 1e0
      @fact length(structure) --> 0
      @fact ndims(structure) --> ndim
  end
end

# Does not make use of indexing interface, for better test separation
# This is not a good example on how to use Structure
facts("Pushing values") do

  for ndim in [2, 3]

    context("Push atom without attributes") do
      structure = Structure(eye(ndim))
      push!(structure, ones(ndim), "Al")
      @fact size(structure.positions, 2) --> 1
      @fact size(structure.positions, 1) --> ndim
      @fact nrow(structure.properties) --> 1
      @fact ncol(structure.properties) --> 1
      @fact ndims(structure) --> ndim
      @fact length(structure) --> 1
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
      @fact length(structure) --> 1

      @fact nrow(structure.properties) --> 1
      @fact ncol(structure.properties) --> 2
      @fact structure.positions[:, 1] --> roughly(ones(ndim))
      @fact structure.properties[1, :specie] --> "Al"
      @fact structure.properties[1, :fruit] --> "abricot"
    end
  end
end

exitstatus()
end

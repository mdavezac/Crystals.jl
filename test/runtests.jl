module CrystalTest
using Crystals
using FactCheck: @fact, facts, context, exitstatus, roughly, exactly
using DataFrames: nrow, ncol, NA, DataArray

contains(x) = y -> x ∈ y

facts("Positions") do
  positions = Positions([1 2; 3 4])
  @fact length(positions) --> 2
  @fact positions[1] --> [1, 3]
  @fact positions[2] --> [2, 4]
end

facts("Construction") do
  context("Empty 2d") do
    crystal = Crystal(eye(2), 2e0)
    @fact nrow(crystal.atoms) --> 0
    @fact names(crystal.atoms) --> [:position]
  end

  context("3d with only positions") do
    crystal = Crystal(eye(3), position=transpose([1 1 1; 2 3 4]))
    @fact nrow(crystal.atoms) --> 2
    @fact names(crystal.atoms) --> [:position]
    @fact crystal.atoms[1, :position] --> [1, 1, 1]
    @fact crystal.atoms[2, :position].x --> 2
    @fact crystal.atoms[2, :position].y --> 3
    @fact crystal.atoms[2, :position].z --> 4
  end

  context("2d with positions and species") do
    crystal = Crystal(eye(2), position=[1 1 1; 2 3 4], specie=["Al", "O", "O"])
    @fact nrow(crystal.atoms) --> 3
    @fact names(crystal.atoms) --> [:position, :specie]

    # Index of position column depends on input
    crystal = Crystal(eye(2), specie=["Al", "O", "O"], position=[1 1 1; 2 3 4])
    @fact names(crystal.atoms) --> [:specie, :position]
  end

  context("4d, constructing via arguments") do
    crystal = Crystal(eye(4), Any[["Al", "O"], [1 1; 2 2; 3 3; 4 4]],
                      [:specie, :position])
    @fact names(crystal.atoms) --> [:specie, :position]
    @fact nrow(crystal.atoms) --> 2
    @fact eltype(crystal.atoms[:position]) --> exactly(Crystals.Position4D{Int64})
  end
end

facts("Check direct indexing") do
  context("getindex") do
    crystal = Crystal(eye(2), specie=["Al", "O", "O"],
                      position=[1 1 1; 2 3 4], label=[:+, :-, :-])
    @fact crystal[:label] --> [:+, :-, :-]
    @fact crystal[1, :position] --> [1, 2]
    @fact crystal[2, :position] --> [1, 3]
    @fact crystal[end, :position] --> [1, 4]
    @fact crystal[:, end] --> [:+, :-, :-]
    @fact crystal[1, [:position, :label]] -->
        Crystal(eye(2), position=[1, 2], label=[:+]).atoms
    @fact crystal[[1, 3], [:position, :label]] -->
        Crystal(eye(2), position=[1 1; 2 4], label=[:+, :-]).atoms
    @fact crystal[:] --> crystal.atoms
  end

  context("setindex!") do
    crystal = Crystal(eye(3), specie=["Al", "O"],
                      position=transpose([1 1 1; 2 3 4]),
                      label=[:+, :-])
    context("Single column") do
      crystal[:label] = [:z, :a]
      @fact crystal[:label] --> [:z, :a]
      crystal[3] = [:Z, :A]
      @fact crystal[:label] --> [:Z, :A]

      crystal[:position] = transpose([1 3 4; 2 2 6])
      @fact crystal[1, :position] --> [1, 3, 4]
      @fact crystal[2, :position] --> [2, 2, 6]

      crystal[2] = transpose([6 2 8; 4 3 2])
      @fact crystal[1, :position] --> [6, 2, 8]
      @fact crystal[2, :position] --> [4, 3, 2]
    end

    context("Multi-column") do
      original = deepcopy(crystal)
      other = Crystal(eye(3), specie=["Ru", "Ta"],
                        position=transpose([2 4 6; 4 1 2]),
                        label=[:a, :b])
      crystal[[:specie, :label]] = other
      @fact crystal[:specie] --> exactly(other[:specie])
      @fact crystal[:label] --> exactly(other[:label])

      crystal[[:specie, :label]] = original.atoms
      @fact crystal[:specie] --> exactly(original[:specie])
      @fact crystal[:label] --> exactly(original[:label])

      crystal[[false, true]] = other
      @fact crystal[:specie] --> exactly(original[:specie])
      @fact crystal[:label] --> exactly(original[:label])
      @fact crystal[:position] --> exactly(other[:position])

      crystal[[false, true]] = transpose([1 1 2; 3 3 2])
      @fact crystal[1, :position] --> Crystals.Position3D(1, 1, 2)
      @fact crystal[2, :position] --> Crystals.Position3D(3, 3, 2)

      crystal[[:this, :specie]] = "Al"
      @fact crystal[:this] --> ["Al", "Al"]
      @fact crystal[:specie] --> ["Al", "Al"]
    end
  end
end

exitstatus()
end

facts("Convertion vector <--> position") do
   array = [1, 2, 3]
   position = convert(Position, array)
   @fact typeof(position) --> exactly(Crystals.Position3D{eltype(array)})
   @fact [position...] --> array
   @fact typeof(convert(Array, position)) --> typeof(array)
   @fact convert(Array, position) --> array
end

facts("Convertions matrix <--> array of positions") do
  context("Automatic conversion") do
    matrix = transpose([1 2 3; 4 5 6])
    positions = convert(PositionArray, matrix)

    @fact eltype(positions) --> exactly(Crystals.Position3D{eltype(matrix)})
    @fact typeof(positions) -->
    exactly(Vector{Crystals.Position3D{eltype(matrix)}})
    @fact length(positions) --> size(matrix, 2)
    @fact length(positions[1]) --> size(matrix, 1)
    for i = 1:length(positions)
      @fact positions[i] --> matrix[:, i]
    end

    back = convert(Array, positions)
    @fact typeof(back) --> typeof(matrix)
    @fact size(back) --> size(matrix)
    @fact back --> roughly(matrix)

    to_integer = convert(PositionArray, [1 2; 2 3])
    @fact eltype(eltype(to_integer)) <: Integer --> true
    to_float = convert(PositionArray, [1 2.; 2 3])
    @fact eltype(eltype(to_float)) <: AbstractFloat  --> true
  end

  context("Explicit element type conversion") do
    positions = convert(Vector{Crystals.Position2D{Int8}}, [1 2; 3 4])
    @fact eltype(eltype(positions)) --> Int8
  end
end

facts("Convertions dataarray <--> array of positions") do
  context("Automatic conversion") do
    matrix = transpose([1 2 3; 4 5 6])
    positions = convert(PositionDataArray, matrix)

    @fact eltype(positions) --> exactly(Crystals.Position3D{eltype(matrix)})
    @fact typeof(positions) -->
        exactly(DataArray{Crystals.Position3D{eltype(matrix)}, 1})
    @fact length(positions) --> size(matrix, 2)
    @fact length(positions[1]) --> size(matrix, 1)
    for i = 1:length(positions)
      @fact positions[i] --> matrix[:, i]
    end

    back = convert(Array, positions)
    @fact typeof(back) --> typeof(matrix)
    @fact size(back) --> size(matrix)
    @fact back --> roughly(matrix)

    to_integer = convert(PositionArray, [1 2; 2 3])
    @fact eltype(eltype(to_integer)) <: Integer --> true
    to_float = convert(PositionArray, [1 2.; 2 3])
    @fact eltype(eltype(to_float)) <: AbstractFloat  --> true
  end
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
    crystal = Crystal(eye(2), position=[1 1 1; 2 3 4], species=["Al", "O", "O"])
    @fact nrow(crystal.atoms) --> 3
    @fact names(crystal.atoms) --> [:position, :species]

    # Index of position column depends on input
    crystal = Crystal(eye(2), species=["Al", "O", "O"], position=[1 1 1; 2 3 4])
    @fact names(crystal.atoms) --> [:species, :position]
  end

  context("4d, constructing via arguments") do
    crystal = Crystal(eye(4), Any[["Al", "O"], [1 1; 2 2; 3 3; 4 4]],
                      [:species, :position])
    @fact names(crystal.atoms) --> [:species, :position]
    @fact nrow(crystal.atoms) --> 2
    @fact eltype(crystal.atoms[:position]) --> exactly(Crystals.Position4D{Int64})
  end
end

facts("Check direct indexing") do
  context("getindex") do
    crystal = Crystal(eye(2), species=["Al", "O", "O"],
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
    crystal = Crystal(eye(3), species=["Al", "O"],
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
      other = Crystal(eye(3), species=["Ru", "Ta"],
                        position=transpose([2 4 6; 4 1 2]),
                        label=[:a, :b])
      crystal[[:species, :label]] = other[[:species, :label]]
      @fact crystal[:species] --> exactly(other[:species])
      @fact crystal[:label] --> exactly(other[:label])

      crystal[[:species, :label]] = original.atoms[[:species, :label]]
      @fact crystal[:species] --> exactly(original[:species])
      @fact crystal[:label] --> exactly(original[:label])

      crystal[[false, true]] = other[[:position]]
      @fact crystal[:species] --> exactly(original[:species])
      @fact crystal[:label] --> exactly(original[:label])
      @fact crystal[:position] --> exactly(other[:position])

      crystal[[false, true]] = transpose([1 1 2; 3 3 2])
      @fact crystal[1, :position] --> Crystals.Position3D(1, 1, 2)
      @fact crystal[2, :position] --> Crystals.Position3D(3, 3, 2)

      crystal[[:extra_column, :species]] = "Al"
      @fact crystal[:extra_column] --> ["Al", "Al"]
      @fact crystal[:species] --> ["Al", "Al"]
    end

    context("Single-row, Single-Column") do
      crystal[1, :label] = :aa
      @fact crystal[1, :label] --> :aa

      crystal[1, :position] = [1, 2, 3]
      @fact crystal[1, :position] --> [1, 2, 3]
      @fact typeof(crystal[1, :position]) --> x -> x <: Crystals.Position

      crystal[1, :position] = 1
      @fact crystal[1, :position] --> [1, 1, 1]

      crystal[1, :position] = :2
      @fact crystal[1, :position] --> [2, 2, 2]
      crystal[1, :position] = :2, :3, :4
      @fact crystal[1, :position] --> [2, 3, 4]

      @fact_throws MethodError crystal[1, :position] = :a
      @fact_throws ErrorException crystal[1, :nonexistent] = "Al"
    end

    context("Single-row, Multi-Column") do
      crystal[1, [:species, :label]] = "aha"
      @fact crystal[1, :species] --> "aha"
      @fact crystal[1, :label] --> :aha

      crystal[2, [true, false]] = "zha"
      @fact crystal[2, :species] --> "zha"
      @fact_throws ErrorException crystal[1, [:species, :nonexistent]] = "Al"
    end

    context("Multi-row, single-Column") do
      crystal = Crystal(eye(3), species=["Al", "O", "O"],
                      position=transpose([1 1 1; 2 3 4; 4 5 2]),
                      label=[:+, :-, :0])
      crystal[[1, 3], :species] = "Ala"
      @fact crystal[:species] --> ["Ala", "O", "Ala"]

      crystal[[2, 3], :species] = ["H", "B"]
      @fact crystal[:species] --> ["Ala", "H", "B"]

      crystal[[true, false, true], 2] = transpose([1 2 3; 4 5 6])
      @fact crystal[1, :position] --> [1, 2, 3]
      @fact crystal[2, :position] --> [2, 3, 4]
      @fact crystal[3, :position] --> [4, 5, 6]
    end

    context("Multi-row, multi-Column") do
      crystal = Crystal(eye(3), species=["Al", "O", "O"],
                      position=transpose([1 1 1; 2 3 4; 4 5 2]),
                      label=[:+, :-, :0])
      other = Crystal(eye(3), species=["H", "B", "C"],
                      position=transpose([2 1 2; 3 4 3; 5 2 5]),
                      label=[:a, :b, :c])
      crystal[[true, false, true], [:label, :position]] =
          other[1:2, [3, 2]]
      @fact crystal[:label] --> [:a, :-, :b]
      @fact crystal[1, :position] --> [2, 1, 2]
      @fact crystal[2, :position] --> [2, 3, 4]
      @fact crystal[3, :position] --> [3, 4, 3]

      crystal[[false, true, true], [:label, :species]] = ["aa", "bb"]
      @fact crystal[:label] --> [:a, "aa", "bb"]
      @fact crystal[:species] --> ["Al", "aa", "bb"]

      crystal[1:2, [:label, :species]] = "A"
      @fact crystal[:label] --> ["A", "A", "bb"]
      @fact crystal[:species] --> ["A", "A", "bb"]
    end

  end
end

facts("Mutating functions") do
  crystal = Crystal(eye(3), species=["Al", "O", "O"],
                  position=transpose([1 1 1; 2 3 4; 4 5 2]),
                  label=[:+, :-, :0])
  other = Crystal(eye(3), other=["H", "B", "C"],
                  position=transpose([2 1 2; 3 4 3; 5 2 5]),
                  label=[:a, :b, :c])
  df = DataFrame(A=1:3, B=6:8)
  original = deepcopy(crystal)

  merge!(crystal, other.atoms, df)
  @fact names(crystal) --> [:species, :position, :label, :other, :A, :B]
  @fact crystal[:species] --> ["Al", "O", "O"]
  @fact crystal[:position] --> exactly(other[:position])
  @fact crystal[:label] --> exactly(other[:label])
  @fact crystal[:A] --> exactly(df[:A])
  @fact crystal[:B] --> exactly(df[:B])

  delete!(crystal, [:A, :B])
  @fact names(crystal) --> [:species, :position, :label, :other]

  deleterows!(crystal, 3)
  @fact crystal[:species] --> ["Al", "O"]
  @fact crystal[:label] --> [:a, :b]
end

facts("Adding atoms") do
  crystal = Crystal(eye(3), position=transpose([1 1 1; 2 3 4; 4 5 2]))
  push!(crystal, ([1, 4, 2], ))
  @fact size(crystal, 1) --> 4
  @fact crystal[4, :position] --> [1, 4, 2]

  # Add via tuple
  crystal = Crystal(eye(3), species=["Al", "O", "O"],
                  position=transpose([1 1 1; 2 3 4; 4 5 2]),
                  label=[:+, :-, :0])
  push!(crystal, ("B", [1, 4, 2], :aa))
  @fact size(crystal, 1) --> 4
  @fact crystal[4, :species] --> "B"
  @fact crystal[4, :position] --> [1, 4, 2]
  @fact crystal[4, :label] --> :aa

  # Add only position, everything else NA
  push!(crystal, [5, 5, 3])
  @fact size(crystal, 1) --> 5
  @fact crystal[5, :species] --> exactly(NA)
  @fact crystal[5, :position] --> [5, 5, 3]
  @fact crystal[5, :label] --> exactly(NA)

  # Add a new column
  push!(crystal, [6, 6, 2], special=:special, species="BaBa")
  @fact size(crystal, 1) --> 6
  @fact names(crystal) --> [:species, :position, :label, :special]
  @fact crystal[6, :species] --> "BaBa"
  @fact crystal[6, :position] --> [6, 6, 2]
  for i in 1:5
    @fact crystal[i, :special] --> exactly(NA)
  end
  @fact crystal[6, :special] --> :special
  @fact crystal[6, :label] --> exactly(NA)
end

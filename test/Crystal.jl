@testset "> Construction" begin
  @testset ">> Empty 2d" begin
    crystal = Crystal(eye(2)u"nm")
    @test nrow(crystal.atoms) == 0
    @test names(crystal.atoms) == [:position]
  end

  @testset ">> 3d with only positions" begin
    crystal = Crystal(eye(3)u"nm", position=transpose([1 1 1; 2 3 4])u"nm")
    @test nrow(crystal.atoms) == 2
    @test names(crystal.atoms) == [:position]
    @test crystal.atoms[1, :position] == [1, 1, 1]u"nm"
    @test crystal.atoms[2, :position][1] == 2u"nm"
    @test crystal.atoms[2, :position][2] == 3u"nm"
    @test crystal.atoms[2, :position][3] == 4u"nm"
  end

  @testset ">> 3d with fractional" begin
      crystal = Crystal(eye(3)u"nm", position=transpose([1 1 1; 2 3 4]))
      @test nrow(crystal.atoms) == 2
      @test names(crystal.atoms) == [:position]
      @test crystal.atoms[1, :position] == [1, 1, 1]
      @test crystal.atoms[2, :position][1] == 2
      @test crystal.atoms[2, :position][2] == 3
      @test crystal.atoms[2, :position][3] == 4
  end

  @testset ">> 2d with positions and species" begin
    crystal = Crystal(
        eye(2)u"nm", position=[1 1 1; 2 3 4]u"nm", species=["Al", "O", "O"])
    @test nrow(crystal.atoms) == 3
    @test names(crystal.atoms) == [:position, :species]

    # Index of position column depends on input
    crystal = Crystal(
        eye(2)u"nm", species=["Al", "O", "O"], position=[1 1 1; 2 3 4]u"nm")
    @test names(crystal.atoms) == [:species, :position]
  end

  @testset ">> 4d, constructing via arguments" begin
    crystal = Crystal(eye(4)u"nm", Any[["Al", "O"], [1 1; 2 2; 3 3; 4 4]u"nm"],
                      [:species, :position])
    @test names(crystal.atoms) == [:species, :position]
    @test nrow(crystal.atoms) == 2
    @test eltype(crystal.atoms[:position]) === Position{eltype([1, 1]u"nm"), 4}
  end
end

@testset "> Check direct indexing" begin
  @testset ">> getindex" begin
    crystal = Crystal(eye(2)u"nm", species=["Al", "O", "O"],
                      position=[1 1 1; 2 3 4], label=[:+, :-, :-])
    @test crystal[:label] == [:+, :-, :-]
    @test crystal[1, :position] == [1, 2]
    @test crystal[2, :position] == [1, 3]
    @test crystal[end, :position] == [1, 4]
    @test crystal[:, end] == [:+, :-, :-]
    @test crystal[1, [:position, :label]] ==
        Crystal(eye(2)u"nm", position=[1, 2], label=[:+]).atoms
    @test crystal[[1, 3], [:position, :label]] ==
        Crystal(eye(2)u"nm", position=[1 1; 2 4], label=[:+, :-]).atoms
    @test crystal[:] == crystal.atoms
  end

  @testset ">> setindex!" begin
    crystal = Crystal(eye(3)u"nm", species=["Al", "O"],
                      position=transpose([1 1 1; 2 3 4])u"nm",
                      label=[:+, :-])
    @testset ">>> Single column" begin
      crystal[:label] = [:z, :a]
      @test crystal[:label] == [:z, :a]
      crystal[3] = [:Z, :A]
      @test crystal[:label] == [:Z, :A]

      crystal[:position] = transpose([1 3 4; 2 2 6])u"nm"
      @test crystal[1, :position] == [1, 3, 4]u"nm"
      @test crystal[2, :position] == [2, 2, 6]u"nm"

      crystal[2] = transpose([6 2 8; 4 3 2])
      @test crystal[1, :position] == [6, 2, 8]
      @test crystal[2, :position] == [4, 3, 2]
    end

    @testset ">>> Multi-column" begin
      original = deepcopy(crystal)
      other = Crystal(eye(3)u"nm", species=["Ru", "Ta"],
                        position=transpose([2 4 6; 4 1 2]),
                        label=[:a, :b])
      crystal[[:species, :label]] = other[[:species, :label]]
      @test crystal[:species] === other[:species]
      @test crystal[:label] === other[:label]

      crystal[[:species, :label]] = original.atoms[[:species, :label]]
      @test crystal[:species] === original[:species]
      @test crystal[:label] === original[:label]

      crystal[[false, true]] = other[[:position]]
      @test crystal[:species] === original[:species]
      @test crystal[:label] === original[:label]
      @test crystal[:position] === other[:position]

      crystal[[false, true]] = transpose([1 1 2; 3 3 2])
      @test crystal[1, :position] == [1, 1, 2]
      @test crystal[2, :position] == [3, 3, 2]

      crystal[[:extra_column, :species]] = "Al"
      @test crystal[:extra_column] == ["Al", "Al"]
      @test crystal[:species] == ["Al", "Al"]
    end

    @testset ">>> Single-row, Single-Column" begin
      crystal[1, :label] = :aa
      @test crystal[1, :label] == :aa

      crystal[1, :position] = [1, 2, 3]
      @test crystal[1, :position] == [1, 2, 3]
      @test typeof(crystal[1, :position]) <: Position

      crystal[1, :position] = 1
      @test crystal[1, :position] == [1, 1, 1]

      crystal[1, :position] = :2
      @test crystal[1, :position] == [2, 2, 2]
      crystal[1, :position] = :2, :3, :4
      @test crystal[1, :position] == [2, 3, 4]

      @test_throws MethodError crystal[1, :position] = :a
      @test_throws ErrorException crystal[1, :nonexistent] = "Al"
    end

    @testset "Single-row, Multi-Column" begin
      crystal[1, [:species, :label]] = "aha"
      @test crystal[1, :species] == "aha"
      @test crystal[1, :label] == :aha

      crystal[2, [true, false]] = "zha"
      @test crystal[2, :species] == "zha"
      @test_throws ErrorException crystal[1, [:species, :nonexistent]] = "Al"
    end

    @testset ">>> Multi-row, single-Column" begin
      crystal = Crystal(eye(3)u"nm", species=["Al", "O", "O"],
                      position=transpose([1 1 1; 2 3 4; 4 5 2]),
                      label=[:+, :-, :0])
      crystal[[1, 3], :species] = "Ala"
      @test crystal[:species] == ["Ala", "O", "Ala"]

      crystal[[2, 3], :species] = ["H", "B"]
      @test crystal[:species] == ["Ala", "H", "B"]

      crystal[[true, false, true], 2] = transpose([1 2 3; 4 5 6])
      @test crystal[1, :position] == [1, 2, 3]
      @test crystal[2, :position] == [2, 3, 4]
      @test crystal[3, :position] == [4, 5, 6]
    end

    @testset "Multi-row, multi-Column" begin
      crystal = Crystal(eye(3)u"nm", species=["Al", "O", "O"],
                      position=transpose([1 1 1; 2 3 4; 4 5 2]),
                      label=[:+, :-, :0])
      other = Crystal(eye(3)u"nm", species=["H", "B", "C"],
                      position=transpose([2 1 2; 3 4 3; 5 2 5]),
                      label=[:a, :b, :c])
      crystal[[true, false, true], [:label, :position]] = other[1:2, [3, 2]]
      @test crystal[:label] == [:a, :-, :b]
      @test crystal[1, :position] == [2, 1, 2]
      @test crystal[2, :position] == [2, 3, 4]
      @test crystal[3, :position] == [3, 4, 3]

      crystal[[false, true, true], [:label, :species]] = ["aa", "bb"]
      @test crystal[:label] == [:a, "aa", "bb"]
      @test crystal[:species] == ["Al", "aa", "bb"]

      crystal[1:2, [:label, :species]] = "A"
      @test crystal[:label] == ["A", "A", "bb"]
      @test crystal[:species] == ["A", "A", "bb"]
    end
  end
end

@testset "> Mutating functions" begin
  crystal = Crystal(eye(3)u"nm", species=["Al", "O", "O"],
                  position=transpose([1 1 1; 2 3 4; 4 5 2]),
                  label=[:+, :-, :0])
  other = Crystal(eye(3)u"nm", other=["H", "B", "C"],
                  position=transpose([2 1 2; 3 4 3; 5 2 5]),
                  label=[:a, :b, :c])
  df = DataFrame(A=1:3, B=6:8)
  original = deepcopy(crystal)

  merge!(crystal, other.atoms, df)
  @test names(crystal) == [:species, :position, :label, :other, :A, :B]
  @test crystal[:species] == ["Al", "O", "O"]
  @test crystal[:position] === other[:position]
  @test crystal[:label] === other[:label]
  @test crystal[:A] === df[:A]
  @test crystal[:B] === df[:B]

  delete!(crystal, [:A, :B])
  @test names(crystal) == [:species, :position, :label, :other]

  deleterows!(crystal, 3)
  @test crystal[:species] == ["Al", "O"]
  @test crystal[:label] == [:a, :b]
end

@testset "> Adding atoms" begin
  crystal = Crystal(eye(3)u"nm", position=transpose([1 1 1; 2 3 4; 4 5 2]))
  push!(crystal, ([1, 4, 2], ))
  @test size(crystal, 1) == 4
  @test crystal[4, :position] == [1, 4, 2]

  # Add via tuple
  crystal = Crystal(eye(3)u"nm", species=["Al", "O", "O"],
                  position=transpose([1 1 1; 2 3 4; 4 5 2]),
                  label=[:+, :-, :0])
  push!(crystal, ("B", [1, 4, 2], :aa))
  @test size(crystal, 1) == 4
  @test crystal[4, :species] == "B"
  @test crystal[4, :position] == [1, 4, 2]
  @test crystal[4, :label] == :aa

  # Add only position, everything else NA
  push!(crystal, [5, 5, 3])
  @test size(crystal, 1) == 5
  @test crystal[5, :species] === NA
  @test crystal[5, :position] == [5, 5, 3]
  @test crystal[5, :label] === NA

  # Add a new column
  push!(crystal, [6, 6, 2], special=:special, species="BaBa")
  @test size(crystal, 1) == 6
  @test names(crystal) == [:species, :position, :label, :special]
  @test crystal[6, :species] == "BaBa"
  @test crystal[6, :position] == [6, 6, 2]
  for i in 1:5
    @test crystal[i, :special] === NA
  end
  @test crystal[6, :special] == :special
  @test crystal[6, :label] === NA

  # add a row using a dataframe
  append!(crystal, crystal[end, :])
  @test size(crystal, 1) == 7
  @test crystal[end, :position] == crystal[end - 1, :position]

  # add row using iterator
  for row in eachrow(crystal)
      push!(crystal, row)
      break
  end
  @test size(crystal, 1) == 8
  @test crystal[end, :position] == crystal[1, :position]
end

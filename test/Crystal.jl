@testset "> Construction" begin
    @testset ">> Empty 2d" begin
        crystal = Crystal(eye(2)u"nm")
        @test size(crystal.positions) == (2, 0)
        @test nrow(crystal.properties) == 0
        @test length(crystal) == 0
        @test is_fractional(crystal) == false
    end

    @testset ">> 3d with real positions" begin
        crystal = Crystal(eye(3)u"nm", tpositions=[1 1 1; 2 3 4]u"nm")
        @test length(crystal) == 2
        @test nrow(crystal.properties) == 0
        @test crystal.positions[:, 1] == [1, 1, 1]u"nm"
        @test crystal.positions[:, 2][1] == 2u"nm"
        @test crystal.positions[:, 2][2] == 3u"nm"
        @test crystal.positions[:, 2][3] == 4u"nm"
        @test is_fractional(crystal) == false
    end

    @testset ">> 4d with fractional positions" begin
        crystal = Crystal(eye(4)u"nm", tpositions=[1 1 1 1; 2 3 4 5])
        @test length(crystal) == 2
        @test nrow(crystal.properties) == 0
        @test crystal.positions[:, 1] == [1, 1, 1, 1]
        @test crystal.positions[:, 2] == [2, 3, 4, 5]
        @test is_fractional(crystal) == true
    end

    @testset ">> 2d with atomic properties" begin
        crystal = Crystal(eye(2)u"nm", tpositions=[1 1; 2 3]u"nm",
                          species=["Al", "O"])
        @test length(crystal) == 2
        @test nrow(crystal.properties) == 2
        @test crystal.properties[:species] == ["Al", "O"]
        @test crystal.positions[:, 1] == [1, 1]u"nm"
        @test crystal.positions[:, 2] == [2, 3]u"nm"
        @test is_fractional(crystal) == false
    end

    @testset ">> Convert position to crystal's type" begin
        position_for_crystal = Crystals.Structures.position_for_crystal
        cell = [0 1 1; 1 0 1; 1 1 0]u"nm"
        real = Crystal(cell)
        @test !is_fractional(real)
        @test position_for_crystal(real, [1, 1, 1]) == cell * [1, 1, 1]
        @test position_for_crystal(real, [1, 1, 1]u"nm") == [1, 1, 1]u"nm"
        frac = Crystal(Float64[0 1 1; 1 0 1; 1 1 0]u"nm", position=[0, 0, 0])
        @test position_for_crystal(frac, [1, 1, 1]) == [1, 1, 1]
        @test position_for_crystal(frac, [1, 1, 1]u"nm") == inv(frac.cell) * [1, 1, 1]u"nm"
    end

end

@testset "> Pushing" begin
    @testset "> simple line" begin
        crystal = Crystal(Float64[0 1 1; 1 0 1; 1 1 0]u"nm")

        push!(crystal, [0.25, 0.25, 0.25], species="Al")
        @test length(crystal) == 1
        @test all(crystal.positions[:, end] .== [0.5, 0.5, 0.5]u"nm")
        @test Set(names(crystal.properties)) == Set([:species])
        @test crystal.properties[:species] == ["Al"]

        push!(crystal, [0.25, 0.25, 0.25]u"nm", species="α")
        @test length(crystal) == 2
        @test all(crystal.positions[:, end] .== [0.25, 0.25, 0.25]u"nm")
        @test Set(names(crystal.properties)) == Set([:species])
        @test nrow(crystal.properties) == 2
        @test crystal.properties[end, :species] == "α"
    end
end

@testset "> Check direct indexing" begin
    @testset ">> getindex" begin
        crystal = Crystal(eye(2)u"nm", species=["Al", "O", "O"],
                          position=[1 1 1; 2 3 4], label=[:+, :-, :-])
        @testset ">>> integer" begin
            @test typeof(crystal[1]) === typeof(crystal)
            @test crystal[1].cell === crystal.cell
            @test crystal[1].positions == crystal.positions[:, 1:1]
            @test crystal[1].properties == crystal.properties[1, :]
        end

        @testset ">>> symbol" begin
            @test crystal[:position] === crystal.positions
            @test crystal[:species] === crystal.properties[:species]
        end

        @testset ">>> range" begin
            @test typeof(crystal[2:3]) === typeof(crystal)
            @test crystal[2:3].cell === crystal.cell
            @test crystal[2:3].positions == crystal.positions[:, 2:3]
            @test crystal[2:3].properties == crystal.properties[2:3, :]
        end

        @testset ">>> array of integers" begin
            @test typeof(crystal[[3, 2]]) === typeof(crystal)
            @test crystal[[3, 2]].cell === crystal.cell
            @test crystal[[3, 2]].positions == crystal.positions[:, [3, 2]]
        end

        @testset ">>> list of symbols" begin
            @test typeof(crystal[[:species]]) <: DataFrame
            @test names(crystal[[:species]]) == [:species]
            @test crystal[[:species]][:species] == crystal.properties[:species]
            @test typeof(crystal[[:species, :position]]) === typeof(crystal)
            @test names(crystal[[:species, :position]]) == [:species, :position]
            @test crystal[[:species, :position]][:species] == crystal.properties[:species]
        end

        @testset ">>> integer|range|array of integers|Colon, symbol" begin
            @test crystal[2, :position] == crystal.positions[:, 2]
            @test crystal[2:3, :species] == crystal.properties[2:3, :species]
            @test crystal[[3, 2], :position] == crystal.positions[:, [3, 2]]
            @test crystal[[3, 2], :species] == crystal.properties[[3, 2], :species]
            @test crystal[:, :species] === crystal.properties[:, :species]
            @test crystal[:, :position] === crystal.positions
        end

        @testset ">>> integer|range|array of integers|Colon, list of symbol" begin
            @test typeof(crystal[2, [:position, :label]]) === typeof(crystal)
            @test names(crystal[2, [:position, :label]]) == [:label, :position]
            @test nrow(crystal[2, [:position, :label]]) == 1
            @test crystal[2, [:position, :label]][:position] == crystal.positions[:, 2:2]
            @test crystal[2, [:position, :label]][:label] == [crystal.properties[2, :label]]

            @test typeof(crystal[2, [:species, :label]]) === typeof(crystal.properties)
            @test crystal[2, [:species, :label]][:label] == [crystal.properties[2, :label]]
            @test names(crystal[2, [:species, :label]]) == [:species, :label]
            @test nrow(crystal[2, [:species, :label]]) == 1

            subcrystal = crystal[[3, 2], [:position, :label]]
            @test typeof(subcrystal) === typeof(crystal)
            @test subcrystal[:position] == crystal.positions[:, [3, 2]]
            @test subcrystal[:label] == crystal.properties[[3, 2], :label]

            @test typeof(crystal[[3, 2], [:species, :label]]) === typeof(crystal.properties)
        end

        @testset ">>> position components" begin
            @test crystal[2, :position, 1] == crystal.positions[1, 2]
            @test crystal[:position, 1] == crystal.positions[1, :]
        end
    end

    @testset ">> setindex!" begin
        original = Crystal([0 5 5; 5 0 5; 5 5 0]u"nm", species=["Al", "O", "O"],
                           position=[1, 1, 1]u"nm",
                           position=[1, 2, 1]u"nm",
                           position=[0, 0, 1]u"nm",
                           label=[:+, :-, :-])

        @testset ">>> symbol" begin
            crystal = deepcopy(original)

            crystal[:label] = :u
            @test crystal[:label] == [:u, :u, :u]

            crystal[:label] = [:v, :w, :x]
            @test crystal[:label] == [:v, :w, :x]

            crystal[:position] = [1 0 0; 0 2 1; 4 0 1]u"nm"
            @test crystal[:position] == [1 0 0; 0 2 1; 4 0 1]u"nm"
        end

        @testset ">>> atom number, symbol" begin
            crystal = deepcopy(original)
            crystal[2, :label] = :z
            @test crystal[:label] == [:+, :z, :-]

            crystal[2, :position] = [3, 4, 5]u"nm"
            @test crystal[2, :position] == [3, 4, 5]u"nm"

            crystal[2, :position] = [1, 1, -1]
            @test crystal[2, :position] == [0, 0, 10]u"nm"
        end

        @testset ">>> multiple indices, symbol" begin
            crystal = copy(original)

            crystal[2:3, :label] = :a
            @test crystal[:label] == [:+, :a, :a]

            crystal[[1, 3], :position] = transpose([1 2 3; 4 5 6])u"nm"
            @test crystal.positions == transpose([1 2 3; 1 2 1; 4 5 6])u"nm"

            crystal[[1, 3], :position] = transpose([1 0 0; 1 -1 0])
            @test crystal.positions == transpose([0 5 5; 1 2 1; -5 5 0])u"nm"
        end

        #
        #     @testset ">>> Multi-column" begin
        #       original = deepcopy(crystal)
        #       other = Crystal(eye(3)u"nm", species=["Ru", "Ta"],
        #                         position=transpose([2 4 6; 4 1 2]),
        #                         label=[:a, :b])
        #       crystal[[:species, :label]] = other[[:species, :label]]
        #       @test crystal[:species] === other[:species]
        #       @test crystal[:label] === other[:label]
        #
        #       crystal[[:species, :label]] = original.atoms[[:species, :label]]
        #       @test crystal[:species] === original[:species]
        #       @test crystal[:label] === original[:label]
        #
        #       crystal[[false, true]] = other[[:position]]
        #       @test crystal[:species] === original[:species]
        #       @test crystal[:label] === original[:label]
        #       @test crystal[:position] === other[:position]
        #
        #       crystal[[false, true]] = transpose([1 1 2; 3 3 2])
        #       @test crystal[1, :position] == [1, 1, 2]
        #       @test crystal[2, :position] == [3, 3, 2]
        #
        #       crystal[[:extra_column, :species]] = "Al"
        #       @test crystal[:extra_column] == ["Al", "Al"]
        #       @test crystal[:species] == ["Al", "Al"]
        #     end
        #
        #     @testset ">>> Single-row, Single-Column" begin
        #       crystal[1, :label] = :aa
        #       @test crystal[1, :label] == :aa
        #
        #       crystal[1, :position] = [1, 2, 3]
        #       @test crystal[1, :position] == [1, 2, 3]
        #       @test typeof(crystal[1, :position]) <: Position
        #
        #       crystal[1, :position] = 1
        #       @test crystal[1, :position] == [1, 1, 1]
        #
        #       crystal[1, :position] = :2
        #       @test crystal[1, :position] == [2, 2, 2]
        #       crystal[1, :position] = :2, :3, :4
        #       @test crystal[1, :position] == [2, 3, 4]
        #
        #       @test_throws MethodError crystal[1, :position] = :a
        #       @test_throws ErrorException crystal[1, :nonexistent] = "Al"
        #     end
        #
        #     @testset "Single-row, Multi-Column" begin
        #       crystal[1, [:species, :label]] = "aha"
        #       @test crystal[1, :species] == "aha"
        #       @test crystal[1, :label] == :aha
        #
        #       crystal[2, [true, false]] = "zha"
        #       @test crystal[2, :species] == "zha"
        #       @test_throws ErrorException crystal[1, [:species, :nonexistent]] = "Al"
        #     end
        #
        #     @testset ">>> Multi-row, single-Column" begin
        #       crystal = Crystal(eye(3)u"nm", species=["Al", "O", "O"],
        #                       position=transpose([1 1 1; 2 3 4; 4 5 2]),
        #                       label=[:+, :-, :0])
        #       crystal[[1, 3], :species] = "Ala"
        #       @test crystal[:species] == ["Ala", "O", "Ala"]
        #
        #       crystal[[2, 3], :species] = ["H", "B"]
        #       @test crystal[:species] == ["Ala", "H", "B"]
        #
        #       crystal[[true, false, true], 2] = transpose([1 2 3; 4 5 6])
        #       @test crystal[1, :position] == [1, 2, 3]
        #       @test crystal[2, :position] == [2, 3, 4]
        #       @test crystal[3, :position] == [4, 5, 6]
        #     end
        #
        #     @testset "Multi-row, multi-Column" begin
        #       crystal = Crystal(eye(3)u"nm", species=["Al", "O", "O"],
        #                       position=transpose([1 1 1; 2 3 4; 4 5 2]),
        #                       label=[:+, :-, :0])
        #       other = Crystal(eye(3)u"nm", species=["H", "B", "C"],
        #                       position=transpose([2 1 2; 3 4 3; 5 2 5]),
        #                       label=[:a, :b, :c])
        #       crystal[[true, false, true], [:label, :position]] = other[1:2, [3, 2]]
        #       @test crystal[:label] == [:a, :-, :b]
        #       @test crystal[1, :position] == [2, 1, 2]
        #       @test crystal[2, :position] == [2, 3, 4]
        #       @test crystal[3, :position] == [3, 4, 3]
        #
        #       crystal[[false, true, true], [:label, :species]] = ["aa", "bb"]
        #       @test crystal[:label] == [:a, "aa", "bb"]
        #       @test crystal[:species] == ["Al", "aa", "bb"]
        #
        #       crystal[1:2, [:label, :species]] = "A"
        #       @test crystal[:label] == ["A", "A", "bb"]
        #       @test crystal[:species] == ["A", "A", "bb"]
        #     end
    end
end
#
# @testset "> Mutating functions" begin
#   crystal = Crystal(eye(3)u"nm", species=["Al", "O", "O"],
#                   position=transpose([1 1 1; 2 3 4; 4 5 2]),
#                   label=[:+, :-, :0])
#   other = Crystal(eye(3)u"nm", other=["H", "B", "C"],
#                   position=transpose([2 1 2; 3 4 3; 5 2 5]),
#                   label=[:a, :b, :c])
#   df = DataFrame(A=1:3, B=6:8)
#   original = deepcopy(crystal)
#
#   merge!(crystal, other.atoms, df)
#   @test names(crystal) == [:species, :position, :label, :other, :A, :B]
#   @test crystal[:species] == ["Al", "O", "O"]
#   @test crystal[:position] === other[:position]
#   @test crystal[:label] === other[:label]
#   @test crystal[:A] === df[:A]
#   @test crystal[:B] === df[:B]
#
#   delete!(crystal, [:A, :B])
#   @test names(crystal) == [:species, :position, :label, :other]
#
#   deleterows!(crystal, 3)
#   @test crystal[:species] == ["Al", "O"]
#   @test crystal[:label] == [:a, :b]
# end
#
# @testset "> Adding atoms" begin
#   crystal = Crystal(eye(3)u"nm", position=transpose([1 1 1; 2 3 4; 4 5 2]))
#   push!(crystal, ([1, 4, 2], ))
#   @test size(crystal, 1) == 4
#   @test crystal[4, :position] == [1, 4, 2]
#
#   # Add via tuple
#   crystal = Crystal(eye(3)u"nm", species=["Al", "O", "O"],
#                   position=transpose([1 1 1; 2 3 4; 4 5 2]),
#                   label=[:+, :-, :0])
#   push!(crystal, ("B", [1, 4, 2], :aa))
#   @test size(crystal, 1) == 4
#   @test crystal[4, :species] == "B"
#   @test crystal[4, :position] == [1, 4, 2]
#   @test crystal[4, :label] == :aa
#
#   # Add only position, everything else NA
#   push!(crystal, [5, 5, 3])
#   @test size(crystal, 1) == 5
#   @test crystal[5, :species] === NA
#   @test crystal[5, :position] == [5, 5, 3]
#   @test crystal[5, :label] === NA
#
#   # Add a new column
#   push!(crystal, [6, 6, 2], special=:special, species="BaBa")
#   @test size(crystal, 1) == 6
#   @test names(crystal) == [:species, :position, :label, :special]
#   @test crystal[6, :species] == "BaBa"
#   @test crystal[6, :position] == [6, 6, 2]
#   for i in 1:5
#     @test crystal[i, :special] === NA
#   end
#   @test crystal[6, :special] == :special
#   @test crystal[6, :label] === NA
#
#   # add a row using a dataframe
#   append!(crystal, crystal[end, :])
#   @test size(crystal, 1) == 7
#   @test crystal[end, :position] == crystal[end - 1, :position]
#
#   # add row using iterator
#   for row in eachrow(crystal)
#       push!(crystal, row)
#       break
#   end
#   @test size(crystal, 1) == 8
#   @test crystal[end, :position] == crystal[1, :position]
# end

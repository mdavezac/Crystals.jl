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
        realpos = Crystal(cell)
        @test !is_fractional(realpos)
        @test position_for_crystal(realpos, [1, 1, 1]) == cell * [1, 1, 1]
        @test position_for_crystal(realpos, [1, 1, 1]u"nm") == [1, 1, 1]u"nm"
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
            @test crystal[:x] == crystal.positions[1, :]
            @test crystal[:y] == crystal.positions[2, :]
            @test crystal[:fractional] === crystal.positions
            @test crystal[:cartesian] == crystal.cell * crystal.positions
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

            cartesian = crystal[[:species, :cartesian]]
            @test !is_fractional(cartesian)
            @test cartesian[:species] == crystal.properties[:species]
            @test cartesian[:position] == crystal.cell * crystal.positions

            @test crystal[[:species, :y, :x]][:x] == crystal.positions[1, :]
            @test crystal[[:species, :y, :x]][:y] == crystal.positions[2, :]

            fractional = crystal[[:species, :fractional]]
            @test is_fractional(fractional)
            @test fractional[:species] == crystal.properties[:species]
            @test fractional[:position] == crystal.positions
        end

        @testset ">>> integer|range|array of integers|Colon, symbol" begin
            @test crystal[2, :position] == crystal.positions[:, 2]
            @test crystal[2, :fractional] == crystal.positions[:, 2]
            @test crystal[2, :cartesian] == crystal.cell * crystal.positions[:, 2]
            @test crystal[2:3, :species] == crystal.properties[2:3, :species]
            @test crystal[[3, 2], :position] == crystal.positions[:, [3, 2]]
            @test crystal[[3, 2], :species] == crystal.properties[[3, 2], :species]
            @test crystal[[3, 2], :x] == crystal.positions[1, [3, 2]]
            @test crystal[[3, 2], :x] == crystal.positions[1, [3, 2]]
            @test crystal[:, :species] === crystal.properties[:, :species]
            @test crystal[:, :position] === crystal.positions
            @test crystal[:, :fractional] === crystal.positions
            @test crystal[:, :cartesian] == crystal.cell * crystal.positions
        end

        @testset ">>> integer|range|array of integers|Colon, list of symbol" begin
            @test typeof(crystal[2, [:position, :label]]) === typeof(crystal)
            @test names(crystal[2, [:position, :label]]) == [:label, :position]
            @test nrow(crystal[2, [:position, :label]]) == 1
            @test crystal[2, [:position, :label]][:position] == crystal.positions[:, 2:2]
            @test crystal[2, [:position, :label]][:label] == [crystal.properties[2, :label]]
            @test crystal[2, [:x, :y, :label]][:x] == [crystal.positions[1, 2]]
            @test crystal[2, [:x, :y, :label]][:y] == [crystal.positions[2, 2]]

            @test typeof(crystal[2, [:species, :label]]) === typeof(crystal.properties)
            @test crystal[2, [:species, :label]][:label] == [crystal.properties[2, :label]]
            @test names(crystal[2, [:species, :label]]) == [:species, :label]
            @test nrow(crystal[2, [:species, :label]]) == 1

            subcrystal = crystal[[3, 2], [:position, :label]]
            @test typeof(subcrystal) === typeof(crystal)
            @test subcrystal[:position] == crystal.positions[:, [3, 2]]
            @test subcrystal[:label] == crystal.properties[[3, 2], :label]

            @test typeof(crystal[[3, 2], [:species, :label]]) === typeof(crystal.properties)

            subcrystal = crystal[[3, 2], [:fractional, :label]]
            @test typeof(subcrystal) === typeof(crystal)
            @test subcrystal[:position] == crystal.positions[:, [3, 2]]

            subcrystal = crystal[[3, 2], [:cartesian, :label]]
            @test !is_fractional(subcrystal)
            @test subcrystal[:position] == crystal.cell * crystal.positions[:, [3, 2]]
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
            crystal = deepcopy(original)

            crystal[2:3, :label] = :a
            @test crystal[:label] == [:+, :a, :a]

            crystal[[1, 3], :position] = transpose([1 2 3; 4 5 6])u"nm"
            @test crystal.positions == transpose([1 2 3; 1 2 1; 4 5 6])u"nm"

            crystal[[1, 3], :position] = transpose([1 0 0; 1 -1 0])
            @test crystal.positions == transpose([0 5 5; 1 2 1; -5 5 0])u"nm"
        end

        @testset ">>> From dataframe" begin
            crystal = deepcopy(original)

            df = DataFrame(species=["I", "N", "B"], μ=[1, 2, 4])
            crystal[[:species, :μ]] = df
            @test :μ ∈ names(crystal)
            @test crystal[:μ] === df[:μ]
            @test crystal[:species] === df[:species]

            df = DataFrame(species=["A", "B"], μ=[5, 3])
            crystal[[3, 1], [:species, :μ]] = df
            @test crystal[:μ] == [3, 2, 5]
            @test crystal[:species] == ["B", "N", "A"]
        end

        @testset ">>> From crystal" begin
            crystal = deepcopy(original)
            other = deepcopy(original)
            other[:label] = [:a, :b, :c]

            crystal[[:label]] = other
            # follows same logic as for dataframes
            @test crystal[:label] === other[:species]

            crystal = deepcopy(original)
            crystal[[:position, :label]] = other
            # except for positions, which are always correctly extracted
            @test crystal[:position] ≈ other[:position]

            # lets try with a fractional crystal
            crystal = Crystal(1.0 * original.cell, 1.0 * original.positions,
                              deepcopy(original.properties))
            other = Crystal([0 5 5; 5 0 5; 5 5 0]u"nm", species=["Al", "O", "O"],
                            position=[0.2  0.4   0.2; 0.1  0.0   0.1; 0.2  0.4  -0.2],
                            label=[:a, :b, :c])

            crystal[[:position, :label]] = other
            @test crystal[:label] === other[:species]
            @test crystal[:position] ≈ other.cell * other[:position]

            # try with specific rows
            crystal = Crystal(1.0 * original.cell, 1.0 * original.positions,
                              deepcopy(original.properties))
            crystal[[3, 1], [:position]] = other[[1, 3], :]
            @test crystal[[1, 3], :position] ≈ other.cell * other[[3, 1], :position]
            @test crystal[2, :position] ≈ original[2, :position]
        end

        @testset ">>> Set position components" begin
            crystal = deepcopy(original)
            crystal[1, :position, 1] = 8u"nm"
            @test crystal[1, :position, 1] ≈ 8u"nm"

            crystal[2, :position, [3, 1]] = [7, 8]u"nm"
            @test crystal[2, :position, [1, 3]] ≈ [8, 7]u"nm"

            crystal[2, :position, [3, 1]] = [7, 8]u"nm"
            @test crystal[2, :position, [1, 3]] ≈ [8, 7]u"nm"
        end
    end
end

@testset "> Deleting stuff" begin
    original = Crystal([0 5 5; 5 0 5; 5 5 0]u"nm", species=["Al", "O", "O"],
                       position=[1, 1, 1]u"nm",
                       position=[1, 2, 1]u"nm",
                       position=[0, 0, 1]u"nm",
                       label=[:+, :-, :-])
    @testset ">> columns" begin
        crystal = deepcopy(original)
        delete!(crystal, :label)
        @test :label ∉ names(crystal)
        @test :species ∈ names(crystal)

        crystal = deepcopy(original)
        delete!(crystal, [:label, :species])
        @test :label ∉ names(crystal)
        @test :species ∉ names(crystal)

        @test_throws ErrorException delete!(crystal, :position)

        crystal = deepcopy(original)
        empty!(crystal)
        @test nrow(crystal) == 0
        @test ncol(crystal) == 1
    end

    @testset ">> rows" begin
        crystal = deepcopy(original)
        deleterows!(crystal, 1)
        @test nrow(crystal) == nrow(original) - 1
        @test crystal.properties == original.properties[2:end, :]
        @test crystal.positions == original.positions[:, 2:end]

        crystal = deepcopy(original)
        deleterows!(crystal, 1:2:3)
        @test nrow(crystal) == 1
        @test crystal.properties == original.properties[2, :]
        @test crystal.positions == original.positions[:, 2:2]

        crystal = deepcopy(original)
        deleterows!(crystal, :)
        @test nrow(crystal) == 0
        @test ncol(crystal) == ncol(original)

        crystal = deepcopy(original)
        delete!(crystal, :)
        @test nrow(crystal) == 0
        @test ncol(crystal) == ncol(original)
    end
end

@testset "> vcat and append" begin
    crysA = Crystal([0 5 5; 5 0 5; 5 5 0]u"nm", species=["Al", "O", "O"],
                    position=[1, 1, 1]u"nm",
                    position=[1, 2, 1]u"nm",
                    position=[0, 0, 1]u"nm",
                    label=[:+, :-, :-])
    crysB = Crystal([0 5 5; 5 0 5; 5 5 0]u"nm",
                    species=[missing, missing],
                    position=[2, -1, 1]u"nm",
                    position=[0, 0, -1]u"nm",
                    label=[:-, :a])

    crystal = vcat(crysA, crysB)
    @test nrow(crystal) == 5
    @test crystal.positions == transpose([1 1 1; 1 2 1; 0 0 1; 2 -1 1; 0 0 -1])u"nm"
    @test crystal[1:3, :species] == ["Al", "O", "O"]
    @test all(ismissing.(crystal[4:5, :species]))
    @test crystal[:label] == [:+, :-, :-, :-, :a]
    @test_throws ArgumentError vcat(crysA, crysB[[:position, :label]])

    crystal = deepcopy(crysA[[:position, :label]])
    append!(crystal, crysB[[:position, :label]])
    @test nrow(crystal) == 5
    @test crystal.positions == transpose([1 1 1; 1 2 1; 0 0 1; 2 -1 1; 0 0 -1])u"nm"
    @test crystal[:label] == [:+, :-, :-, :-, :a]

    crystal = deepcopy(crysA)
    # eltypes do not match (one is Missing, the other String)
    @test_throws ErrorException append!(crystal, crysB)
    # mismatching columns
    @test_throws ErrorException append!(crystal, crysB[[:position, :label]])
end

@testset "> round crystal cell and positions" begin
    crystal = Crystal([0 0.501 0.501; 0.496 0.001 0.497; 0.497 0.497 0]u"nm",
                      position=[0.001, -0.001, -0.001]u"nm",
                      position=[0.25, 0.251, -0.247]u"nm")
    round!(crystal, 2)
    @test crystal.cell == [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"
    @test crystal[:position] == transpose([0 0 0; 0.25 0.25 -0.25]u"nm")
end

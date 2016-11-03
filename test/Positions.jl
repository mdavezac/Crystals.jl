using DataFrames: DataArray, @data, NA
@testset "> Construction" begin
    @test typeof(Position(1, 1)) <: Position
end

@testset "> Convertion vector <--> position" begin
   array = [1, 2, 3]
   position = Position(array)
   @test typeof(position) === Position{eltype(array), 3}
   @test [position...] == array
   @test typeof(convert(Array, position)) === typeof(array)
   @test convert(Array, position) == array
   position = convert(Position, array)
   @test eltype(position) === eltype(array)

   position = Position{Int8}(array)
   @test eltype(position) === Int8
   position = convert(Position{Int8}, array)
   @test eltype(position) === Int8
end

@testset "> Convertions matrix <--> array of positions" begin
  @testset ">> Automatic conversion" begin
    matrix = transpose([1 2 3; 4 5 6])
    positions = convert(PositionArray, matrix)

    @test eltype(positions) === Position{eltype(matrix), 3}
    @test typeof(positions) === Vector{Position{eltype(matrix), 3}}
    @test length(positions) == size(matrix, 2)
    @test length(positions[1]) == size(matrix, 1)
    for i = 1:length(positions)
      @test positions[i] == matrix[:, i]
    end

    back = convert(Array, positions)
    @test typeof(back) === typeof(matrix)
    @test size(back) == size(matrix)
    @test back ≈ matrix

    to_integer = convert(PositionArray, [1 2; 2 3])
    @test eltype(eltype(to_integer)) <: Integer
    to_float = convert(PositionArray, [1 2.; 2 3])
    @test eltype(eltype(to_float)) <: AbstractFloat
  end

  @testset ">> Explicit element type conversion" begin
    positions = convert(PositionArray{Int8}, [1 2; 3 4])
    @test eltype(eltype(positions)) === Int8
  end
end

@testset "> Convertions dataarray <--> array of positions" begin
  @testset ">> Automatic conversion" begin
    matrix = transpose([1 2 3; 4 5 6])
    positions = convert(PositionDataArray, matrix)

    @test eltype(positions) === Position{eltype(matrix), 3}
    @test typeof(positions) === DataArray{Position{eltype(matrix), 3}, 1}
    @test length(positions) == size(matrix, 2)
    @test length(positions[1]) == size(matrix, 1)
    for i = 1:length(positions)
      @test positions[i] == matrix[:, i]
    end

    back = convert(Array, positions)
    @test typeof(back) === typeof(matrix)
    @test size(back) == size(matrix)
    @test back ≈ matrix

    to_integer = convert(PositionArray, [1 2; 2 3])
    @test eltype(eltype(to_integer)) <: Integer
    to_float = convert(PositionArray, [1 2.; 2 3])
    @test eltype(eltype(to_float)) <: AbstractFloat
  end
end

@testset "> Matrix - Position(Data)Array multiplication" begin
   cell = [-1 1 1; 1 -1 1; 1 1 -1]
   positions = @data [Position(1, 1, 2), Position(3, 4, 5)]

   const expected = cell * convert(Array, positions)
   @test convert(Array, cell * positions) == expected
   @test convert(Array, cell * positions.data) == expected

   positions[1] = NA
   @test (cell * positions)[1] === NA
   @test (cell * positions)[2] !== NA
end

@testset "Sum/differences" begin
    @test Position(1, 2) + 1 == Position(2, 3)
    @test Position(1, 2) + 1.5 == Position(2.5, 3.5)

    @test Position(1, 2) + [-1, 2] == Position(0, 4)
    @test Position(1, 2) + Position(-1, 2) == Position(0, 4)
    @test Position(1, 2) + [-1.5, 2.5] == Position(-0.5, 4.5)
    @test Position(1, 2) + Position(-1.5, 2.5) == Position(-0.5, 4.5)

    @test Position(1, 2) - Position(-1.5, 2.5) == Position(2.5, -0.5)
    @test Position(1, 2) - 1 == Position(0, 1)

    array = convert(PositionArray, [1 2 3; 4 5 6])
    @test array .+ Position(1, 2) == convert(PositionArray, [2 3 4; 6 7 8])
    @test array .+ Position(1.5, 2) ==
        convert(PositionArray, [2.5 3.5 4.5; 6 7 8])
    @test array .+ [1.5, 2] == convert(PositionArray, [2.5 3.5 4.5; 6 7 8])

    @test array + array == convert(PositionArray, [2 4 6; 8 10 12])
    @test array + [1 2 3; 4 5 6] == convert(PositionArray, [2 4 6; 8 10 12])

    array = convert(PositionDataArray, [1 2 3; 4 5 6])
    @test array .+ Position(1, 2) == convert(PositionDataArray, [2 3 4; 6 7 8])
    @test array + Position(1, 2) == convert(PositionDataArray, [2 3 4; 6 7 8])
    @test array .+ Position(1.5, 2) ==
        convert(PositionDataArray, [2.5 3.5 4.5; 6 7 8])
    @test array + [1.5, 2] == convert(PositionDataArray, [2.5 3.5 4.5; 6 7 8])
    @test typeof(array + [1.5, 2]) <: PositionDataArray{Float64, 2}

    @test array + array == convert(PositionDataArray, [2 4 6; 8 10 12])
    @test array + [1 2 3; 4 5 6] ==
        convert(PositionDataArray, [2 4 6; 8 10 12])
end

@testset "Positions with units" begin
    position = Position(1, 2)u"m"
    @test position[1] == 1u"m"
    @test position[2] == 2u"m"
end

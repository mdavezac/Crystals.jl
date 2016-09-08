facts("Construction") do
    @fact typeof(Position(1, 1)) --> is_subtype(Position)
end

facts("Convertion vector <--> position") do
   array = [1, 2, 3]
   position = Position(array)
   @fact typeof(position) --> exactly(Position{eltype(array), 3})
   @fact [position...] --> array
   @fact typeof(convert(Array, position)) --> typeof(array)
   @fact convert(Array, position) --> array
   position = convert(Position, array)
   @fact eltype(position) --> eltype(array)

   position = Position{Int8}(array)
   @fact eltype(position) --> Int8
   position = convert(Position{Int8}, array)
   @fact eltype(position) --> Int8
end

facts("Convertions matrix <--> array of positions") do
  context("Automatic conversion") do
    matrix = transpose([1 2 3; 4 5 6])
    positions = convert(PositionArray, matrix)

    @fact eltype(positions) --> exactly(Position{eltype(matrix), 3})
    @fact typeof(positions) -->
    exactly(Vector{Position{eltype(matrix), 3}})
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
    positions = convert(PositionArray{Int8}, [1 2; 3 4])
    @fact eltype(eltype(positions)) --> Int8
  end
end

facts("Convertions dataarray <--> array of positions") do
  context("Automatic conversion") do
    matrix = transpose([1 2 3; 4 5 6])
    positions = convert(PositionDataArray, matrix)

    @fact eltype(positions) --> exactly(Position{eltype(matrix), 3})
    @fact typeof(positions) -->
        exactly(DataArray{Position{eltype(matrix), 3}, 1})
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

facts("Matrix - Position(Data)Array multiplication") do
   cell = [-1 1 1; 1 -1 1; 1 1 -1]
   positions = @data [Position(1, 1, 2), Position(3, 4, 5)]

   const expected = cell * convert(Array, positions)
   @fact convert(Array, cell * positions) --> expected
   @fact convert(Array, cell * positions.data) --> expected

   positions[1] = NA
   @fact (cell * positions)[1] --> exactly(NA)
   @fact (cell * positions)[2] --> not(exactly(NA))
end

facts("Sum/differences") do
    @fact Position(1, 2) + 1 --> Position(2, 3)
    @fact Position(1, 2) + 1.5 --> Position(2.5, 3.5)

    @fact Position(1, 2) + [-1, 2] --> Position(0, 4)
    @fact Position(1, 2) + Position(-1, 2) --> Position(0, 4)
    @fact Position(1, 2) + [-1.5, 2.5] --> Position(-0.5, 4.5)
    @fact Position(1, 2) + Position(-1.5, 2.5) --> Position(-0.5, 4.5)

    @fact Position(1, 2) - Position(-1.5, 2.5) --> Position(2.5, -0.5)
    @fact Position(1, 2) - 1 --> Position(0, 1)

    array = convert(PositionArray, [1 2 3; 4 5 6])
    @fact array .+ Position(1, 2) --> convert(PositionArray, [2 3 4; 6 7 8])
    @fact array .+ Position(1.5, 2) -->
        convert(PositionArray, [2.5 3.5 4.5; 6 7 8])
    @fact array .+ [1.5, 2] --> convert(PositionArray, [2.5 3.5 4.5; 6 7 8])

    @fact array + array --> convert(PositionArray, [2 4 6; 8 10 12])
    @fact array + [1 2 3; 4 5 6] --> convert(PositionArray, [2 4 6; 8 10 12])

    array = convert(PositionDataArray, [1 2 3; 4 5 6])
    @fact array .+ Position(1, 2) --> convert(PositionDataArray, [2 3 4; 6 7 8])
    @fact array + Position(1, 2) --> convert(PositionDataArray, [2 3 4; 6 7 8])
    @fact array .+ Position(1.5, 2) -->
        convert(PositionDataArray, [2.5 3.5 4.5; 6 7 8])
    @fact array + [1.5, 2] --> convert(PositionDataArray, [2.5 3.5 4.5; 6 7 8])

    @fact array + array --> convert(PositionDataArray, [2 4 6; 8 10 12])
    @fact array + [1 2 3; 4 5 6] -->
        convert(PositionDataArray, [2 4 6; 8 10 12])
end

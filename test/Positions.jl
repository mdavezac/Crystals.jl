facts("Convertion vector <--> position") do
   array = [1, 2, 3]
   position = convert(Position, array)
   @fact typeof(position) --> exactly(Position3D{eltype(array)})
   @fact [position...] --> array
   @fact typeof(convert(Array, position)) --> typeof(array)
   @fact convert(Array, position) --> array

   position = convert(Position{Int8}, array)
   @fact eltype(position) --> Int8
end

facts("Convertions matrix <--> array of positions") do
  context("Automatic conversion") do
    matrix = transpose([1 2 3; 4 5 6])
    positions = convert(PositionArray, matrix)

    @fact eltype(positions) --> exactly(Position3D{eltype(matrix)})
    @fact typeof(positions) -->
    exactly(Vector{Position3D{eltype(matrix)}})
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

    @fact eltype(positions) --> exactly(Position3D{eltype(matrix)})
    @fact typeof(positions) -->
        exactly(DataArray{Position3D{eltype(matrix)}, 1})
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

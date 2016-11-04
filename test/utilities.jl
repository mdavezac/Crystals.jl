@testset "> Periodic images" begin
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]

  for i in 1:10
    position = cell * rand(3)
    image = position + cell * rand(Int8, (3, ))
    @test is_periodic(position, image, cell)
    @test !is_periodic(position, image + rand(3) * 1e-4, cell)

    positions = cell * rand((3, 10))
    images = positions + cell * rand(Int8, (3, 10))
    positions[:, [1, 3, 4]] += rand((3, 3))
    @test is_periodic(positions, images, cell)[2]
    @test all(is_periodic(positions, images, cell)[5:end])
    @test !any(is_periodic(positions, images, cell)[[1, 3, 4]])
  end
end

# @testset "> Fold back into cell" begin
#   cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]
#
#   for i in 1:10
#     position = cell * (rand(3) + rand(Int8, (3, )))
#     folded = into_cell(position, cell)
#     @test is_periodic(position, folded, cell) --> true
#     @test 0 .≤ inv(cell) * folded .< 1 --> all
#
#     position = cell * (rand((3, 10)) + rand(Int8, (3, 10)))
#     folded = into_cell(position, cell)
#     @test is_periodic(position, folded, cell) --> all
#     @test 0 .≤ inv(cell) * folded .< 1 --> all
#   end
# end
#
# @testset "> Fold back around origin" begin
#   cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]
#
#   for i in 1:10
#     position = cell * (rand(3) + rand(Int8, (3, )))
#     folded = origin_centered(position, cell)
#     @test is_periodic(position, folded, cell) --> true
#     @test -0.5 .< inv(cell) * folded .≤ 0.5 --> all
#
#     position = cell * (rand((3, 10)) + rand(Int8, (3, 10)))
#     folded = origin_centered(position, cell)
#     @test is_periodic(position, folded, cell) --> all
#     @test -0.5 .< inv(cell) * folded .≤ 0.5 --> all
#   end
# end
#
# function none_smaller(position::Vector, cell::Matrix)
#   const d = norm(position)
#   for i = -2:2, j = -2:2, k = -2:2
#     i == j == k == 0 && continue
#     norm(position + cell * [i, j, k]) < d - 1e-12  && return false
#   end
#   true
# end
# none_smaller(cell::Matrix) = x -> none_smaller(x, cell)
#
# @testset "> Fold back into voronoi" begin
#   cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]
#
#   for i in 1:20
#     position = cell * (rand(3) + rand(Int8, (3, )))
#     folded = into_voronoi(position, cell)
#     @test is_periodic(position, folded, cell) --> true
#     @test folded --> none_smaller(cell)
#
#     position = cell * (rand((3, 10)) + rand(Int8, (3, 10)))
#     folded = into_voronoi(position, cell)
#     @test is_periodic(position, folded, cell) --> all
#     for i = 1:size(position, 2)
#       @test folded[:, i] --> none_smaller(cell)
#     end
#   end
# end
#
# @testset "> Supercell" begin
#   lattice = Crystal(
#     [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0], 2.0,
#     position=[0 0.25; 0 0.25; 0 0.25], species=["In", "Ga"]
#   )
#
#   result = supercell(lattice, lattice.cell * [-1 1 1; 1 -1 1; 1 1 -1])
#   @test result.cell --> roughly(eye(3))
#   @test result.scale --> roughly(lattice.scale)
#   @test nrow(result) % nrow(lattice) --> 0
#   const ncells =  nrow(result) // nrow(lattice)
#   @test ncells --> 4
#   @test names(result) --> ∪(names(result), (:site_id, :cell_id))
#   @test 1 .≤ result[:site_id] .≤ nrow(lattice) --> all
#   actual = convert(Array, result[:position])
#   expected = convert(Array, lattice[result[:site_id], :position])
#   @test is_periodic(actual, expected, lattice.cell) --> all
#   @test length(unique(result[:cell_id])) --> ncells
#   for i in unique(result[:cell_id])
#     @test countnz(result[:cell_id] .== i) --> nrow(lattice)
#   end
#   @test length(unique(result[:site_id])) --> nrow(lattice)
#   for i in unique(result[:site_id])
#     @test countnz(result[:site_id] .== i) --> ncells
#   end
#   @test length(unique(result[:position])) --> nrow(result)
# end
#
# @testset "> cell parameters" begin
#     cells = Matrix[
#       [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0],
#     ]
#     parameters = Vector[
#       [√0.5, √0.5, √0.5, 60, 60, 60]
#     ]
#     for (cell, params) in zip(cells, parameters)
#         actual_params = cell_parameters°(cell)
#         @test [actual_params...] --> roughly([params...])
#         @test [cell_parameters°(cell_parameters°(actual_params...))...] -->
#                 roughly([params...])
#         @test det(inv(cell) * cell_parameters°(actual_params...)) --> roughly(1)
#     end
# end

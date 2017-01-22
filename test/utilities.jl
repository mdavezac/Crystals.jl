using Unitful: ustrip, dimension
using Crystals.Utilities: eldimension, CellParameters

@testset "> Hart-Forcade" begin
    lattice = Crystal([0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm", position=[0, 0, 0]u"nm")
    crystal = Crystal([1 0 0; 0 1 0; 0 0 1]u"nm",
                      position=[0, 0, 0]u"nm",
                      position=[0.5, 0.5, 0]u"nm",
                      position=[0.5, 0, 0.5]u"nm",
                      position=[0, 0.5, 0.5]u"nm")
    hf = @inferred(hart_forcade(lattice.cell, crystal.cell))
    @test hf.quotient == [1, 2, 2]
    for i in eachindex(crystal)
        @test isinteger(hf.transform * crystal[i, :position])
    end
    for i=1:hf.quotient[1], j=1:hf.quotient[2], k=1:hf.quotient[3]
        @test isinteger(inv(lattice.cell) * (inv(hf.transform) * [i, j, k]))
    end
end

@testset "> eldimension" begin
    @test @inferred(eldimension(1u"m")) === typeof(dimension(u"m"))
    @test @inferred(eldimension(1)) === typeof(dimension(1))
    @test @inferred(eldimension([1, 1]u"m")) === typeof(dimension(u"m"))
    @test @inferred(eldimension([1, 1])) === typeof(dimension(1))
end

# @testset "> is_position_array" begin
#     @test !is_position_array([1, 1])
#     @test !is_position_array(Position(1, 1))
#     @test is_position_array([Position(1, 1)])
#     @test is_position_array([1 1; 1 1])
# end
#
@testset "> to fractional" begin
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"
  position = [0.1, 0.1, 0.2]u"nm"
  @test @inferred(to_fractional(position, cell)) ≈ inv(cell) * position
  uposition = ustrip(position)
  @test @inferred(to_fractional(uposition, cell)) === uposition
end

@testset "> to cartesian" begin
    cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"
    position = [0.1, 0.1, 0.2]u"nm"
    @test @inferred(to_cartesian(position, cell)) === position
    uposition = ustrip(position)
    @test @inferred(to_cartesian(uposition, cell)) ≈ cell * uposition
end

@testset "> Periodic images" begin
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"

  for i in 1:1 #10
    position = cell * rand(3)
    image = position + cell * rand(Int8, (3, ))
    @test @inferred(is_periodic(position, image, cell))
    @test !@inferred(is_periodic(position, image + (rand(3) * 1e-4) * u"nm", cell))

    positions = cell * rand((3, 10))
    images = positions + cell * rand(Int8, (3, 10))
    positions[:, [1, 3, 4]] += rand((3, 3))u"nm"
    @test @inferred(is_periodic(positions, images, cell))[2]
    @test all(is_periodic(positions, images, cell)[5:end])
    @test !any(is_periodic(positions, images, cell)[[1, 3, 4]])
  end
end

@testset "> Fold back into cell" begin
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"

  for i in 1:10
    position = cell * (rand(3) + rand(Int8, (3, )))
    folded = into_cell(position, cell)
    @test @inferred(is_periodic(position, folded, cell))
    @test all(0 .≤ inv(cell) * folded .< 1)

    position = cell * (rand((3, 10)) + rand(Int8, (3, 10)))
    folded = into_cell(position, cell)
    @test all(@inferred(is_periodic(position, folded, cell)))
    @test all(0 .≤ inv(cell) * folded .< 1)
  end
end

@testset "> Fold back around origin" begin
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"

  for i in 1:10
    position = cell * (rand(3) + rand(Int8, (3, )))
    folded = @inferred(origin_centered(position, cell))
    @test is_periodic(position, folded, cell)
    @test all(-0.5 .< inv(cell) * folded .≤ 0.5)

    position = cell * (rand((3, 10)) + rand(Int8, (3, 10)))
    folded = @inferred(origin_centered(position, cell))
    @test all(is_periodic(position, folded, cell))
    @test all(-0.5 .< inv(cell) * folded .≤ 0.5)
  end
end

function none_smaller(position::Vector, cell::Matrix)
  const d = norm(position)
  for i = -2:2, j = -2:2, k = -2:2
    i == j == k == 0 && continue
    norm(position + cell * [i, j, k]) < d * (1 - 1e-12)  && return false
  end
  true
end

@testset "> Fold back into voronoi" begin
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"

  for i in 1:10
    position = cell * (rand(3) + rand(Int8, (3, )))
    folded = into_voronoi(position, cell)
    @test is_periodic(position, folded, cell)
    @test none_smaller(folded, cell)

    position = cell * (rand((3, 10)) + rand(Int8, (3, 10)))
    folded = into_voronoi(position, cell)
    @test all(is_periodic(position, folded, cell))
    for i = 1:size(position, 2)
      @test none_smaller(folded[:, i], cell)
    end
  end
end

@testset "> Supercell" begin
  lattice = Crystal(
    [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm",
    position=[0 0.25; 0 0.25; 0 0.25]u"nm", species=["In", "Ga"]
  )

  result = supercell(lattice, lattice.cell * [-1 1 1; 1 -1 1; 1 1 -1])
  @test result.cell ≈ eye(3)u"nm"
  @test nrow(result) % nrow(lattice) == 0
  const ncells =  nrow(result) // nrow(lattice)
  @test ncells == 4
  @test names(result) == ∪(names(result), (:site_id, ))
  @test all(1 .≤ result[:site_id] .≤ nrow(lattice))
  actual = convert(Array, result[:position])
  expected = convert(Array, lattice[result[:site_id], :position])
  @test all(is_periodic(actual, expected, lattice.cell))
  @test length(unique(result[:site_id])) == nrow(lattice)
  for i in unique(result[:site_id])
    @test countnz(result[:site_id] .== i) == ncells
  end
  for atom in eachatom(result)
      @test countnz(is_periodic(atom[:cartesian], lattice[:cartesian], lattice.cell)) == 1
      @test countnz(is_periodic(atom[:cartesian], result[:cartesian], result.cell)) == 1
  end
end

@testset "> cell parameters" begin
    cells = Matrix[
        [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm",
    ]
    parameters = CellParameters[
      CellParameters(sqrt(0.5)u"nm", sqrt(0.5)u"nm", sqrt(0.5)u"nm", 60u"°", 60u"°", 60u"°")
    ]
    for (cell, params) in zip(cells, parameters)
        actual_params = cell_parameters(cell)
        @test actual_params.a ≈ params.a
        @test actual_params.b ≈ params.b
        @test actual_params.c ≈ params.c
        @test actual_params.α ≈ params.α
        @test actual_params.β ≈ params.β
        @test actual_params.γ ≈ params.γ

        round_trip = cell_parameters(cell_parameters(actual_params))
        @test round_trip.a ≈ params.a
        @test round_trip.b ≈ params.b
        @test round_trip.c ≈ params.c
        @test round_trip.α ≈ params.α
        @test round_trip.β ≈ params.β
        @test round_trip.γ ≈ params.γ
        @test det(inv(cell) * cell_parameters(actual_params)) ≈ 1
    end
end

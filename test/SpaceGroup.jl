cells = Any[
  [1.5 -0.5 -0.5;  0  1 -0.5; 0.5 -0.5  1]u"nm",
  [1.5 -0.5 -0.5;  0  1 -0.5; 0.5 -0.5  1]
]
@testset "> Potential equivalents regression $(unit(eltype(cell)))" for cell in cells
    cell = [1.5 -0.5 -0.5;  0  1 -0.5; 0.5 -0.5  1]u"nm"
    avecs, bvecs, cvecs = Crystals.SpaceGroup.potential_equivalents(cell)
    norms = mapreducedim(x -> x*x, +, cell, 1)
    @test size(avecs) == (3, 4)
    @test all(mapreducedim(x -> x * x, +, avecs, 1) .== norms[1])

    @test size(bvecs) == (3, 8)
    @test all(mapreducedim(x -> x * x, +, bvecs, 1) .== norms[2])

    @test size(cvecs) == (3, 8)
    @test all(mapreducedim(x -> x * x, +, cvecs, 1) .== norms[3])
end

@testset "> Point group operations" begin
    parameters = [
        ([0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm", 48),
        ([-0.5 0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 -0.5], 48),
        ([-0.6 0.5 0.5; 0.6 -0.5 0.5; 0.6 0.5 -0.5], 4),
        ([-0.7 0.7 0.7; 0.6 -0.5 0.5; 0.6 0.5 -0.5], 8),
        ([-0.765 0.7 0.7; 0.665 -0.5 0.5; 0.6 0.5 -0.5], 2)
    ]
    allops = point_group_operations([0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0])

    @testset "Cell: $cell" for (cell, numops) in parameters
        ops = point_group_operations(cell)
        @test length(ops) == numops
        for op in ops
            @test size(op.scalefwd) == (3, 3)
            @test op.offset ≈ [0, 0, 0]
            transformation = inv(cell) * op.scalefwd * cell
            @test isinteger(round(transformation, 8))
            @test volume(transformation) ≈ 1.0

            if numops != 48
                found = 0
                for op in allops
                    transformation = inv(cell) * op.scalefwd * cell
                    isinteger(round(transformation, 8)) || continue
                    isapprox(volume(transformation), 1, rtol=1e-8) || continue
                    found += 1
                end
                @test found == numops
            end
        end
    end
end

# @testset "> Inner translations" begin do
#     diamond = Lattices.diamond()
#     @testset ">> No translations" begin
#         translations = inner_translations(diamond)
#         @test length(translations) --> 0
#     end
#     @testset ">> Supercell" begin
#         cell = -1
#         while det(cell) ≤ 0
#           cell = rand(-3:3, (3, 3))
#         end
#         large = supercell(diamond, diamond.cell * cell)
#         @test nrow(large) --> round(Integer, det(cell)) * nrow(diamond)
#         translations = inner_translations(large)
#         @test length(translations) --> round(Integer, det(cell) - 1)
#     end
# end

# @testset "> Make primitive" begin do
#     zinc_blende = Lattices.zinc_blende()
#     @test is_primitive(zinc_blende) --> true
#     @test primitive(zinc_blende) --> exactly(zinc_blende)

#     cell = -1
#     while det(cell) ≤ 0
#         cell = rand(-3:3, (3, 3))
#     end
#     large = supercell(zinc_blende, zinc_blende.cell * cell)
#     @test is_primitive(large) --> false
#     prim = primitive(large)
#     @test niggly(prim.cell) --> roughly(niggly(zinc_blende.cell))
#     @test nrow(prim) --> nrow(zinc_blende)
#     for i in 1:2
#         prim_pos = prim[i, :position]
#         d_pos = zinc_blende[i, :position]
#         @test is_periodic(prim_pos, d_pos, prim.cell) --> true
#         @test is_periodic(prim_pos, d_pos, zinc_blende.cell) --> true
#     end
# end

# @testset "> Space group operations" begin do
#     @testset ">> fcc" begin
#       fcc = Lattices.fcc()
#       ops = space_group(fcc)
#       @test length(ops) --> 48
#       for op in ops
#           @test op * [0, 0, 0] --> [0, 0, 0]
#       end
#     end

#     @testset ">> b5" begin
#         b5 = Lattices.b5()
#         ops = space_group(b5)
#         @test length(ops) --> 48
#         @test count(ops) do op; all(abs(op.offset) .< 1e-12) end --> 12
#     end
# end

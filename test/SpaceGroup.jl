function valid_super_cell(size::Integer=3)
    cell = rand(-3:3, (size, size))
    while det(cell) ≤ 1
        cell = rand(-3:3, (3, 3))
    end
    cell
end

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
    allops = point_group([0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0])

    @testset "Cell: $cell" for (cell, numops) in parameters
        ops = point_group(cell)
        @test length(ops) == numops
        for op in ops
            @test size(op) == (3, 3)
            transformation = inv(cell) * op * cell
            @test all(isinteger, round.(transformation, 8))
            @test volume(transformation) ≈ 1.0

            if numops != 48
                found = 0
                for op in allops
                    transformation = inv(cell) * op * cell
                    all(isinteger, round.(transformation, 8)) || continue
                    isapprox(volume(transformation), 1, rtol=1e-8) || continue
                    found += 1
                end
                @test found == numops
            end
        end
    end
end

@testset "> Inner translations" begin
    @testset ">> Implementation details" begin
        a = DataFrame(species=["Al", "Al", "O", "O", "O"], momo=[1, 2, 3, 3, 1],
                      mama=[:+, :+, :+, :-, :-])

        species = Crystals.SpaceGroup.species_ids(a, [:species])
        @test species[1] == species[2]
        @test species[3] == species[4] == species[5]
        @test species[1] ≠ species[3]

        species = Crystals.SpaceGroup.species_ids(a, [:species, :momo])
        @test length(Set(species[[1, 2, 3, 5]])) == 4
        @test species[3] == species[4]
    end

    diamond = Lattices.diamond()
    diamond[:species] = ["A", "A"]
    @testset ">> No translations" begin
        translations = Crystals.SpaceGroup.inner_translations(diamond)
        @test size(translations, 2) == 0
        @test size(translations, 1) == size(diamond[:position], 1)
    end
    @testset ">> Supercell" begin
        cell = valid_super_cell(3)
        large = supercell(diamond, diamond.cell * cell)
        @test nrow(large) == round(Integer, det(cell)) * nrow(diamond)
        translations = Crystals.SpaceGroup.inner_translations(large)
        @test size(translations, 2) == round(Integer, det(cell) - 1)
        @test size(translations, 1) == size(diamond[:position], 1)
    end
end

@testset "> Make primitive" begin
    zinc_blende = Lattices.zinc_blende()
    @test is_primitive(zinc_blende)
    @test primitive(zinc_blende).cell ≈ zinc_blende.cell
    @test nrow(primitive(zinc_blende)) == nrow(zinc_blende)

    cell = valid_super_cell(3)
    large = supercell(zinc_blende, zinc_blende.cell * cell)
    @test !is_primitive(large)
    prim = primitive(large)
    @test niggly(prim.cell) ≈ niggly(zinc_blende.cell)
    @test nrow(prim) == nrow(zinc_blende)
    for i in 1:2
        prim_pos = prim[i, :position]
        d_pos = zinc_blende[i, :position]
        @test is_periodic(prim_pos, d_pos, prim.cell)
        @test is_periodic(prim_pos, d_pos, zinc_blende.cell)
    end
end

@testset "> Space group operations" begin
    @testset ">> fcc" begin
      fcc = Lattices.fcc()
      ops = space_group(fcc)
      @test length(ops) == 48
      for op in ops
          # No translations
          @test op([0, 0, 0]u"nm") ≈ [0, 0, 0]u"nm"
      end
    end

    @testset ">> b5" begin
        b5 = Lattices.b5()
        ops = space_group(b5)
        @test length(ops) == 48
        @test count(o -> all(abs.(ustrip.(o([0, 0, 0]u"nm"))) .< 1e-12), ops) == 12
    end
end

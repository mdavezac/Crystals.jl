facts("Potential equivalents regression") do
    cell = [1.5 -0.5 -0.5;  0  1 -0.5; 0.5 -0.5  1]
    avecs, bvecs, cvecs = Crystals.SpaceGroup.potential_equivalents(cell)
    norms = mapreducedim(x -> x*x, +, cell, 1)
    @fact size(avecs) --> (3, 4)
    @fact mapreducedim(x -> x * x, +, avecs, 1) .== norms[1] --> all

    @fact size(bvecs) --> (3, 8)
    @fact mapreducedim(x -> x * x, +, bvecs, 1) .== norms[2] --> all

    @fact size(cvecs) --> (3, 8)
    @fact mapreducedim(x -> x * x, +, cvecs, 1) .== norms[3] --> all
end

facts("Point group operations") do
    parameters = [
    ([0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0], 48),
    ([-0.5 0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 -0.5], 48),
    ([-0.6 0.5 0.5; 0.6 -0.5 0.5; 0.6 0.5 -0.5], 4),
    ([-0.7 0.7 0.7; 0.6 -0.5 0.5; 0.6 0.5 -0.5], 8),
    ([-0.765 0.7 0.7; 0.665 -0.5 0.5; 0.6 0.5 -0.5], 2)
    ]
    allops = point_group_operations([0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0])

    for (cell, numops) in parameters
        ops = point_group_operations(cell)
        @fact length(ops) --> numops
        for op in ops
            @fact size(op.scalefwd) --> (3, 3)
            @fact op.offset --> roughly([0, 0, 0])
            transformation = inv(cell) * op.scalefwd * cell
            @fact transformation --> all_integers
            @fact abs(det(transformation)) --> roughly(1)

            if numops != 48
                failed = 0
                for op in allops
                    transformation = inv(cell) * op.scalefwd * cell
                    all_integers(transformation) || continue
                    abs(abs(det(transformation)) - 1) < 1e-8 || continue
                    failed += 1
                end
                @fact failed --> numops
            end
        end
    end
end

facts("Inner translations") do
    diamond = Lattices.diamond()
    context("No translations") do
        translations = inner_translations(diamond)
        @fact length(translations) --> 0
    end
    context("Supercell") do
        cell = -1
        while det(cell) ≤ 0
          cell = rand(-3:3, (3, 3))
        end
        large = supercell(diamond, diamond.cell * cell)
        @fact nrow(large) --> round(Integer, det(cell)) * nrow(diamond)
        translations = inner_translations(large)
        @fact length(translations) --> round(Integer, det(cell) - 1)
    end
end

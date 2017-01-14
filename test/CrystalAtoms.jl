crystal = Crystal(eye(2)u"nm", species=["Al", "O", "O"],
                  position=[1 1 1; 2 3 4], label=[:+, :-, :-])

@testset "> CrystalAtom, row $i" for i in 1:size(crystal, 1)
    atom = Crystals.CrystalAtoms.CrystalAtom(crystal, i)
    @test atom[:position] == crystal[i, :position]
    @test atom[:species] == crystal[i, :species]
    @test atom[:label] == crystal[i, :label]
end

@testset "> Iteration, row $i" for (i, atom) in enumerate(eachatom(crystal))
    @test atom[:position] == crystal[i, :position]
    @test atom[:species] == crystal[i, :species]
    @test atom[:label] == crystal[i, :label]
end

@testset "> map" begin
    iterated = Any[]
    map(x->push!(iterated, iterated), eachatom(crystal))
    @test length(iterated) == length(crystal)
end

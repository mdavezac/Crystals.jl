crystal = Crystal(eye(2)u"nm", species=["Al", "O", "O"],
                  position=[1 1 1; 2 3 4], label=[:+, :-, :-])

@testset "> CrystalAtom, row $i" for i in 1:size(crystal, 1)
    atom = Crystals.CrystalAtoms.CrystalAtom(crystal, i)
    @test atom[:position] == crystal[i, :position]
    @test atom[:position, 1] == crystal[i, :position, 1]
    @test atom[:species] == crystal[i, :species]
    @test atom[:label] == crystal[i, :label]
    atom[:label] = Symbol("a$i")
    @test atom[:label] == crystal[i, :label] == Symbol("a$i")
    atom[:position, 1] = 2i - 15
    @test atom[:position, 1] == crystal[i, :position, 1] == 2i - 15
    @test names(atom) == names(crystal)
end

crystal = Crystal(eye(2)u"nm", species=["Al", "O", "O"],
                  position=[1 1 1; 2 3 4], label=[:+, :-, :-])
@test eltype(eachrow(crystal)) === Crystals.CrystalAtoms.CrystalAtom{typeof(crystal)}
@testset "> Iteration, row $i" for (i, atom) in enumerate(eachrow(crystal))
    @test atom[:position] == crystal[i, :position]
    @test atom[:species] == crystal[i, :species]
    @test atom[:label] == crystal[i, :label] != Symbol("a$i")
    atom[:label] = Symbol("a$i")
    @test atom[:label] == crystal[i, :label] == Symbol("a$i")
    atom[:position, 1] = 2i - 15
    @test atom[:position, 1] == crystal[i, :position, 1] == 2i - 15
end

crystal = Crystal(eye(2)u"nm", species=["Al", "O", "O"],
                  position=[1 1 1; 2 3 4], label=[:+, :-, :-])
@testset "> map" begin
    iterated = Any[]
    map(x->push!(iterated, iterated), eachatom(crystal))
    @test length(iterated) == length(crystal)
    @test size(iterated) == (length(crystal), )
end

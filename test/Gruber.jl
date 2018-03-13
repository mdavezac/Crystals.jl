@testset "> Gruber" begin
    @testset ">>Is integer matrices (no units)" begin
        cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]
        const lim = 5

        for a00=-1:1, a10 =-lim:lim + 1, a11=-1:1, a20=-lim:lim + 1, a21=-lim:lim + 1, a22=-1:1
            a = [a00 0 0; a10 a11 0; a20 a21 a22]
            abs(det(a)) > 0 || continue
            g = gruber(cell * a)
            @test abs(det(cell)) > 1e-12
            @test volume(cell * a) ≈ volume(g)
            @test all(isinteger, inv(cell) * g)
            @test all(isinteger, round.(inv(g) * cell, 8))
        end
    end

    @testset ">>Is integer matrices (with units)" begin
        cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]u"nm"
        const lim = 5

        for a00=-1:1, a10 =-lim:lim + 1, a11=-1:1, a20=-lim:lim + 1, a21=-lim:lim + 1, a22=-1:1
            a = [a00 0 0; a10 a11 0; a20 a21 a22]
            abs(det(a)) > 0 || continue
            g = gruber(cell * a)
            @test abs(det(ustrip.(cell))) > 1e-12
            @test volume(cell * a) ≈ volume(g)
            @test all(isinteger, inv(cell) * g)
            @test all(isinteger, round.(inv(g) * cell, 8))
            break
        end
    end

    @testset ">> Regressions" begin
          cells = Any[
            [0.10658442  0.85196214  0.11820731;
              0.08762888  0.06597675  0.33166737;
              0.90830071  0.79275198  0.28060048],
            [0.315935703201265183 -0.491002658691486404 -0.450428983281787709;
              0.143383123347384012 -0.208941528033948387  0.315007354520186578;
              0.4831297491880433   -0.127532144165201089  0.433845061278055],
            -[0.5  -0.5   1.0; 0.0  -0.5  -0.5; -0.5   0.0  -0.5]
          ]
          expecteds = Any[
              [0.118207309999999996 0.01162288999999999684 -0.74537772000000002137;
              0.33166737000000001689 0.2440384899999999968 0.02165213000000000554;
              0.28060047999999998547 -0.62770023000000008118 0.11554872999999998839],
              [0.175066955490221221  0.275362027791566488 0.491002658691486404;
              0.065558404686564375 -0.380565759206750953 0.208941528033948387;
              -0.355597605022842211 -0.078247456255212788 0.127532144165201089],
              [-0.5 -0.5  0;  0 -0.5 -0.5; 0.5 0  0.5]
          ]
          for (cell, expected) in zip(cells, expecteds)
            actual = gruber(cell)
            @test actual ≈ expected
          end
    end

    @testset ">> Error if singular" begin
        @test_throws ErrorException gruber([1 0 0; 0 1 0; 0 2 0])
    end

    @testset ">> Error if max iter" begin
        @test_throws ErrorException gruber([1 0 0; 1 1  0; 4 5 1]; itermax=2)
    end
end

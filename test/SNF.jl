@testset "> Smith Normal Form" begin
  @testset ">> choose pivot" begin

    @testset ">>> With zero columns" begin
      smith = [0 4 4; 0 6 12; 0 -4 -16]
      left = eye(eltype(smith), 3)
      jt = Crystals.SNF.choose_pivot!(left, smith, 1, 1)
      @test jt == 2
      @test smith == [0 4 4; 0 6 12; 0 -4 -16]
      @test left == eye(eltype(smith), 3)
    end

    @testset ">>> Switch rows" begin
      smith = [0 0 4; 0 6 12; 0 -4 -16]
      original = deepcopy(smith)
      left = eye(eltype(smith), 3)
      jt = Crystals.SNF.choose_pivot!(left, smith, 1, 1)
      @test jt == 2
      @test smith == [0 6 12; 0 0 4; 0 -4 -16]
      @test left == [0 1 0; 1 0 0; 0 0 1]
      @test left * original == smith
    end
  end

  @testset ">> Improve col pivot" begin
    @testset ">>> All multiples" begin
      smith = [1 0 4; 0 2 12; 0 -4 -16]
      left = eye(eltype(smith), 3)
      Crystals.SNF.improve_col_pivot!(left, smith, 1, 1)
      @test left == eye(eltype(smith), 3)
      @test smith == [1 0 4; 0 2 12; 0 -4 -16]   

      Crystals.SNF.improve_col_pivot!(left, smith, 2, 2)
      @test left == eye(eltype(smith), 3)
      @test smith == [1 0 4; 0 2 12; 0 -4 -16]   
    end

    @testset ">>> Create multiple" begin
      smith = [1 0 4; 0 3 12; 0 -4 -16]
      original = deepcopy(smith)
      left = eye(eltype(smith), 3)
      Crystals.SNF.improve_col_pivot!(left, smith, 2, 2)
      @test all(smith[:, 2] .% smith[2, 2] .== 0)
      @test left * original == smith
      @test abs(det(left)) ≥ 0
    end
  end

  @testset ">> Improve row pivot" begin
    @testset ">>> All multiples" begin
      smith = transpose([1 0 4; 0 2 12; 0 -4 -16])
      right = eye(eltype(smith), 3)
      Crystals.SNF.improve_row_pivot!(smith, right, 1, 1)
      @test right == eye(eltype(smith), 3)
      @test smith == transpose([1 0 4; 0 2 12; 0 -4 -16])

      Crystals.SNF.improve_row_pivot!(smith, right, 2, 2)
      @test right == eye(eltype(smith), 3)
      @test smith == transpose([1 0 4; 0 2 12; 0 -4 -16]   )
    end

    @testset ">>> Create multiple" begin
      smith = transpose([1 0 4; 0 3 12; 0 -4 -16])
      original = deepcopy(smith)
      right = eye(eltype(smith), 3)
      Crystals.SNF.improve_row_pivot!(smith, right, 2, 2)
      @test all(smith[2, :] .% smith[2, 2] .== 0)
      @test original * right == smith
      @test abs(det(right)) ≥ 0
    end
  end

  @testset ">> Diagnonalize one row-column" begin
    smith = [1 2 4; 3 3 12; 5 -4 -16]
    original = deepcopy(smith)
    right = eye(eltype(smith), 3)
    left = eye(eltype(smith), 3)
    Crystals.SNF.diagonalize_at_entry!(left, smith, right, 2, 2)
    @test smith[2, 2] ≠ 0
    @test countnz(smith[:, 2]) == 1
    @test countnz(smith[2, :]) == 1
    @test left * original * right == smith
    @test abs(det(left)) ≥ 0
    @test abs(det(right)) ≥ 0
  end

  @testset ">> Diaginalize all entries" begin
    original = [1 2 4; 3 3 12; 5 -4 -16]
    smith, left, right = Crystals.SNF.diagonalize_all_entries(original)
    @test isdiag(smith)
    @test left * original * right == smith
    @test abs(det(left)) ≥ 0
    @test abs(det(right)) ≥ 0
  end

  @testset ">> Move zero columns left" begin
    matrices = Any[
      [0 1 0; 0 2 2; 0 3 3],
      [0 1 0; 0 0 1; 0 0 0]
    ]
    for smith in matrices
      original = deepcopy(smith)
      right = eye(eltype(smith), 3)
      Crystals.SNF.move_zero_entries!(smith, right)
      @test smith[:, 1] == original[:, 2]
      @test smith[:, 2] == original[:, 3]
      @test smith[:, 3] == original[:, 1]
      @test original * right == smith
      @test abs(det(right)) ≥ 0
    end
  end

  @testset ">> The whole shebang" begin
    matrices = vcat(
      Any[
        BigInt[1 2; -4 5],
        BigInt[0 -4 -5; 0 0 -4; 0 3 -1],
        BigInt[2 4 4; -6 6 12; 10 -4 -16],
        BigInt[-5 -5 2 -2; -4 -4 0 -1; -1 -1 4 -1; 0 -4 4 -4]
      ], 
      Any[rand(-BigInt(5):5, tuple(repeat(rand(2:5, (1, )), inner=[2])...)) for u in 1:10]
    )
    @testset ">>> try $matrix" for matrix in matrices
      smith, left, right = Crystals.SNF.smith_normal_form(matrix)
      @test isdiag(smith)
      @test left * matrix * right == smith
      k = findfirst(x -> x == 0, diag(smith))
      if k ≠ 0
      @test countnz(diag(smith)[k:end]) == 0
      else
        k = size(smith, 2) + 1
      end
      for i in 2:(k - 1)
        @test smith[i, i] % smith[i - 1, i - 1] == 0
      end
      @test all(diag(smith) .≥ 0)
    end
  end

  @testset ">> Units" begin
      matrix = BigInt[1 2; -4 5]
      usmith, uleft, uright = Crystals.SNF.smith_normal_form(matrix * u"nm")
      smith, left, right = Crystals.SNF.smith_normal_form(matrix)
      @test  uleft == left
      @test  uright == right
      @test eltype(usmith) == eltype(matrix * u"nm")
      usmith = reshape([ustrip.(u) for u in usmith], size(usmith))
      @test usmith == smith
  end
end

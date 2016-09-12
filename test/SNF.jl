facts("Smith Normal Form") do
  context("choose pivot") do

    context("With zero columns") do
      smith = [0 4 4; 0 6 12; 0 -4 -16]
      left = eye(eltype(smith), 3)
      jt = Crystals.SNF.choose_pivot!(left, smith, 1, 1)
      @fact jt --> 2
      @fact smith --> [0 4 4; 0 6 12; 0 -4 -16]
      @fact left --> eye(eltype(smith), 3)
    end

    context("Switch rows") do
      smith = [0 0 4; 0 6 12; 0 -4 -16]
      original = deepcopy(smith)
      left = eye(eltype(smith), 3)
      jt = Crystals.SNF.choose_pivot!(left, smith, 1, 1)
      @fact jt --> 2
      @fact smith --> [0 6 12; 0 0 4; 0 -4 -16]
      @fact left --> [0 1 0; 1 0 0; 0 0 1]
      @fact left * original --> smith
    end
  end

  context("Improve col pivot") do
    context("All multiples") do
      smith = [1 0 4; 0 2 12; 0 -4 -16]
      left = eye(eltype(smith), 3)
      Crystals.SNF.improve_col_pivot!(left, smith, 1, 1)
      @fact left --> eye(eltype(smith), 3)
      @fact smith --> [1 0 4; 0 2 12; 0 -4 -16]   

      Crystals.SNF.improve_col_pivot!(left, smith, 2, 2)
      @fact left --> eye(eltype(smith), 3)
      @fact smith --> [1 0 4; 0 2 12; 0 -4 -16]   
    end

    context("Create multiple") do
      smith = [1 0 4; 0 3 12; 0 -4 -16]
      original = deepcopy(smith)
      left = eye(eltype(smith), 3)
      Crystals.SNF.improve_col_pivot!(left, smith, 2, 2)
      @fact smith[:, 2] .% smith[2, 2] .== 0 --> all
      @fact left * original --> smith
      @fact abs(det(left)) --> greater_than(0)
    end
  end

  context("Improve row pivot") do
    context("All multiples") do
      smith = transpose([1 0 4; 0 2 12; 0 -4 -16])
      right = eye(eltype(smith), 3)
      Crystals.SNF.improve_row_pivot!(smith, right, 1, 1)
      @fact right --> eye(eltype(smith), 3)
      @fact smith --> transpose([1 0 4; 0 2 12; 0 -4 -16])

      Crystals.SNF.improve_row_pivot!(smith, right, 2, 2)
      @fact right --> eye(eltype(smith), 3)
      @fact smith --> transpose([1 0 4; 0 2 12; 0 -4 -16]   )
    end

    context("Create multiple") do
      smith = transpose([1 0 4; 0 3 12; 0 -4 -16])
      original = deepcopy(smith)
      right = eye(eltype(smith), 3)
      Crystals.SNF.improve_row_pivot!(smith, right, 2, 2)
      @fact smith[2, :] .% smith[2, 2] .== 0 --> all
      @fact original * right --> smith
      @fact abs(det(right)) --> greater_than(0)
    end
  end

  context("Diagnonalize one row-column") do
    smith = [1 2 4; 3 3 12; 5 -4 -16]
    original = deepcopy(smith)
    right = eye(eltype(smith), 3)
    left = eye(eltype(smith), 3)
    Crystals.SNF.diagonalize_at_entry!(left, smith, right, 2, 2)
    @fact smith[2, 2] --> not(0)
    @fact countnz(smith[:, 2]) --> 1
    @fact countnz(smith[2, :]) --> 1
    @fact left * original * right --> smith
    @fact abs(det(left)) --> greater_than(0)
    @fact abs(det(right)) --> greater_than(0)
  end

  context("Diaginalize all entries") do
    original = [1 2 4; 3 3 12; 5 -4 -16]
    smith, left, right = Crystals.SNF.diagonalize_all_entries(original)
    @fact smith --> isdiag
    @fact left * original * right --> smith
    @fact abs(det(left)) --> greater_than(0)
    @fact abs(det(right)) --> greater_than(0)
  end

  context("Move zero columns left") do
    matrices = Any[
      [0 1 0; 0 2 2; 0 3 3],
      [0 1 0; 0 0 1; 0 0 0]
    ]
    for smith in matrices
      original = deepcopy(smith)
      right = eye(eltype(smith), 3)
      Crystals.SNF.move_zero_entries!(smith, right)
      @fact smith[:, 1] --> original[:, 2]
      @fact smith[:, 2] --> original[:, 3]
      @fact smith[:, 3] --> original[:, 1]
      @fact original * right --> smith
      @fact abs(det(right)) --> greater_than(0)
    end
  end

  context("The whole shebang") do
    matrices = vcat(
      Any[
        BigInt[1 2; -4 5],
        BigInt[0 -4 -5; 0 0 -4; 0 3 -1],
        BigInt[2 4 4; -6 6 12; 10 -4 -16],
        BigInt[-5 -5 2 -2; -4 -4 0 -1; -1 -1 4 -1; 0 -4 4 -4]
      ], 
      Any[rand(-BigInt(5):5, tuple(repeat(rand(2:5, (1, )), inner=[2])...)) for u in 1:10]
    )
    for matrix in matrices
      context("try $matrix") do 
        smith, left, right = Crystals.SNF.smith_normal_form(matrix)
        @fact smith --> isdiag
        @fact left * matrix * right --> smith
        k = findfirst(x -> x == 0, diag(smith))
        if k ≠ 0
        @fact countnz(diag(smith)[k:end]) --> 0
        else
          k = size(smith, 2) + 1
        end
        for i in 2:(k - 1)
          @fact smith[i, i] % smith[i - 1, i - 1] --> 0
        end
        @fact diag(smith) .≥ 0 --> all
      end
    end
  end
end

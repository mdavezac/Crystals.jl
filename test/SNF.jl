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
      left = eye(eltype(smith), 3)
      jt = Crystals.SNF.choose_pivot!(left, smith, 1, 1)
      @fact jt --> 2
      @fact smith --> [0 6 12; 0 0 4; 0 -4 -16]
      @fact left --> [0 1 0; 1 0 0; 0 0 1]
      @fact left * smith --> [0 0 4; 0 6 12; 0 -4 -16]
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
      left = eye(eltype(smith), 3)
      Crystals.SNF.improve_col_pivot!(left, smith, 2, 2)
      @fact smith[:, 2] .% smith[2, 2] .== 0 --> all
      @fact left * smith --> [1 0 4; 0 3 12; 0 -4 -16]   
    end
  end

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
    right = eye(eltype(smith), 3)
    Crystals.SNF.improve_row_pivot!(smith, right, 2, 2)
    @fact smith[2, :] .% smith[2, 2] .== 0 --> all
    @fact smith * right --> transpose([1 0 4; 0 3 12; 0 -4 -16])
  end
end

facts("Smith Normal Form") do
  context("Implementation functions") do
    get_min_max = Crystals.SNF.get_min_max
    for i in 1:10
      sequence = shuffle(collect(0:10))
      min = findfirst(x -> x == 1, sequence)
      max = findfirst(x -> x == 10, sequence)
      @fact get_min_max(sequence) --> (min, max)
    end
  end
  facts("Random Matrices") do
    for test in 1:1
      matrix = rand(-4:4, (3, 3))
      left, S, right = smith_normal_form(matrix)
      @fact diagm(diag(S)) --> S
      @fact left * S * transpose(right) --> matrix
    end
  end
end

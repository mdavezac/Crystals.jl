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

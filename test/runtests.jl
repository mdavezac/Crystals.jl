module CrystalTest
using Crystal
using FactCheck

facts("Construction") do
  context("Empty") do
    structure = Structure()
    @fact all(abs(structure.cell - eye(3)) .< 1e-8)
  end
end

FactCheck.exitstatus()
end

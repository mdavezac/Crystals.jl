function is_integer(x::Matrix; tolerance=1e-12)
  for u in x
    abs(u - round(u)) < tolerance || return false
  end
  true
end
function is_integer(; tolerance=1e-12)
  x -> is_integer(x; tolerance=tolerance)
end

facts("Gruber cell") do
  context("Is integer matrices") do
    cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]
    const lim = 5
 
    for a00=-1:1, a10 =-lim:lim + 1, a11=-1:1, a20=-lim:lim + 1, a21=-lim:lim + 1, a22=-1:1
      a = [a00 0 0; a10 a11 0; a20 a21 a22]
      abs(det(a)) >0 || continue
      g = gruber(cell * a)
      @fact abs(det(cell)) --> greater_than(1e-12)
      @fact inv(cell) * g --> is_integer
      @fact inv(g) * cell --> is_integer
    end
  end

  context("Regressions") do
    cells = Any[
       [0.10658442  0.85196214  0.11820731;
        0.08762888  0.06597675  0.33166737;
        0.90830071  0.79275198  0.28060048],
       [0.315935703201265183 -0.491002658691486404 -0.450428983281787709;
        0.143383123347384012 -0.208941528033948387  0.315007354520186578;
        0.4831297491880433   -0.127532144165201089  0.433845061278055]
    ]
    expecteds = Any[
        [0.118207309999999996 0.01162288999999999684 -0.74537772000000002137;
         0.33166737000000001689 0.2440384899999999968 0.02165213000000000554;
         0.28060047999999998547 -0.62770023000000008118 0.11554872999999998839],
        [0.175066955490221221  0.275362027791566488 0.491002658691486404;
         0.065558404686564375 -0.380565759206750953 0.208941528033948387;
        -0.355597605022842211 -0.078247456255212788 0.127532144165201089]
    ]
    for (cell, expected) in zip(cells, expecteds)
      actual = gruber(cell)
      @fact actual --> roughly(expected)
    end
  end

  context("Error if singular") do
    @fact_throws ErrorException gruber([1 0 0; 0 1 0; 0 2 0])
  end

  context("Error if max iter") do
    @fact_throws ErrorException gruber([1 0 0; 1 1  0; 4 5 1]; itermax=2)
  end
end

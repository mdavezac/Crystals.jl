function is_integer(x::Matrix; tolerance=1e-12)
  for u in x
    abs(u - round(u)) < tolerance || return false
  end
  true
end
function is_integer(; tolerance=1e-12)
  x -> is_integer(x; tolerance=tolerance)
end

facts("Gruber") do
  cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]
  lim = 5

  for a00=-1:1, a10 =-lim:lim + 1, a11=-1:1, a20=-lim:lim + 1, a21=-lim:lim + 1, a22=-1:1
    a = [a00 0 0; a10 a11 0; a20 a21 a22]
    abs(det(a)) >0 || continue
    g = gruber(cell * a)
    @fact abs(det(cell)) --> greater_than(1e-12)
    @fact inv(cell) * g --> is_integer
    @fact inv(g) * cell --> is_integer
  end
end

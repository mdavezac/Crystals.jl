"""
Hart-Forcade transform

Computes the cyclic group of a supercell with respect to a lattice. It makes it
possible to identify the class of periodically equivalent cell that a given
position within the supercell belongs to.

Returns the transform and the quotient.
"""
function hart_forcade(lattice::Matrix, supercell::Matrix; digits=8)
  fractional = convert(Matrix{Int64}, round(inv(lattice) * supercell, digits))

  snf, left, right = smith_normal_form(fractional)

  left * inv(cell), diag(snf)
end

"""
True if the two positions are periodic
  
Returns a boolean if the input are two positions, and an array of booleans if
the input are arrays of positions.
"""
is_periodic(a::Matrix, b::Matrix, cell::Matrix; tolerance=default_tolerance) =
  all(abs(origin_centered(a - b, cell)) .< tolerance, 1)

is_periodic(a::Vector, b::Vector, cell::Matrix; tolerance=default_tolerance) =
  all(abs(origin_centered(a - b, cell)) .< tolerance)

""" Folds periodic positions into cell """
into_cell(pos::Array, cell::Matrix) = cell * mod(inv(cell) * pos, 1)

""" Folds vector back to origin """
origin_centered(pos::Array, cell::Matrix) =
    cell * (mod(inv(cell) * pos .+ 0.5, -1) .+ 0.5)

"""
Folds vector into first Brillouin zone of the input cell

Returns the periodic image with the smallest possible norm.
"""
function into_voronoi(pos, cell)
  result = origin_centered(pos, cell)
  zcentered = deepcopy(result)
  norms = [norm(result[:, i]) for i in 1:size(result, 2)]
  for n in 1:length(norms)
    for i = -1:1:2, j = -1:1:2, k = -1:1:2
      translation = cell * [i, j, k]
      position = zcentered[:, n] + translation 
      d = norm(position)
      if d < norms[n]
        result[:, n] = position
        norms[n] = d
      end
    end
  end
  result
end

# """ Creates a supercell of an input lattice """
# function supercell(lattice::Crystal, supercell):
#   length(lattice) == 0 && error("Lattice is empty")
#   result = Crystal(lattice.cell, lattice.scale)
#
#   transform, quotient = hart_forcade(lattice.cell, result.cell)
#   itransform = inv(transform)
#   icell = inv(cell)
#   for i = 1:quotient[1], j = 1:quotient[2], k = 1:quotient[3]
#     atoms = deepcopy(lattice.atoms)
#     atoms[:position] = 
#   end
#
#     invtransform = inv(transform.transform)
#     invcell = inv(result.cell)
#
#     for i in range(transform.quotient[0]):
#       for j in range(transform.quotient[1]):
#         for k in range(transform.quotient[2]):
#           pos = dot(invtransform, array([i, j, k], dtype='float64'))
#
#           for l, site in enumerate(lattice):
#             atom = site.copy()
#             atom.pos = into_cell(pos + site.pos, result.cell, invcell)
#             atom.site = l
#             result.append(atom)
#             return result;

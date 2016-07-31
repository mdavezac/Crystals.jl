module SpaceGroup
  """
  gvectors with equivalent norms

  Figures out vectors in the lattice defined by the cell which have the same
  lengths as the column vectors defining the cell. these new vectors are
  potentially symmetrically equivalent.
  """
  function potential_equivalents(cell::Matrix; tolerance::Real=1e-8)
    const volume = abs(det(cell))
    const a0 = cell[:, 1]
    const a1 = cell[:, 2]
    const a2 = cell[:, 3]

    lengths = reducedim(+, cell .* cell, 1)
    max_norm = mapreduce(i -> norm(cell[:, i]), max, 0, 1:size(cell, 2))

    const n0 = Int(ceil(max_norm * norm(cross(a1, a2)) / volume))
    const n1 = Int(ceil(max_norm * norm(cross(a2, a0)) / volume))
    const n2 = Int(ceil(max_norm * norm(cross(a0, a1)) / volume))

    gvectors = Any[Array{eltype(cell), 1}[] for u in 1:length(lengths)]
    for i in -n0:n0, j in -n1:n1, k in -n2:n2
      g = cell * eltype(cell)[i, j, k]
      glength = sum(g .* g)
      for (length, result) in zip(lengths, gvectors)
        if abs(length - glength) < tolerance
          push!(result, g)
        end
      end
    end

    [hcat(gs...) for gs in gvectors]
  end
end

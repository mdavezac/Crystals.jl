module SpaceGroup
  " Computes all gvectors in prolate defined by cell "
  function gvectors(cell::Matrix; tolerance::Real=1e-8):
    const volume = abs(det(cell))
    const a0 = cell[:, 0]
    const a1 = cell[:, 1]
    const a2 = cell[:, 2]

    lengths = reducedim(+, cell .* cell, 2)
    max_norm = mapreduce(i -> norm(cell[:, i]), max, 0, 1:size(cell, 2))

    const n0 = ceil(max_norm * norm(cross(a1, a2)) / volume)
    const n1 = ceil(max_norm * norm(cross(a2, a0)) / volume)
    const n2 = ceil(max_norm * norm(cross(a0, a1)) / volume)

    gvectors = repeat(Any[], length(lengths))
    for i in -n0:n0:
        for j in -n1:n1:
            for k in -n2:n2:
                g = cell * [i, j, k]
                glength = sum(g .* g)
                for length, result in zip(lengths, gvectors):
                    if abs(length - glength) < tolerance:
                        push!(result, g)
    return gvectors



end

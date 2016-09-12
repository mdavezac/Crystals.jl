module SNF
export smith_normal_form

function choose_pivot!{T <: Integer}(left::Matrix{T}, smith::Matrix{T},
                                     step::Integer, index::Integer=1)
  index > size(smith, 2) && return 0
  while all(smith[:, index] .== 0)
    index += 1
    index ≤ size(smith, 2) || return 0;
  end
  if smith[step, index] == 0
    k = findfirst(x -> x ≠ 0, smith[:, index])
    left[step, :], left[k, :] = deepcopy(left[k, :]), deepcopy(left[step, :])
    smith[step, :], smith[k, :] = deepcopy(smith[k, :]), deepcopy(smith[step, :])
  end
  index
end

function improve_col_pivot!{T <: Integer}(left::Matrix{T}, smith::Matrix{T},
                                          row::Integer, col::Integer)
  @assert size(left, 1) == size(left, 2) == size(smith, 1) == size(smith, 2)
  for k in 1:size(smith, 1)
    smith[k, col] % smith[row, col] == 0 && continue

    β, σ, τ = gcdx(smith[row, col], smith[k, col])
    α = smith[row, col] / β
    γ = smith[k, col] / β

    Lp = eye(T, size(left, 1))
    Lp[row, [row, k]] = [σ, τ]
    Lp[k, [row, k]] = [-γ, α]

    left[:] = Lp * left
    smith[:] = Lp * smith
  end
end

function improve_row_pivot!{T <: Integer}(smith::Matrix{T}, right::Matrix{T},
                                          row::Integer, col::Integer)
  @assert size(right, 1) == size(right, 2) == size(smith, 1) == size(smith, 2)
  for k in 1:size(smith, 1)
    smith[row, k] % smith[row, col] == 0 && continue

    β, σ, τ = gcdx(smith[row, col], smith[row, k])
    α = smith[row, col] / β
    γ = smith[row, k] / β

    Rp = eye(T, size(right, 1))
    Rp[[col, k], col] = [σ, τ]
    Rp[[col, k], k] = [-γ, α]

    right[:] = right * Rp
    smith[:] = smith * Rp
  end
end

function diagonalize_at_entry!{T <: Integer}(left::Matrix{T}, smith::Matrix{T},
                                             right::Matrix{T}, row::Integer, col::Integer)
  @assert size(right, 1) == size(right, 2) == size(smith, 1) == size(smith, 2)
  while countnz(smith[:, col]) > 1 || countnz(smith[row, :]) > 1
    improve_col_pivot!(left, smith, row, col)
    for i in 1:size(left, 2)
      i == row && continue
      smith[i, col] == 0 && continue
      const β = smith[i, col] ÷ smith[row, col]
      left[i, :] -= left[row, :] * β
      smith[i, :] -= smith[row, :] * β
    end

    improve_row_pivot!(smith, right, row, col)
    for i in 1:size(left, 1)
      i == col && continue
      smith[row, i] == 0 && continue
      const β = smith[row, i] ÷ smith[row, col]
      right[:, i] -= right[:, col] * β
      smith[:, i] -= smith[:, col] * β
    end
  end
end

function diagonalize_all_entries!{T <: Integer}(left::Matrix{T}, smith::Matrix{T}, right::Matrix{T})
  @assert size(smith, 1) == size(smith, 2)
  step = 1
  col = choose_pivot!(left, smith, step)
  while 0 < col ≤ size(smith, 2)
    diagonalize_at_entry!(left, smith, right, step, col)

    step += 1
    col += 1
    col = choose_pivot!(left, smith, step, col)
  end
  return smith, left, right
end

function diagonalize_all_entries{T <: Integer}(smith::Matrix{T})
  smith = deepcopy(smith)
  left = eye(eltype(smith), size(smith, 1))
  right = eye(eltype(smith), size(smith, 2))
  diagonalize_all_entries!(left, smith, right)
  return smith, left, right
end

function move_zero_entries!{T <: Integer}(smith::Matrix{T}, right::Matrix{T})
  nonzero = find(1:size(smith, 2)) do i; any(smith[:, i] .≠ 0) end
  length(nonzero) == size(smith, 2) && return
  for (i, j) in enumerate(nonzero)
    i == j && continue
    smith[:, i], smith[:, j] = deepcopy(smith[:, j]), deepcopy(smith[:, i])
    right[:, i], right[:, j] = deepcopy(right[:, j]), deepcopy(right[:, i])
  end
end

function make_divisible!{T <: Integer}(left::Matrix{T}, smith::Matrix{T}, right::Matrix{T})
  divides = x -> (smith[x, x] ≠ 0 && smith[x, x] % smith[x - 1, x - 1] ≠ 0)
  while (i = findfirst(divides, 2:size(smith, 2))) ≠ 0
    (smith[i + 1, i + 1] == 0 || smith[i + 1, i + 1] % smith[i, i] == 0 ) && continue
    smith[:, i] += smith[:, i + 1]
    right[:, i] += right[:, i + 1]
    diagonalize_all_entries!(left, smith, right)
  end
end

function smith_normal_form{T <: Integer}(matrix::Matrix{T})
  smith, left, right = diagonalize_all_entries(matrix)
  move_zero_entries!(smith, right)
  make_divisible!(left, smith, right)
  left = convert(Matrix{T}, diagm([u ≥ 0 ? 1: -1 for u in diag(smith)])) * left
  smith = abs(smith)
  smith, left, right
end

end

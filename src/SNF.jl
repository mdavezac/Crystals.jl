module SNF
export smith_normal_form

function choose_pivot!{T <: Integer}(left::Matrix{T}, smith::Matrix{T},
                                     step::Integer, index::Integer=1)
  while all(smith[:, index] .== 0)
    index += 1
    index ≤ size(smith, 2) || return 0;
  end
  if smith[step, index] == 0
    k = findfirst(x -> x ≠ 0, smith[:, index])
    left[step, :], left[k, :] = left[k, :], left[step, :]
    smith[step, :], smith[k, :] = smith[k, :], smith[step, :]
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

    Lpinv = eye(T, size(left, 1))
    Lpinv[row, [row, k]] = [α, -τ]
    Lpinv[k, [row, k]] = [γ, σ]

    left[:] = left * Lpinv
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

    Rpinv = eye(T, size(right, 1))
    Rpinv[[col, k], col] = [α, -τ]
    Rpinv[[col, k], k] = [γ, σ]

    right[:] = Rpinv * right
    smith[:] = smith * Rp
  end
end

end

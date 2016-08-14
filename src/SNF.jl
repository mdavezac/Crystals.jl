module SNF
export smith_normal_form

const default_tolerance = 1e-12

function one_nonzero(vector::Array)
  found = false
  for v in vector
    if v != 0
      found && return false
      found = true
    end
  end
  true
end

function check_nonzero(matrix::Matrix, index::Int)
  found = false
  for i in 1:size(matrix, 1)
    if matrix[i, index] != 0
      found && return false
      found = true
    end
  end
  for i in 1:size(matrix, 2)
    if i != index && matrix[index, i] != 0
      found && return false
      found = true
    end
  end
  true
end

function get_min_max(sequence)
  min, max = nothing, 1
  for k in 2:length(sequence)
    abs(sequence[k]) > abs(sequence[max]) && (max = k)
    sequence[k] != 0 &&
      (min === nothing || abs(sequence[k]) < abs(sequence[min])) &&
      (min = k)
  end
  @assert min === nothing || sequence[min] != 0
  min, max
end

function check_same_magnitude(matrix::Matrix, index::Int)
  minelem, maxelem = get_min_max(matrix[:, index])
  @assert minelem !== nothing
  abs(matrix[minelem, index]) != abs(matrix[maxelem, index]) &&
      return minelem, maxelem
  n0, n1 = countnz(matrix[maxelem, :]), countnz(matrix[minelem, :])
  if n0 < n1 || (n0 == n1 && matrix[maxelem, index] < matrix[minelem, index])
    return maxelem, minelem
  end
  minelem, maxelem
end

function smith_col_impl(left::Matrix, smith::Matrix, index::Int)
  L, S = copy(left), copy(smith)
  while !one_nonzero(S[:, index])
    minelem, maxelem = check_same_magnitude(S, index)
    const multiple = S[maxelem, index] ÷ S[minelem, index]
    S[maxelem, :] -= multiple * S[minelem, :]
    L[maxelem, :] -= multiple * L[minelem, :]
  end
  if S[index, index] == 0
    k = findfirst(x -> x != 0, S[:, index])
    k == 0 && error("Internal error ", S[:, index])
    S[k, :], S[index, :] = copy(S[index, :]), copy(S[k, :])
    L[k, :], L[index, :] = copy(L[index, :]), copy(L[k, :])
  end
  if S[index, index] < 0
    S[index, :] *= -1
    L[index, :] *= -1
  end
  L, S
end

function smith_row_impl(smith::Matrix, right::Matrix, index::Int)
  R, S = smith_col_impl(transpose(right), transpose(smith), index)
  transpose(S), transpose(R)
end

function smith_normal_form{T <: Integer}(matrix::Matrix{T})
  det(matrix) == 0 && error("Input matrix is singular")
  size(matrix, 1) != size(matrix, 2) && error("Input matrix is not square")
  smith = copy(matrix)
  left = eye(eltype(smith), 3)
  right = eye(eltype(smith), 3)
  for index in 1:size(smith, 1)
    once = true
    while once || !check_nonzero(smith, index)
      once = false
      left, smith = smith_col_impl(left, smith, index)
      smith, right = smith_row_impl(smith, right, index)

      const diag = smith[index, index]
      maxrow, maxmod = 0, 0
      for  i = index + 1:size(smith, 1), j = index + 1: size(smith, 2)
        smith[i, j] % diag == 0 && continue
        if maxmod == 0 || abs(smith[i, j] % diag) > maxmod
          maxrow = i
          maxmod = abs(smith[i, j] % diag)
        end
      end
      if maxmod != 0
        smith[index, :] += smith[maxrow, :]
        left[index, :] += left[maxrow, :]
      end
    end
  end
  if smith[end, end] < 0
    smith[end, :] *= -1
    left[end, :] *= -1
  end
  left, smith, right
end
end

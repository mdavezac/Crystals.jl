module Gruber
export gruber

const default_tolerance = 1e-8
const max_no_change = 2

function no_opt_change_test(new, last)
  const m_new = 16e0 * new;
  const diff = new - last;
  const m_new_plus_diff = m_new + diff;
  const m_new_plus_diff_minus_m_new = m_new_plus_diff - m_new;
  m_new_plus_diff_minus_m_new != 0;
end

function def_test(real::Real; tolerance::Real=default_tolerance)
  real ≥ tolerance && return [0, 1]
  real > -tolerance && return [1, 0]
  return [0, 0]
end

function def_test(args; tolerance=default_tolerance)
  result = [0, 0]
  for u in args
    result += def_test(u, tolerance=tolerance)
  end
  result
end

function def_gt_0(args...; tolerance=default_tolerance)
  zero, positive = def_test(args; tolerance=tolerance)
  positive == 3 || (zero == 0 && positive == 1)
end


function n1_action(params::Vector, rinv::Matrix)
  rinv[:, :] = rinv * [0 -1 0; -1 0 0; 0 0 -1]
  params[1], params[2] = params[2], params[1]
  params[4], params[5] = params[5], params[4]
end

function n2_action(params::Vector, rinv::Matrix)
  rinv[:, :] = rinv * [-1 0 0; 0 0 -1; 0 -1 0]
  params[2], params[3] = params[3], params[2]
  params[5], params[6] = params[6], params[5]
end

function n3_action(params::Vector, rinv::Matrix; tolerance=default_tolerance)
  const i = params[4] ≤ -tolerance ? -1 : 1
  const j = params[5] ≤ -tolerance ? -1 : 1
  const k = params[6] ≤ -tolerance ? -1 : 1
  rinv[:, :] = rinv * [i 0 0; 0 j 0; 0 0 k]
  params[4:end] = abs(params[4:end])
end

function n4_action(params::Vector, rinv::Matrix; tolerance=default_tolerance)
  const i = params[4] ≥ tolerance ? -1 : 1
  const j = params[5] ≥ tolerance ? -1 : 1
  const k = params[6] ≥ tolerance ? -1 : 1
  update = diagm([i, j, k])
  if k == 1 && params[6] > -tolerance
    update[3, 3] = -1
  elseif j == 1 && params[5] > -tolerance 1
    update[3, 3] = -1
  elseif i == 1 && params[4] > -tolerance 0
    update[1, 1] = -1
  elseif i * j * k == -1
    error("Internal error")
  end
  rinv[:, :] = rinv * update
  params[4:end] = -abs(params[4:end])
end

function n5_action(params::Vector, rinv::Matrix; tolerance=default_tolerance)
  const sign = params[4] > tolerance ? -1 : 1
  rinv[:, :] = rinv * [1 0 0; 0 1 sign; 0 0 1]
  params[3] += params[2] + sign * params[4]
  params[4] += 2sign * params[2]
  params[5] += sign * params[6]
end

function n6_action(params::Vector, rinv::Matrix; tolerance=default_tolerance)
  const sign = params[5] > tolerance ? -1 : 1
  rinv[:, :] = rinv * [1 0 sign; 0 1 0; 0 0 1]
  params[3] += params[1] + sign * params[5]
  params[4] += sign * params[6]
  params[5] += 2sign * params[1]
end

function n7_action(params::Vector, rinv::Matrix; tolerance=default_tolerance)
  const sign = params[6] > tolerance ? -1 : 1
  rinv[:, :] = rinv * [1 sign 0; 0 1 0; 0 0 1]
  params[2] += params[1] + sign * params[6]
  params[4] += sign * params[5]
  params[6] += 2sign * params[1]
end

function n8_action(params::Vector, rinv::Matrix)
  rinv[:, :] = rinv * [1 0 1; 0 1 1; 0 0 1]
  params[3] += sum(params[1:2]) + sum(params[4:end])
  params[4] += 2params[2] + params[6]
  params[5] += 2params[1] + params[6]
end

""" Determines Gruber cell of an input cell

The Gruber cell is an optimal parameterization of a lattice, eg shortest
cell-vectors and angles closest to 90 degrees.

# Arguments
* `cell::Matrix`: the input lattice cell-vectors. Cannot be singular.
* `itermax::Int`: maximum number of iterations before bailing out
* `tolerance::Real`: tolerance parameter when comparing real numbers
"""
function gruber(cell::Matrix; tolerance=default_tolerance, itermax=50, max_no_change=10)
  size(cell, 1) == size(cell, 2) || error("Matrix not rectangular")
  abs(det(cell)) > tolerance || error("Singular matrix");
  const metric = transpose(cell) * cell
  params =
    vcat(diag(metric), [2metric[2, 3], 2metric[1, 3], 2metric[1, 2]])
  rinv = eye(size(metric, 1))
  nochange, previous = 0, -params[1:3]
  iteration::Int = 0
  for iteration in 1:itermax
    const ε = tolerance
    condition0 =
        (a, b, d, e) -> a ≥ b + ε || (abs(a - b) < ε && abs(d) ≥ abs(e) + ε)

    condition0(params[1], params[2], params[4], params[5]) &&
        n1_action(params, rinv)
    condition0(params[2], params[3], params[5], params[6]) &&
        (n2_action(params, rinv); continue)

    if def_gt_0(params[4:end]...; tolerance=tolerance)
      n3_action(params, rinv; tolerance=tolerance)
    else
      n4_action(params, rinv; tolerance=tolerance)
      if all(abs(previous - params[1:3]) .< tolerance)
        no_change += 1
      else
        no_change = 0
      end
      no_change < max_no_change || break
    end

    condition1 =
       (d, b, e, f) -> abs(d) ≥ b + ε ||
          (abs(d - b) < ε && 2e ≤ f - ε) ||
          (abs(d + b) < ε && f ≤ -ε)
    condition1(params[4], params[2], params[5], params[6]) &&
        (n5_action(params, rinv; tolerance=tolerance); continue)
    condition1(params[5], params[1], params[4], params[6]) &&
        (n6_action(params, rinv; tolerance=tolerance); continue)
    condition1(params[6], params[1], params[4], params[5]) &&
        (n7_action(params, rinv; tolerance=tolerance); continue)

    const sum_no_c = sum(params[1:2]) + sum(params[4:end])
    (
        sum_no_c ≤ -ε ||
        (abs(params[6] - params[1]) < ε && (2params[4] ≤ params[5] - ε)) ||
        (abs(sum_no_c) < ε  && 2params[1] + 2params[5] + params[6] ≥ ε)
    ) || break
    n8_action(params, rinv)
  end
  iteration == itermax && error("Reached end of iteration without converging")
  cell * rinv
end

end

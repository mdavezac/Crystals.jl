module Gruber
export gruber, niggly

using Crystals.Constants: default_tolerance
using Crystals.Utilities: cell_parameters
using Crystals.Structures: Crystal, volume
using Crystals: Log
using Unitful: unit, ustrip, Quantity

function no_opt_change_test(new, last)
    const m_new = 16e0 * new;
    const difference = new - last;
    const m_new_plus_diff = m_new + difference;
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
    z₀, positive = def_test(args; tolerance=tolerance)
    positive == 3 || (z₀ == 0 && positive == 1)
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
    params[4:end] = abs.(params[4:end])
end

function n4_action(params::Vector, rinv::Matrix; tolerance=default_tolerance)
    const i = params[4] ≥ tolerance ? -1 : 1
    const j = params[5] ≥ tolerance ? -1 : 1
    const k = params[6] ≥ tolerance ? -1 : 1
    update = diagm([i, j, k])
    if i * j * k < 0
        if k == 1 && params[6] > -tolerance
            update[3, 3] = -1
        elseif j == 1 && params[5] > -tolerance
            update[3, 3] = -1
        elseif i == 1 && params[4] > -tolerance
            update[1, 1] = -1
        elseif i * j * k == -1
            Log.error("Internal error")
        end
    end
    rinv[:, :] = rinv * update
    params[4:end] = -abs.(params[4:end])
end

function n5_action(params::Vector, rinv::Matrix; tolerance=default_tolerance)
    const pos_or_neg = params[4] > tolerance ? -1 : 1
    rinv[:, :] = rinv * [1 0 0; 0 1 pos_or_neg; 0 0 1]
    params[3] += params[2] + pos_or_neg * params[4]
    params[4] += 2pos_or_neg * params[2]
    params[5] += pos_or_neg * params[6]
end

function n6_action(params::Vector, rinv::Matrix; tolerance=default_tolerance)
    const pos_or_neg = params[5] > tolerance ? -1 : 1
    rinv[:, :] = rinv * [1 0 pos_or_neg; 0 1 0; 0 0 1]
    params[3] += params[1] + pos_or_neg * params[5]
    params[4] += pos_or_neg * params[6]
    params[5] += 2pos_or_neg * params[1]
end

function n7_action(params::Vector, rinv::Matrix; tolerance=default_tolerance)
    const pos_or_neg = params[6] > tolerance ? -1 : 1
    rinv[:, :] = rinv * [1 pos_or_neg 0; 0 1 0; 0 0 1]
    params[2] += params[1] + pos_or_neg * params[6]
    params[4] += pos_or_neg * params[5]
    params[6] += 2pos_or_neg * params[1]
end

function n8_action(params::Vector, rinv::Matrix)
    rinv[:, :] = rinv * [1 0 1; 0 1 1; 0 0 1]
    params[3] += sum(params[1:2]) + sum(params[4:end])
    params[4] += 2params[2] + params[6]
    params[5] += 2params[1] + params[6]
end

"""
    gruber(cell::Matrix;
           tolerance::Real=default_tolerance, itermax::Unsigned=50,
           max_no_change::Unsigned=10)

Determines Gruber cell of an input cell.

The Gruber cell is an optimal parameterization of a lattice, e.g. shortest
cell-vectors and angles closest to 90 degrees. The new cell is in the same basis
as the origin cell: no rotation has been incurred. The cell parameters are
uniquely determined, even though the cell itself is not (certain symmetry
operations leaving the structure unchanged may yield a more recognizable cell).
If you want a unique Cartesian cell (in a different Cartesian basis), use
the `niggly` algorithm.

# Arguments
* `cell::Matrix`: the input lattice cell-vectors. Cannot be singular.
* `itermax::Integer`: maximum number of iterations before bailing out
* `tolerance::Number`: tolerance parameter when comparing real numbers
* `max_no_change::Integer`: Maximum number of times to go through algorithm
  without changes before bailing out
"""
function gruber{T <: Number}(cell::Matrix{T};
                             tolerance::Real=default_tolerance, itermax::Integer=50,
                             max_no_change::Integer=10)
    size(cell, 1) == size(cell, 2) || Log.error("Matrix not rectangular")
    volume(cell) > tolerance || Log.error("Singular matrix");
    if itermax ≤ 0
        itermax = typemax(itermax)
    end
    if max_no_change ≤ 0
        max_no_change = typemax(max_no_change)
    end

    const ε = tolerance
    const metric = transpose(cell) * cell
    params = vcat(diag(metric), [2metric[2, 3], 2metric[1, 3], 2metric[1, 2]])
    rinv = eye(size(metric, 1))
    no_change, previous = 0, -params[1:3]
    iteration::Int = 0
    for iteration in 1:itermax
        condition0 =
            (a, b, d, e) -> a ≥ b + ε || (abs(a - b) < ε && abs(d) ≥ abs(e) + ε)

        condition0(params[1], params[2], params[4], params[5]) &&
            n1_action(params, rinv)
        condition0(params[2], params[3], params[5], params[6]) &&
            (n2_action(params, rinv); continue)

        if def_gt_0(params[4:end]...; tolerance=ε)
            n3_action(params, rinv; tolerance=ε)
        else
            n4_action(params, rinv; tolerance=ε)
            if all(abs.(previous .- params[1:3]) .< ε)
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
            (n5_action(params, rinv; tolerance=ε); continue)
        condition1(params[5], params[1], params[4], params[6]) &&
            (n6_action(params, rinv; tolerance=ε); continue)
        condition1(params[6], params[1], params[4], params[5]) &&
            (n7_action(params, rinv; tolerance=ε); continue)

        const sum_no_c = sum(params[1:2]) + sum(params[4:end])
        (
            sum_no_c ≤ -ε ||
            (abs(params[6] - params[1]) < ε && (2params[4] ≤ params[5] - ε)) ||
            (abs(sum_no_c) < ε  && 2params[1] + 2params[5] + params[6] ≥ ε)
        ) || break
        n8_action(params, rinv)
    end
    iteration == itermax &&
        Log.error("Reached end of iteration without converging")
    cell * rinv
end

function gruber{T, D, U}(cell::Matrix{Quantity{T, D, U}};
                         tolerance::Real=default_tolerance, itermax::Integer=50,
                         max_no_change::Integer=10)
    gruber(ustrip(cell);
           tolerance=tolerance,
           itermax=itermax,
           max_no_change=max_no_change) * unit(Quantity{T, D, U})
end

"""
    niggly(cell::Matrix; kwargs...)

Determines a unique Cartesian cell equivalent to the input, with the shortest
possible vectors and squarest angles. For an explanation of the parameters, see
`gruber`. In practice, this function computes the cell-parameters of a `gruber` cell and
then reconstructs the cell matrix. Hence, the result should be quite unique for any lattice
representation, including any rotation of the underlying Cartesian coordinates.
"""
niggly(cell::Matrix, kwargs...) = cell_parameters(cell_parameters(gruber(cell; kwargs...)))

end

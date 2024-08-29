import JSON

# Some enumerations and correspondance dicts
# Dict comprehension isn't defined in 0.4 and isn't handled by Compat
# See https://github.com/JuliaLang/Compat.jl/issues/231
const objtypes = ["none", "constant", "linear", "quadratic", "sum_of_squares", "other"]
const classdb_objtype = Dict(
  "N" => objtypes[1],
  "C" => objtypes[2],
  "L" => objtypes[3],
  "Q" => objtypes[4],
  "S" => objtypes[5],
  "O" => objtypes[6],
)

const contypes = ["unc", "fixed_vars", "bounds", "network", "linear", "quadratic", "general"]
const classdb_contype = Dict(
  "U" => contypes[1],
  "X" => contypes[2],
  "B" => contypes[3],
  "N" => contypes[4],
  "L" => contypes[5],
  "Q" => contypes[6],
  "O" => contypes[7],
)

const origins = ["academic", "modelling", "real"]
const classdb_origin = Dict("A" => origins[1], "M" => origins[2], "R" => origins[3])

"""
    select(; min_var=1, max_var=Inf, min_con=0, max_con=Inf,
             objtype=*, contype=*,
             only_free_var=false, only_bnd_var=false,
             only_linear_con=false, only_nonlinear_con=false,
             only_equ_con=false, only_ineq_con=false,
             custom_filter=*)

Returns a subset of the CUTEst problems based on the classification file `classf.json`.
This file is included with the package.

## Arguments

- `min_var` and `max_var`: Define the range for the number of variables in the problem.
- `min_con` and `max_con`: Define the range for the number of constraints in the problem (e.g., use `max_con=0` for unconstrained problems or `min_con=1` for constrained problems).
- `only_*` flags:
  - `only_free_var`: Include only problems with free variables.
  - `only_bnd_var`: Include only problems with bounded variables.
  - `only_linear_con`: Include only problems with linear constraints.
  - `only_nonlinear_con`: Include only problems with nonlinear constraints.
  - `only_equ_con`: Include only problems with equality constraints.
  - `only_ineq_con`: Include only problems with inequality constraints.
  Note: These flags are mutually exclusive; only one can be true at a time, or both can be false.

- `objtype`: Classification of the objective function according to the [MASTSIF classification file](https://www.cuter.rl.ac.uk/Problems/classification.shtml). It can be a number, a symbol, a string, or an array of those:
  - `1`, `:none` or `"none"`: No objective function.
  - `2`, `:constant` or `"constant"`: Objective function is a constant.
  - `3`, `:linear` or `"linear"`: Objective function is linear.
  - `4`, `:quadratic` or `"quadratic"`: Objective function is quadratic.
  - `5`, `:sum_of_squares` or `"sum_of_squares"`: Objective function is a sum of squares.
  - `6`, `:other` or `"other"`: Objective function does not fit the above categories.

- `contype`: Classification of the constraints according to the same MASTSIF classification file:
  - `1`, `:unc` or `"unc"`: No constraints.
  - `2`, `:fixed_vars` or `"fixed_vars"`: Only fixed variables.
  - `3`, `:bounds` or `"bounds"`: Only bounded variables.
  - `4`, `:network` or `"network"`: Constraints represent a network adjacency matrix.
  - `5`, `:linear` or `"linear"`: Linear constraints.
  - `6`, `:quadratic` or `"quadratic"`: Quadratic constraints.
  - `7`, `:other` or `"other"`: More general constraints.

- `custom_filter`: A function to apply additional filtering to the problem data, which is a dictionary with the following fields:
  - `"objtype"`: String representing the objective function type.
  - `"contype"`: String representing the constraint type.
  - `"regular"`: Boolean indicating whether the problem is regular.
  - `"derivative_order"`: Integer for the highest derivative order available.
  - `"origin"`: String indicating the origin of the problem: `"academic"`, `"modelling"`, or `"real"`.
  - `"has_interval_var"`: Boolean indicating if there are interval variables.
  - `"variables"`: Dictionary with fields related to variables.
  - `"constraints"`: Dictionary with fields related to constraints.

To select only problems with a fixed number of variables, use:

```julia
custom_filter = x -> x["variables"]["can_choose"] == false
```
"""
function select(;
  min_var = 1,
  max_var = Inf,
  min_con = 0,
  max_con = Inf,
  objtype = objtypes,
  contype = contypes,
  only_free_var = false,
  only_bnd_var = false,
  only_linear_con = false,
  only_nonlinear_con = false,
  only_equ_con = false,
  only_ineq_con = false,
  custom_filter::Function = x -> true,
)
  # Checks for conflicting option
  @assert !only_free_var || !only_bnd_var
  @assert !only_linear_con || !only_nonlinear_con
  @assert !only_equ_con || !only_ineq_con

  objtype = canonicalize_ftype(objtype, objtypes)
  contype = canonicalize_ftype(contype, contypes)

  if !(objtype ⊆ objtypes)
    error("objtypes $objtype not supported")
  end
  if !(contype ⊆ contypes)
    error("contypes $contype not supported")
  end

  data = JSON.parsefile(joinpath(dirname(@__FILE__), "classf.json"))
  problems = keys(data)
  selection = Vector{String}()
  for p in problems
    pv = data[p]["variables"]
    pc = data[p]["constraints"]
    nvar, ncon = pv["number"], pc["number"]
    if nvar < min_var ||
       nvar > max_var ||
       ncon < min_con ||
       ncon > max_con ||
       (only_free_var && pv["free"] < nvar) ||
       (only_bnd_var && pv["free"] > 0) ||
       (only_linear_con && pc["linear"] < ncon) ||
       (only_nonlinear_con && pc["linear"] > 0) ||
       (only_equ_con && pc["equality"] < ncon) ||
       (only_ineq_con && pc["equality"] > 0) ||
       !(data[p]["objtype"] in objtype) ||
       !(data[p]["contype"] in contype) ||
       !(custom_filter(data[p]))
      continue
    end
    push!(selection, p)
  end
  return selection
end

canonicalize_ftype(reqtype::Integer, allowedtypes) = [allowedtypes[reqtype]]
canonicalize_ftype(reqtype::Symbol, allowedtypes) = [string(reqtype)]
canonicalize_ftype(reqtype::AbstractString, allowedtypes) = [reqtype]
canonicalize_ftype(reqtype::AbstractVector{T}, allowedtypes) where {T <: Integer} =
  allowedtypes[reqtype]
canonicalize_ftype(reqtype::AbstractVector{Symbol}, allowedtypes) = map(string, reqtype)
canonicalize_ftype(reqtype, allowedtypes) = reqtype

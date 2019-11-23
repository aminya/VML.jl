using VML
using Distributions, Plots, BenchmarkTools

include(joinpath(dirname(dirname(@__FILE__)), "test", "common.jl"))
complex = !isempty(ARGS) && ARGS[1] == "complex"


"""
    bench(fns, input)

benchmark function for VML.jl. Calls both Base and VML functions and stores the benchmarks in two nested Dict. First layer specifies type, and second layer specifies the function name. The result is a Tuple, 1st element being benchmark for Base/SpecialFunctions and 2nd element being for VML.

# Examples
```julia
times = bench(fns, input)

times[Float64][:acos][1] # Base.acos benchmark for Float64
times[Float64][:acos][2] # VML.acos benchmark for Float64
```
"""
function bench(fns, input)
    Dict(t => begin
        Dict( fn[2] => begin
            base_fn = eval(:($(fn[1]).$(fn[2])))
            vml_fn = eval(:(VML.$(fn[2])))
            println("benchmarking $vml_fn")
            timesBase = @benchmark $base_fn.($inp...)
            timesVML = @benchmark $vml_fn($inp...)

            timesBase, timesVML
        end for (fn, inp) in zip(fns, input[t]) )
    end for t in types)
end

function ratioci(y, x, alpha = 0.05)
    tq² = abs2(quantile(TDist(length(x) + length(y) - 2), alpha))
    μx = mean(x)
    σx² = varm(x, μx)
    μy = mean(y)
    σy² = varm(y, μy)
    a = sqrt((μx * μy)^2 - (μx^2 - tq² * σx²) * (μy^2 - tq² * σy²))
    b = μx^2 - tq² * σx²
    (((μx * μy) - a) / b, ((μx * μy) + a) / b)
end

# First generate some random data and test functions in Base on it
const NVALS = 1_000_000
base_unary = complex ? base_unary_complex : base_unary_real
base_binary = complex ? base_binary_complex : base_binary_real
types = complex ? (Complex64, Complex128) : (Float32, Float64)

input = Dict( t =>
[
 [(randindomain(t, NVALS, domain),) for (_, _, domain) in base_unary];
 [(randindomain(t, NVALS, domain1), randindomain(t, NVALS, domain2)) for (_, _, domain1, domain2) in base_binary];
 (randindomain(t, NVALS, (0, 100)), randindomain(t, 1, (-1, 20))[1])
]
    for t in types)

fns = [[x[1:2] for x in base_unary_real];
       [x[1:2] for x in base_binary_real]]


# do benchmark
builtint, vmlt = bench(fns, input)


# Print ratio
clf()
colors = ["r", "y"]
for itype = 1:length(types)
    builtint = builtin[types[itype]]
    vmlt = vml[types[itype]]
    μ = vec(map(mean, builtint)./map(mean, vmlt))
    ci = zeros(Float64, 2, length(fns))
    for ifn = 1:length(builtint)
        lower, upper = ratioci(builtint[ifn], vmlt[ifn])
        ci[1, ifn] = μ[ifn] - lower
        ci[2, ifn] = upper - μ[ifn]
    end
    bar(0.2+(0.4*itype):length(fns), μ, 0.4, yerr=ci, color=colors[itype], ecolor="k")
end
ax = gca()
ax[:set_xlim](0, length(fns)+1)
fname = [string(fn.env.name) for fn in fns]
if !complex
    fname[end-1] = "A.^B"
    fname[end] = "A.^b"
end
xticks(1:length(fns)+1, fname, rotation=70, fontsize=10)
title("VML Performance")
ylabel("Relative Speed (Base/VML)")
legend([string(x) for x in types])
ax[:axhline](1; color="black", linestyle="--")
savefig("performance$(complex ? "_complex" : "").png")

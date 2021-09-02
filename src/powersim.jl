using Random
using Random: default_rng
using .Threads: nthreads

using Distributions, HypothesisTests, StatsBase, ProgressMeter, LoopVectorization,
    ThreadPools, GaitSymmetry, SafeBuffers
using GaitSymmetry: Ratio

export gensimdata, calcsymfun, estimatepower, estimatepower!, powersim

function gensimdata(rng, xdist, ydist, samples, T=promote_type(eltype(xdist), eltype(ydist)))
    stor = Matrix{T}(undef, (samples, 2))
    _gensimdata!(rng, stor, xdist, ydist)
end

function _gensimdata!(rng, stor, xdist, ydist)
    randn!(rng, stor)
    @views _randn_kernel!(stor[:,1], mean(xdist), std(xdist))
    @views _randn_kernel!(stor[:,2], mean(ydist), std(ydist))

    return stor
end

function _randn_kernel!(arr, μ, σ)
    @avx for i in eachindex(arr)
        arr[i] = abs(fma(arr[i], σ, μ))
    end
end

function calcsymfun(symfun, stor::AbstractArray{T}, len=size(stor, 1)) where T
    size(stor, 2) === 2 || throw(DimensionMismatch("expected 2 columns; got $(size(stor, 2))"))
    len ≤ size(stor, 1) || throw(DimensionMismatch("`len` must be ≤ `size(stor, 1)`"))
    res = Vector{T}(undef, len)
    _calcsymfun!(res, symfun, len, stor)
end

function _calcsymfun!(res, symfun, len, stor)
    @inbounds for i in 1:len-1
        x = stor[i,1]
        x′ = stor[i+1,1]
        y = stor[i,2]

        res[i] = symfun(x, x′) - symfun(x, y)
    end
    res[end] = symfun(stor[end,1], stor[1,1]) - symfun(stor[end,1], stor[end,2])

    return res
end

function estimatepower(symfun, xdist, ydist, n, samples)
    stor = gensimdata(default_rng(), xdist, ydist, n*samples)
    res = calcsymfun(symfun, stor)
    _estimatepower(res, n, samples)
end

function estimatepower(symfun, stor, n, samples)
    res = calcsymfun(symfun, stor, n*samples)
    _estimatepower(res, n, samples)
end

function estimatepower!(res, symfun, stor, n, samples)
    fullN = n*samples
    _calcsymfun!(res, symfun, fullN, stor)
    _estimatepower(res, n, samples)
end

function _estimatepower(res, n, samples)
    numsig = _numsignificant(res, n, samples)
    return numsig# / samples
end

function _numsignificant(res, n, samples, α=0.05)
    numsig = 0
    for i in 0:(samples-1)
        @inbounds slice = @view(res[(1:n).+i])
        avg, sd = mean_and_std(slice; corrected=true)
        if pvalue(OneSampleTTest(avg, sd, n)) ≤ α
            numsig += 1
        end
    end

    return numsig
end

function powersim(
    metrics,
    ratios,
    vars,
    Ns,
    samples,
    U=Float32;
    null=false,
    batch=1
)
    p = Progress(length(metrics)*length(ratios)*length(vars)*sum(Ns .* samples)+1;
        dt=5, desc="Batch $(lpad(batch, 2)) progress: ")

    ratios = convert(Array{U}, ratios)
    data = Array{U}(undef, (length(Ns), length(metrics), length(ratios), length(vars)))

    fullN = maximum(Ns)*samples
    respool = BufferPool(nthreads()+1, () -> Vector{U}(undef, fullN))
    storpool = BufferPool(nthreads()+1, () -> Matrix{U}(undef, (fullN, 2)))

    idxs = [ (σi, i) for σi in eachindex(vars) for i in eachindex(ratios) ]
    @qthreads for (σi, i) in idxs
        σ = vars[σi]
        xdist = Normal(1, σ)
        μ = ratios[i]
        μ_sd = μ*σ
        ydist = ifelse(null, xdist, Normal(μ, μ_sd))

        withbuffer(storpool) do stor
            _gensimdata!(default_rng(), stor, xdist, ydist)
            for mi in eachindex(metrics)
                if metrics[mi] === Alv20b
                    metric = Alv20b(U(σ))
                else
                    metric = metrics[mi]
                end

                withbuffer(respool) do res
                    _calcsymfun!(res, metric, fullN, stor)
                    for ni in eachindex(Ns)
                        data[ni, mi, i, σi] = _estimatepower(res, Ns[ni], samples)
                        next!(p; step=Ns[ni]*samples)
                    end
                end
            end
        end
    end
    finish!(p)

    return data
end


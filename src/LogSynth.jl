module LogSynth

using Counters, Random

export SkipListDistribution, AliasTableDistribution

md"""
A `SkipListDistribution` provides an implementation of a multinomial
distribution that has ``O(log(n))`` sample time, but which allows the
underlying probability for any element to be adjusted in ``O(log(n))`` 
time and which supports appending new values to the distribution at 
any time (also in ``O(log(n))`` time).

This is done by keeping a skip-list that contains the cumulative sum 
of a sequence of values ``p_i``. This skip-list can be searched with a
value ``u \in [0, \sum_i p_i)`` to find the largest index ``m`` such 
that ``\sum_{i=1}^m p_i < u``. This is a little bit different from
an ordinary skip list in that each level of the skip list only keeps
the difference to the next value at the same level. By keeping the 
difference instead of the value, we can increase a weight ``p[i]`` 
by increasing a limited number of differences. 

Typically, the weights ``p[i]`` are counts, but that isn't a requirement.

Traditionally, a skip list associates a variable sized vector of skip distances
with each value in the skip list (values are stored in ascending order). The 
length of the skip vector is chosen randomly so short skip vectors are much 
more common. The skip list has a single skip vector that has offsets to the 
first entry at each level. To find a value, we start at the highest level 
in the skip list that is not past the desired value and move down one level
at a time.

The skip list here is done a bit differently. Instead of keeping a vector of
offsets to the next element at different levels (a vertical arrangement), we
keep a vector for each level which contains the difference in value to the next
value (a horizontal arrangement). Advantages of this alternative is that it is
very easy to append to some of the lists and we avoid the use of any pointers.

A common use of this is to store an empirical cumulative distribution. In 
particular, this `SkipListDistribution` is a very efficient for implementing
a Pitman-Yor process which samples from a notionally infinite discrete 
distribution of the integers.
"""
struct SkipListDistribution{Weight <: AbstractFloat, Index <: Integer}
    height::Int32
    weight::Vector{Vector{Weight}}
    child::Vector{Vector{Index}}

    α::Float64
    discount::Float64

    skipProbability::Float64
end

function SkipListDistribution{Weight, Index}(height; alpha=0.0, discount=0.0, skipProbability=0.25) 
    where {Weight <: AbstractFloat, Index <: Integer} 
    init(
        SkipListDistribution{Weight, Index}(height, Vector{Weight}(), Vector{Index}(), 
                                            alpha, discount, skipProbability))
end

function SkipListDistribution(height = 10; alpha=0.0, discount=0.0, skipProbability=0.25) 
    init(
        SkipListDistribution{Float64, Int32}(height, Vector{Float64}(), Vector{Int32}(),
                                             alpha, discount, skipProbability))
end

function init(self::SkipListDistribution{Weight, Index}) where {Weight, Index} 
    for i in 1:self.height
        push!(self.weight, Vector{Weight}())
        push!(self.child, Vector{Index}())
    end
    return self 
end

md"""
Given a probability as an index, return a vector of indexes into the 
skip list
"""
function getindexvector(dist::SkipListDistribution{Weight, Index}, p::Float64) where {Weight <: AbstractFloat, Index <: Int}
    0 ≤ weight ≤ 1  || error("Weight should be in [0,1]")
    n = p * total(dist)
    r = Vector{Index}()
    j = 1
    w = 0.0
    for level = 1:dist.height
        let v = dist.list[level]
            while j <= length(v) && w + v[j].weight < n
                j += 1
                w += v[j].weight
            end
            push!(r, j)
            j = v[j].firstChild
        end
    end
    return r
end

import Base: length
Base.length(dist::SkipListDistribution) = Base.length(dist.weight[end])

function getindexvector(dist::SkipListDistribution{Weight, Index}, i::Integer) where {Weight <: AbstractFloat, Index <: Integer}
    1 ≤ i ≤ length(dist) || error("Index out of bounds. Expected i ∈ [1, $(length(dist))]")
    r = zeros(Index, dist.height)
    j = i
    for level = dist.height:-1:1
        k = 1
        while k < length(dist.child[level]) && dist.child[level][k+1] <= i
            k += 1
        end
        r[level] = k
    end
    return r
end

total(dist::SkipListDistribution{Weight, Index}) where {Weight <: AbstractFloat, Index <: Integer} = sum(dist.weight[1])

function getindexvector(dist::SkipListDistribution{Weight, Index}, p::Float64) 
    where {Weight <: AbstractFloat, Index <: Integer}
    0 ≤ p ≤ 1 || error("Sample probability must be in [0, 1] but was $p")
    p = p * total(dist)
    r = zeros(Index, dist.height)
    w = 0.0
    j = 1
    for level = 1:dist.height
        while j < length(dist.weight[level]) && w + dist.weight[level][j] ≤ p
            w += dist.weight[level][j]
            j += 1
        end
        r[level] = j
        j = dist.child[level][j]
    end
    return r
end

function getindex(dist::SkipListDistribution{Weight, Index}, p::Float64) 
    where {Weight <: AbstractFloat, Index <: Integer} 
    getindexvector(dist, p)[end]
end

function inc!(dist::SkipListDistribution{Weight, Index}, p::Float64, Δw::Weight) 
    where {Weight <: AbstractFloat, Index <: Integer}
    index = getindexvector(dist, p)
    for (level,i) in enumerate(index)
        dist.weight[level][i] += Δw
    end
    return index[end]
end

import Base.push!
function Base.push!(dist::SkipListDistribution{Weight, Index}, w::Weight) 
    where {Weight <: AbstractFloat, Index <: Integer}
    if length(dist) == 0
        for level in 1:dist.height
            push!(dist.child[level], 1)
            push!(dist.weight[level], w)
        end
    else
        # add the new entry to a subset of levels with zero weight
        child = length(dist) + 1
        for level in dist.height:-1:1
            push!(dist.weight[level], 0)
            push!(dist.child[level], child)
            child = length(dist.child[level])
            # probably break out
            rand() < dist.skipProbability || break
        end

        # and then increment all last entries (including old and new)
        for level in 1:dist.height
            dist.weight[level][end] += w
        end
    end
    dist
end

function draw_value(rng, dist::SkipListDistribution{Weight, Index})::Index where {Weight, Index}
    norm = dist.α + total(dist)
    newtable = (dist.α + dist.discount * length(dist)) / norm
    u = rand(rng)
    if u < newtable
        push!(dist, 1 - dist.discount)
        return length(dist)
    else
        return inc!(dist, u / (1-newtable), 1.0)
    end
end
                                

md"""
An `AliasTableDistribution` provides an implementation of a multinomial
distribution that is good for O(1) sampling from static distributions.

If the distribution changes often, a `SkipListDistribution` may be a 
better option.
"""
struct AliasTableDistribution{T}
    u::Vector{Float64}
    k::Vector{Int32}
    v::Vector{T}

    cnt::Counter{T}

    AliasTableDistribution{T}(cnt::Counter{T}) where T = begin
        u, k, v = build_alias_table(cnt)
        new{T}(u, k, v, cnt)
    end
end

build_alias_table(cnt::Counter{T}) where T = begin
    p = float.(values(cnt))
    v = collect(keys(cnt))
    p = p ./ sum(p)
    n = length(v)
    u = p .* n
    k = collect(1:length(p))

    ϵ = 1.0 / n
    z = sortperm(p)
    little = [z[i] for i in 1:n if p[z[i]] < ϵ]
    big = [z[i] for i in n:-1:1 if p[z[i]] > ϵ]
    i = 1
    j = 1
    while i ≤ length(little) && j ≤ length(big)
        deficit = 1 - u[little[i]]
        k[little[i]] = big[j]
        i += 1
        u[big[j]] -= deficit
        if u[big[j]] < 1
            push!(little, big[j])
            j += 1
        elseif u[big[j]] == 1
            j += 1
        end
    end
    return u, k, v
end

draw_value(rng::AbstractRNG, dist::AliasTableDistribution{T}) where T = begin
    n = length(dist.u)
    z = n * rand(rng) + 1
    i = Int(floor(z))
    z = z - i

    if z < dist.u[i]
        return dist.v[i]
    else
        return dist.v[dist.k[i]]
    end
end


end # module

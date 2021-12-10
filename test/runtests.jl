using Test, LogSynth

@test length(SkipListDistribution(5)) == 0

let s = SkipListDistribution(5)
    push!(s, 1.0)
    push!(s, 2.0)
    push!(s, 1.0)
    @test length(s) == 3
    @test LogSynth.total(s) == 4.0
    @test LogSynth.getindexvector(s, 2)[end] == 2
    @test LogSynth.getindexvector(s, 0.24999)[end] == 1
    @test LogSynth.getindexvector(s, 0.25)[end] == 2
    @test LogSynth.getindexvector(s, 0.5)[end] == 2
    @test LogSynth.getindexvector(s, 0.75)[end] == 3

    @test LogSynth.inc!(s, 0.5, 1.0) == 2
    @test length(s) == 3
    @test LogSynth.total(s) == 5.0
    @test LogSynth.getindexvector(s, 0.19999)[end] == 1
    @test LogSynth.getindexvector(s, 0.2)[end] == 2
    @test LogSynth.getindexvector(s, 0.7999)[end] == 2
    @test LogSynth.getindexvector(s, 0.8)[end] == 3
end

let s = SkipListDistribution(5)
    for i in 1:100_000
        push!(s, 1.0)
    end

    # verify that the lengths of the different levels are about right
    # 5σ gives us about 1.5 x 10^-6 chance of random failure
    p = 1.0
    for i in s.height:-1:1
        σ = sqrt(p * (1-p) * length(s))
        @test length(s.child[i]) ≈ length(s) * p atol = 5σ
        p *= s.skipProbability
    end
end




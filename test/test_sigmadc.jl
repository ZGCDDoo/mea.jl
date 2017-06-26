using Mea.SigmaDC
using Mea.Periodize

using Base.Test


@testset "Testing SigmaDC.jl" begin

    modelvec=Periodize.buildmodelvec("./data/self_ctow.dat", "./data/statsparams0.json")
    @test isapprox(SigmaDC.calc_sigmadc(modelvec, 12.0), 0.7637, atol=1e-3)

end  # testset

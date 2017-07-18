using Mea.Transport
using Mea.Periodize

using Base.Test


@testset "Testing Transport.jl" begin


    #setup
        
        cutoff_test = 20.0
        beta = 12.0
        sigmadc_good = 0.7637186262560831
        l11_good = sigmadc_good/beta

    # cutoff_test = 3.0
    # beta = 12.0
    # sigmadc_good = 0.48994580280943345
    # l11_good = sigmadc_good/beta

    @testset "test dfdw" begin
        @test isapprox(Transport.dfdw(12.0, 0.2), -0.91506 , atol=1.0e-4)
        @test isapprox(Transport.dfdw(12.0, -0.1), -2.13473, atol=1.0e-4)
    end

    
    @testset "test coefstrans" begin
        modelvec=Periodize.buildmodelvec("./data/self_ctow0.dat", "./data/statsparams0.json")
        transport = Transport.coefstrans(modelvec, beta, cutoff=cutoff_test, fout_name="dostest.dat", maxevals=80000)        
        @test isapprox(transport["n"], 0.495, atol=1e-3, rtol=1e-3)
        @test isapprox(transport["sigmadc"], sigmadc_good, atol=1.0e-5, rtol=1.0e-5)
        @test isapprox(transport["l11"], sigmadc_good/beta, atol=1.0e-5, rtol=1.0e-5)
        #@test isapprox(transport["l22"], )
        #@test isapprox(transport["l22"],)
        #@test isapprox(transport["seebeck"],)
    end

end  # testset

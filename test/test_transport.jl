using Mea.Transport
using Mea.Periodize

using Base.Test


@testset "Testing Transport.jl" begin


    #setup
        
        # cutoff_test = 20.0
        # beta = 12.0
        # sigmadc_good = 0.7637186262560831
        # l11_good = sigmadc_good/beta

    cutoff_test = 3.0
    beta = 12.0
    sigmadc_good = 0.48994580280943345
    l11_good = sigmadc_good/beta

    @testset "test dfdw" begin
        @test isapprox(Transport.dfdw(12.0, 0.2), -0.91506 , atol=1.0e-4)
        @test isapprox(Transport.dfdw(12.0, -0.1), -2.13473, atol=1.0e-4)
    end

    
    #@testset "test calc_l11" begin
    #    modelvec=Periodize.buildmodelvec("./data/self_ctow.dat", "./data/statsparams0.json")
    #    @test isapprox(Transport.calc_l11(modelvec, beta), l11_good, atol=1.0e-5, rtol=1.0e-5)
    #end

    
    # @testset "test calc_l21" begin
    # end


    # @testset "test calc_l22" begin
    # end


    @testset "test sigmadc" begin
        modelvec=Periodize.buildmodelvec("./data/self_ctow.dat", "./data/statsparams0.json")
        @test isapprox(Transport.calc_sigmadc(modelvec, beta, cutoff=cutoff_test), sigmadc_good, atol=1.0e-5, rtol=1.0e-5)
    end

end  # testset

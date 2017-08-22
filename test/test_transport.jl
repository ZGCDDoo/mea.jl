using Mea.Transport
using Mea.Periodize

using Base.Test


@testset "Testing Transport.jl" begin


    #setup
        
        cutoff_test = 200.0
        beta_ = 12.0
        sigmadc_good = 0.7637186262560831
        l11_good = sigmadc_good/beta_

    # cutoff_test = 3.0
    # beta_ = 12.0
    # sigmadc_good = 0.48994580280943345
    # l11_good = sigmadc_good/beta_

    @testset "test dfdw" begin
        @test isapprox(Transport.dfdw(12.0, 0.2), -0.91506 , atol=1.0e-4)
        @test isapprox(Transport.dfdw(12.0, -0.1), -2.13473, atol=1.0e-4)
    end

    
    @testset "test coefstrans" begin
        modelvec=Periodize.buildmodelvec("./data/self_ctow0.dat", "./data/statsparams0.json")
        transport = Transport.coefstrans(modelvec, beta_, cutoff=cutoff_test, fout_name="dos_green.dat", maxevals=80000, fctper="make_akw2green")
        transport_cum = Transport.coefstrans(modelvec, beta_, cutoff=cutoff_test, fout_name="dos_cum.dat", maxevals=80000, fctper="make_akw2cum")
        transport_trace = Transport.coefstrans(modelvec, beta_, cutoff=cutoff_test, fout_name="dos_trace.dat", maxevals=80000, fctper="make_akw2trace")    
        transport_cuba = Transport.coefstrans(modelvec, beta_, cutoff=cutoff_test, fout_name="dos_cuba.dat", maxevals=80000, libintegrator="cuba", fctper="make_akw2green")            
        println("transport_cuba = ", transport_cuba)
        println("transport = ", transport)
        println("transport_trace = ", transport_trace)
        println("transport_cum = ", transport_cum)
        @test isapprox(transport["n"], 0.495, atol=1e-2, rtol=1e-2)
        @test isapprox(transport["sigmadc"], sigmadc_good, atol=1.0e-5, rtol=1.0e-5)
        @test isapprox(transport_cuba["sigmadc"], sigmadc_good, atol=1.0e-3, rtol=1.0e-3)
        @test isapprox(transport["l11"], sigmadc_good/beta_, atol=1.0e-5, rtol=1.0e-5)
        @test isapprox(transport_cuba["l11"], sigmadc_good/beta_, atol=1.0e-3, rtol=1.0e-3)
        # @test isapprox(transport["l22"], )
        # @test isapprox(transport["l22"],)
        # @test isapprox(transport["seebeck"],)
    end

end  # testset

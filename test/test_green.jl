using Mea.Green
using Base.Test

 @testset "Testing Green.jl" begin

    # setup
    fin_gf_to = "./data/self_short_moy.dat"
    (zn, gf_c) = Green.readgreen_c(fin_gf_to, zn_col=1)
    gf_t = [(7.29148945e+00 -1.25688297e+00im) (4.31791692e+00 -7.37366863e-01im) (-4.37618486e+00 +7.44382438e-01im) (6.02348575e+00 -9.90968651e-01im) (4.05606418e+00 -4.07146712e-01im);
            (6.37534240e+00 -3.08656815e+00im) (3.92112500e+00 -1.86394356e+00im) (-3.80134014e+00 +1.78496589e+00im) (5.15813767e+00 -2.33078630e+00im) (3.66679027e+00 -9.42835945e-01im)]

    @testset "test_readgreen_c" begin:
        """ """
        zn_test = [5.235989e-002, 1.570800e-001]
        @test isapprox(zn, zn_test, atol=1e-6)

        gf_c_test = Green.t_to_c(zn, gf_t)
        gf_t_test = Green.c_to_t(zn, gf_c_test)

        @test isapprox(gf_c, gf_c_test)
        @test isapprox(gf_t, gf_t_test)

    end


    @testset "test_c" begin

        gf_c_test = zeros(Complex{Float64}, (size(zn)[1], 4, 4))
        gf_c_test[:, 1, 1] = gf_t[:, 1] ; gf_c_test[:, 1, 2] = gf_t[:, 3] ; gf_c_test[:, 1, 3] = gf_t[:, 4] ; gf_c_test[:, 1, 4] = gf_t[:, 3]
        gf_c_test[:, 2, 1] = gf_c_test[:, 1, 2]; gf_c_test[:, 2, 2] = gf_t[:, 2] ; gf_c_test[:, 2, 3] = gf_c_test[:, 1, 2] ; gf_c_test[:, 2, 4] = gf_t[:, 5]
        gf_c_test[:, 3, 1] = gf_c_test[:, 1, 3] ; gf_c_test[:, 3, 2] = gf_c_test[:, 2, 3] ; gf_c_test[:, 3, 3] = gf_c_test[:, 1, 1] ; gf_c_test[:, 3, 4] = gf_c_test[:, 1, 2]
        gf_c_test[:, 4, 1] = gf_c_test[:, 1, 4] ; gf_c_test[:, 4, 2] = gf_c_test[:, 2, 4] ; gf_c_test[:, 4, 3] = gf_c_test[:, 3, 4] ; gf_c_test[:, 4, 4] = gf_c_test[:, 2, 2]

        gf_c_test_02 = [(6.02348575e+00 -9.90968651e-01im),  (5.15813767e+00 -2.33078630e+00im)]#, (4.28316747e+00 -2.88633370e+00im), (6.68961495e-02 +1.31238862e-02im)]
        gf_c_test_00 = [(7.29148945e+00 -1.25688297e+00im), (6.37534240e+00 -3.08656815e+00im)]#,  (5.42233260e+00 -4.06248226e+00im), (5.55747767e+00 -2.25696764e+00im)]

        @test isapprox(gf_c, gf_c_test, atol=1e-6)
        @test isapprox(gf_c[:, 1, 3], gf_c_test_02, atol=1e-3)
        @test isapprox(gf_c[:, 1, 1], gf_c_test_00,  atol=1e-7)


   end


    @testset "test_t_to_c" begin

        gf_to_test = Green.c_to_to(zn, gf_c)
        gf_t_test = Green.c_to_t(zn, gf_c)

        gf_test1_c = Green.to_to_c(zn, gf_to_test)
        gf_test2_c = Green.t_to_c(zn, gf_t_test)

        @test isapprox(gf_c, gf_test1_c, rtol=1e-7, atol=1e-7)
        @test isapprox(gf_c, gf_test2_c, rtol=1e-7, atol=1e-7)

    end

 end  # testset

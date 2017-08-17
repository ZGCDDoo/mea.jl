using Base.Test
using Mea.Green
using Mea.Periodize
using JSON


@testset "Testing Periodize.jl" begin


#setup

    fin_gf_to = "./data/self_ctow0.dat"
    paramsfile = "./data/statsparams0.json"
    modelvec = Periodize.buildmodelvec(fin_gf_to, paramsfile)

    ss = div(size(modelvec.sEvec_c_)[1], 2)
    #println("ss = ", ss)
    (t, tp, mu) = (modelvec.t_, modelvec.tp_, modelvec.mu_)
    (wvec, sEvec_c) = (modelvec.wvec_[ss-1:ss], modelvec.sEvec_c_[ss-1:ss, :, :])
    modelvec_short = Periodize.ModelVector(t, tp, mu, wvec, sEvec_c)



    @testset "test_init" begin

         @test isapprox(modelvec.t_, 1.0)
         @test isapprox(modelvec.tp_, 0.40)
         @test isapprox(modelvec.mu_, 3.1736422868580827)
         @test isapprox(modelvec.wvec_[1:2], [-9.737999999999999545e+01, -5.046000000000000085e+01])

         model = Periodize.Model(modelvec, 1)
         cumulant = inv((model.w_ + model.mu_)*eye(Complex{Float64}, 4) - model.sE_)

         @test isapprox(model.t_, 1.0)
         @test isapprox(model.tp_, 0.4)
         @test isapprox(model.mu_, modelvec.mu_)
         @test isapprox(model.w_, -9.737999999999999545e+01)
         @test isapprox(cumulant, model.cumulant_)
         @test isapprox(cumulant, modelvec.cumulants_[1,:,:])

    end


    @testset "test_t_value" begin

        model = Periodize.Model(modelvec, 1)
        (kx , ky) = (rand(), -rand())
        t_value = Periodize.t_value(model, kx, ky)
        hop_test = Periodize.hopping_test(model, kx, ky)
        @test isapprox(t_value, hop_test)
    end


    @testset "test_eps0" begin
        model = Periodize.Model(modelvec, 1)
        (kx, ky) = (0.3, -0.19)
        eps0_value = Periodize.eps_0(model, kx, ky)
        @test isapprox(eps0_value, -4.66985, atol=1e-4)
    end


    @testset "test_exp_k" begin
        model = Periodize.Model(modelvec, 1)
        kx, ky = (0.3, -0.19)
        exp_k_value = Periodize.exp_k(kx, ky)
        real_value = [1.0, (0.955336 + 0.29552im),
                    (0.993956 + 0.109778im), (0.982004 - 0.188859im)]

        @test isapprox(exp_k_value, real_value, rtol=1e-4)
    end


    @testset "test_periodize_akw" begin
        t, tp = (modelvec.t_, modelvec.tp_)
        kx, ky = (0.0, pi/1.0)
        k = [kx, ky]
        ss = div(size(modelvec.sEvec_c_)[1], 2)
        mu = modelvec.mu_
        sE = modelvec.sEvec_c_[ss, :, :]
        ww = modelvec.wvec_[ss]
        model = Periodize.Model(modelvec, ss)
        N_c = 4
        r_sites = [[0.0, 0.0], [1.0,0.0], [1.0,1.0], [0.0,1.0]]
        gf_ktilde = inv((ww + mu)*eye(Complex{Float64}, 4) - Periodize.t_value(model, kx, ky) - sE)
        gf_w_lattice = 0.0 + 0.0im

        for i in 1:N_c
            for j in 1:N_c
                gf_w_lattice += 1.0/N_c * exp(-1.0im*dot(k, r_sites[i] - r_sites[j]) )*gf_ktilde[i, j]
            end
        end

       Akw = -2.0*imag(gf_w_lattice)
       Akw_test = Periodize.make_akwgreen(model)(k)

        Akw2 = (-2.0*imag(gf_w_lattice))^(2.0)
        Akw2_test = Periodize.make_akw2green(model)(k)

        @test isapprox(Akw, Akw_test)
        @test isapprox(Akw2, Akw2_test)
        @test isapprox(Akw, Akw_test)
        @test isapprox(Akw2, Akw2_test)

    end


    @testset "caldos" begin
        dos = Periodize.calcdos(modelvec_short)
        #println("dos = ", dos)
        dos_good = [-1.199999999999999956e-01 8.584885659790794099e-01;
                    -1.000000000000000056e-01 8.741665043156163772e-01
                    ]
        
        @test isapprox(dos, dos_good, atol=1e-4, rtol=1e-4)

    end


    @testset "caldos" begin
        dos2 = Periodize.calcdos2(modelvec_short)
        #println("dos2 = ", dos2)

        dos2_good =  [-0.12 2.20879; -0.1 2.29901]
        @test isapprox(dos2, dos2_good, atol=1e-4, rtol=1e-4)

    end



 end  # testset

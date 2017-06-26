using Base.Test
using Mea.Green
using Mea.Periodize
using JSON


@testset "Testing Periodize.jl" begin

    #setup
# #self_ctow_name = "self_ctow0_b12n0.495tp0.4U6.25.dat"

# (w_vec, sEvec_cw) = Green.readgreen_c(self_ctow_name, 1)

# t = 1.0; tp = 0.4;
# #mu = 3.173511234152317
# mu = JSON.parsefile("statsparams0.json")["mu"][1]
# modelvec = Periodize.ModelVector(t, tp, mu, w_vec, sEvec_cw)

# model = Periodize.Model(t, tp, mu, w_vec[1], sEvec_cw[1,:,:])
# tval = Periodize.t_value(model, 0.0, 0.0)
# tval2 = Periodize.hopping_test(model, 0.0, 0.0)
# @test isapprox(tval, tval2)
# #Periodize.
# Periodize.calcdos(modelvec, fout_name=fout_name)

#setup

    fin_gf_to = "./data/self_ctow.dat"
    paramsfile = "./data/statsparams0.json"
    modelvec = Periodize.buildmodelvec(fin_gf_to, paramsfile)


    
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

    # def test_periodize_Akw(self):
    #     """ """

        
    #     t, tp = (1.0, 0.4)
    #     kx, ky = (0.122, -0.987)
    #     k = np.array([kx,ky])
    #     sE = np.random.rand(4, 4)
    #     sEarr = np.array([sE], dtype=complex)
    #     ww = random()
    #     mu = random()
    #     model = periodize.Model(t, tp, mu, np.array([ww]), sEarr)
    #     N_c = 4
    #     r_sites = np.array([[0.0, 0.0], [1.0,0.0], [1.0,1.0], [0.0,1.0]])
    #     gf_ktilde = linalg.inv((ww + mu)*np.eye(4) - model.t_value(kx, ky) - sE)

    #     gf_w_lattice = 0.0 + 0.0j

    #     for i in range(N_c):
    #         for j in range(N_c):
    #             gf_w_lattice += 1/N_c * np.exp(-1.0j*np.dot(k, r_sites[i] - r_sites[j]) )*gf_ktilde[i, j]
        
    #     Akw = -2.0*gf_w_lattice.imag
    #     Akw_test = model.periodize_Akw(kx, ky, 0)
    #     Akw2 = (-2.0*gf_w_lattice.imag)**(2.0)
    #     Akw2_test = (model.periodize_Akw2(kx, ky ,0))
        
    #     try:
    #         test_tools.compare_arrays(Akw, Akw_test)
    #         test_tools.compare_arrays(Akw2, Akw2_test)
    #         np.testing.assert_allclose(Akw, Akw_test)
    #         np.testing.assert_allclose(Akw2, Akw2_test)
    #     except AssertionError:
    #         self.fail("np all close failed at test_periodize_Akw")         



    # def test_periodize_Gkz(self):
    #     """ """
    #     t, tp = (1.0, 0.4)
    #     kx, ky = (0.122, -0.987)
    #     sE = np.random.rand(4, 4)
    #     sEarr = np.array([sE], dtype=complex)
    #     ww = random()
    #     mu = random()
    #     model = periodize.Model(t, tp, mu, np.array([ww]), sEarr)
    #     Akw = model.periodize_Akw(kx, ky, 0)
    #     Akw_test = -2.0*model.periodize_Gkz_vec(kx, ky).imag

    #     try:
    #         self.assertAlmostEqual(Akw, Akw_test[0])
    #     except AssertionError:
    #         self.fail("np all close failed at test_periodize_Gkw")   



 end  # testset

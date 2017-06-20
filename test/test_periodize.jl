# push!(LOAD_PATH, pwd())
# #println( pwd() )
# #println(LOAD_PATH)
# using Base.Test

# include("green.jl")
# include("periodize.jl")
# using Green
# using Periodize
# import JSON

# self_ctow_name = ARGS[1]
# fout_name = ARGS[2]
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


using Base.Test

@testset "Testing Periodize.jl" begin
 println("Passing for now.")
#     # set up
#     phi = [.95, -.4, -.4]
#     theta = zeros(3)
#     sigma = .15
#     lp = ARMA(phi, theta, sigma)

#     # test simulate
#     sim = simulation(lp, ts_length=250)
#     @test length(sim) == 250

#     # test impulse response
#     imp_resp = impulse_response(lp, impulse_length=75)
#     @test length(imp_resp) == 75

#     @testset "test constructors" begin
#         phi = 0.5
#         theta = 0.4
#         sigma = 0.15

#         a1 = ARMA(phi, theta, sigma)
#         a2 = ARMA([phi;], theta, sigma)
#         a3 = ARMA(phi, [theta;], sigma)

#         for nm in fieldnames(a1)
#             @test getfield(a1, nm) == getfield(a2, nm)
#             @test getfield(a1, nm) == getfield(a3, nm)
#         end
#     end

#     @testset "test autocovariance" begin
#         θ = 0.5
#         σ = 0.15
#         ma1 = ARMA(Float64[], [θ], σ)
#         ac = autocovariance(ma1; num_autocov=5)

#         # first is the variance. equal to (1 + θ^2) sigma^2
#         @test isapprox(ac[1], (1+θ^2)*σ^2; atol=1e-3)

#         # second should be θ σ^2
#         @test isapprox(ac[2], θ*σ^2; atol=1e-3)

#         # all others should be 0
#         @test isapprox(ac[3:end], zeros(ac[3:end]); atol=1e-3)
#     end

 end  # testset

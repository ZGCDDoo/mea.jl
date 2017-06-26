

module Periodize

using Cubature: hcubature
using JSON
using Mea.Green


II = [1.0 + 0.0im 0.0 0.0 0.0
      0.0 1.0 0.0 0.0
      0.0 0.0 1.0 0.0
      0.0 0.0 0.0 1.0]


function buildmodelvec(finsE::String, finparams::String)
    params = JSON.parsefile(finparams)
    t = 1.0
    tp = params["tp"][1]
    mu = params["mu"][1]
    (wvec, sEvec_c) = Green.readgreen_c(finsE)    
    return ModelVector(t, tp, mu, wvec, sEvec_c)
end




type ModelVector
   t_::Float64 ; tp_::Float64 ; mu_::Float64
   wvec_::Array{Float64, 1} ; sEvec_c_::Array{Complex{Float64}, 3}
   cumulants_::Array{Complex{Float64}, 3}

   function ModelVector(t::Float64, tp::Float64, mu::Float64,
                        wvec::Array{Float64, 1}, sEvec_c::Array{Complex{Float64}, 3})

        cumulants = build_cumulants(wvec, mu, sEvec_c)
        return new(t, tp, mu, wvec, sEvec_c, cumulants)
    end


end


function build_cumulants(wvec::Array{Float64, 1}, mu::Float64, sEvec_c::Array{Complex{Float64}, 3})

      cumulants = zeros(Complex{Float64}, size(sEvec_c))

      for (ii, ww) in enumerate(wvec)
          tmp = zeros(Complex{Float64}, (4, 4))
          tmp[1, 1] = tmp[2, 2] = tmp[3, 3] = tmp[4, 4] = (ww + mu)
          tmp -= sEvec_c[ii, :, :]
          cumulants[ii, :, :] = inv(tmp)
      end

      return cumulants
end


type Model
    t_::Float64  ; tp_::Float64 ; mu_::Float64
    w_::Float64 ; sE_::Array{Complex{Float64}, 2}
    cumulant_::Array{Complex{Float64}, 2}
    
    function Model(modelvec::ModelVector, ii::Integer)
        (t, tp, mu, w, sE_c, cumulant) = (modelvec.t_, modelvec.tp_, modelvec.mu_, modelvec.wvec_[ii], modelvec.sEvec_c_[ii, :, :], modelvec.cumulants_[ii, :, :])
        return new(t, tp, mu, w, sE_c, cumulant)
    end

end


function t_value(model::Model, kx::Float64, ky::Float64) # This is t_ij(k_tilde)
    t = model.t_; tp = model.tp_
    t_val = zeros(Complex{Float64}, (4, 4))
    ex = exp(-2.0im*kx) ; emx = conj(ex)
    ey = exp(-2.0im*ky) ; emy = conj(ey)
    tloc = [0.0 -t -tp -t
            -t 0.0 -t 0.0
            -tp -t 0.0 -t
            -t 0.0 -t 0.0]

    t_val += tloc

    t_val[1, 1] += 0.0;               t_val[1, 2] += -t*ex;               t_val[1, 3] += -tp*ex*ey;    t_val[1, 4] += -t*ey
    t_val[2, 1] += -t*emx;            t_val[2, 2] += 0.0;                 t_val[2, 3] += -t*ey;        t_val[2, 4] += -tp*(emx + ey)
    t_val[3, 1] += -tp*emx*emy;       t_val[3, 2] += -t*emy;              t_val[3, 3] += 0.0;          t_val[3, 4] += -t*emx
    t_val[4, 1] += -t*emy;            t_val[4, 2] += -tp*(ex + emy);      t_val[4, 3] += -t*ex;        t_val[4, 4] += 0.0

    return (t_val)
end


function exp_k(kx::Float64, ky::Float64)

    expk = zeros(Complex{Float64}, 4) # Here exp_k[i] is in fact e**(-j*dot(k, r_i)) where r_i is a site of the cluster
    expk[1] = 1.0
    expk[2] = exp(1.0im*kx)
    expk[3] = exp(1.0im*(kx+ky))
    expk[4] = exp(1.0im*ky)
    return expk
end


function build_gf_ktilde(model::Model, kx::Float64, ky::Float64)
    return inv((model.w_ + model.mu_) * II - t_value(model, kx, ky) - model.sE_)
end


function build_gf_ktilde_inverse(model::Model, kx::Float64, ky::Float64)
    return( (model.w_ + model.mu_) * II - t_value(model, kx, ky) - model.sE_)
end


function make_periodize(model::Model)
  function periodize_Akw(kk::Array{Float64, 1}) # periodize the imaginary part (Ak_w)

      (kx, ky) = kk
      N_c = 4.0
      gf_ktilde = build_gf_ktilde(model, kx, ky)
      expk = exp_k(kx, ky)
      return imag(-2.0/N_c * dot(expk, (gf_ktilde * expk)) )
  end
  return periodize_Akw
end


function make_akw2(model::Model)
    function akw2(kk)
        akw = make_periodize(model)
        return (akw(kk)^2.0)
    end
    return akw2
end


function eps_0(model::Model, kx::Float64, ky::Float64)
    return (-2.0*model.t_*(cos(kx) + cos(ky)) - 2.0*model.tp_*cos(kx + ky) )
end


function hopping_test(model::Model, kx::Float64, ky::Float64)
    t = model.t_ ; tp = model.tp_
    k = [kx, ky]
    N_c = 4
    r_sites = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    K_sites = pi*deepcopy(r_sites)
    t_array = zeros(Complex{Float64}, (N_c, N_c))

    for i in 1:N_c
        for j in 1:N_c
            for K in K_sites
                t_array[i, j] += 1.0/N_c * exp(1.0im*dot(K + k, r_sites[i] - r_sites[j])) * eps_0(model, (K + k)...)
            end
        end
    end

    return t_array

end


function calcdos(modelvector::ModelVector; fout_name::String="dos.dat", maxevals::Int64=0)
    
    len_sEvec_c = size(modelvector.sEvec_c_)[1]
    dos = zeros(Float64, len_sEvec_c)

    for n in 1:len_sEvec_c
        #println("IN LOOP of dos # ", n, " out of ", len_sEvec_c)
        model = Model(modelvector, n)
        dos[n] = (2.0*pi)^(-2.0)*hcubature(make_periodize(model), (-pi, -pi), (pi, pi), reltol=1e-8, abstol=1e-8, maxevals=maxevals)[1]
    end
    dos_out = hcat(modelvector.wvec_, dos)
    writedlm(fout_name, dos_out, " ")
    return dos

end


function calcdos2(modelvector::ModelVector; fout_name::String="dos2.dat", maxevals::Int64=0)

    len_sEvec_c = size(modelvector.sEvec_c_)[1]
    dos2 = zeros(Float64, len_sEvec_c)

    for n in 1:len_sEvec_c
        #println("IN LOOP of dos # ", n, " out of ", len_sEvec_c)
        model = Model(modelvector, n)
        dos2[n] = (2.0*pi)^(-2.0)*hcubature(make_akw2(model), (-pi, -pi), (pi, pi), reltol=1.49e-8, abstol=1.49e-8, maxevals=maxevals)[1]
    end
    dos2_out = hcat(modelvector.wvec_, dos2)
    writedlm(fout_name, dos2_out, " ")
    return dos2_out
end


end

module Transport
using JSON


using Mea.Periodize
const cutoffdefault = 20.0

function dfdw(beta::Float64, ww::Float64)
    return (-beta*exp(beta*ww)/(1.0 + exp(beta*ww))^2.0   )
end


function vz2_int(kk::Array{Float64 ,1}) #v_perp^2.0 integrated in z

    (kx, ky) = (kk[1], kk[2])
    result = 2.0
    #result = 2.0*(cos(kx) - cos(ky))^(4.0)
    return result
end


 function make_lab_kintegrand(model::Periodize.Model; fctper::String="make_akw2green")
    #fctper = "make_akw2green, make_akw2cum, make_akw2trace" as defined in Periodize.jl  
    make_akw2 = getfield(Periodize, Symbol(fctper))
    akw2 = make_akw2(model)

    function lab_kintegrand(kk::Array{Float64, 1})
          return(akw2(kk)*2.0)#vz2_int(kk))
      end
       return lab_kintegrand
  end

  function make_lab_kintegrand_cuba(model::Periodize.Model; fctper::String="make_akw2green")
    
    make_akw2 = getfield(Periodize, Symbol(fctper))
    akw2 = make_akw2(model)
    
    function lab_kintegrand_cuba(kk::Array{Float64, 1}, ff::Array{Float64, 1})
        kk_scaled = (-pi + 2*pi*kk)
        ff[1] = (4.0*pi*pi*akw2(kk_scaled)*2.0)#vz2_int(kk_scaled))
    end
    return lab_kintegrand_cuba
   end

function calc_labk(modelvec::Periodize.ModelVector, beta::Float64; cutoff::Float64=cutoffdefault, maxevals::Int64=100000,
                    libintegrator::String="cubature", fctper::String="make_akw2green") #calculate the k integral of the L_ab coeffiecient
    #returns: Array{Float64, 1} size of w_vec that will be integrated in frequency

  modelvec = deepcopy(modelvec)
  cutoff = cutoff
  cutoffidx=0
  for (ii, ww) in enumerate(modelvec.wvec_)
    if abs(ww*beta) < cutoff
      cutoffidx = ii
      break
    end
  end

  modelvec.wvec_ = modelvec.wvec_[cutoffidx:end-cutoffidx]
  modelvec.sEvec_c_ = modelvec.sEvec_c_[cutoffidx:end-cutoffidx, :, :]

  kwargs = Dict(:fctper=>fctper)
  #println("kwargs... = ", kwargs...)
  if libintegrator == "cubature"
    result = Periodize.calcintegral(modelvec, make_lab_kintegrand, maxevals=maxevals; kwargs...)
  elseif libintegrator == "cuba"
    result = Periodize.calcintegral_cuba(modelvec, make_lab_kintegrand_cuba, maxevals=maxevals; kwargs...)
  else
    error("Invalid value for libintegrator")
  end

  integrand_w = zeros(Float64, size(result)[1])
  for ii in 1:size(result)[1]
      integrand_w[ii] = -dfdw(beta, result[ii, 1])*result[ii, 2]
  end


  return (modelvec.wvec_, integrand_w)

end


function calc_l11(integrand_w::Array{Float64, 1}, wvec::Array{Float64, 1}, beta::Float64)

    l11 = 0.0
    assert(size(wvec) == size(integrand_w))

   # frequency integration
   for ii in 1:size(wvec)[1]-1                       # y                                 x
       l11 += 0.5 * (integrand_w[ii] + integrand_w[ii+1]) * (wvec[ii+1] - wvec[ii]) / (2.0*pi)
   end

   l11 /= beta
   #println(l11)
   return l11
end


function calc_l21(integrand_w::Array{Float64, 1}, wvec::Array{Float64, 1}, beta::Float64)

    l21 = 0.0
    assert(size(wvec) == size(integrand_w))

    for jj in 1:size(integrand_w)[1]
        integrand_w[jj] *= wvec[jj]
    end

   # frequency integration
   for ii in 1:size(wvec)[1]-1                       # y                                 x
       l21 += 0.5 * (integrand_w[ii] + integrand_w[ii+1]) * (wvec[ii+1] - wvec[ii]) / (2.0*pi)
   end

   l21 /= beta
   #println(l21)
    return l21
end


function calc_l22(integrand_w::Array{Float64, 1}, wvec::Array{Float64, 1}, beta::Float64)

    l22 = 0.0
    assert(size(wvec) == size(integrand_w))

    for jj in 1:size(integrand_w)[1]
        integrand_w[jj] *= (wvec[jj])^(2.0)
    end

   # frequency integration
   for ii in 1:size(wvec)[1]-1                       # y                                 x
       l22 += 0.5 * (integrand_w[ii] + integrand_w[ii+1]) * (wvec[ii+1] - wvec[ii]) / (2.0*pi)
   end

   l22 /= beta
   #println(l22)
   return l22
end



function calc_sigmadc(modelvec::Periodize.ModelVector, beta::Float64; cutoff::Float64=cutoffdefault)

    (wvec, integrand_w) = calc_labk(modelvec, beta, cutoff=cutoff)
    sigmadc = beta*calc_l11(integrand_w, wvec,  beta)
    #println("sigmadc = ", sigmadc)
    return sigmadc
end



# function calc_seebeck(modelvec::Periodize.ModelVector, beta::Float64; cutoff=cutoffdefault)
#   seebeck = -beta*calc_l21(modelvec, beta, cutoff=cutoff)/calc_l11(modelvec, beta, cutoff=cutoff)
#   #println("seebeck = ", seebeck)
#   return seebeck    
# end


function calc_n(modelvec::Periodize.ModelVector, beta::Float64; fout_name::String="dos.dat", maxevals::Int64=100000)
    dos = Periodize.calcdos(modelvec, fout_name=fout_name, maxevals=maxevals)
    n = 0.0
    for jj in 1:size(dos)[1]
        ww = dos[jj, 1]
        dos[jj, 2] *= 1.0/(exp(beta*ww) + 1.0)
    end

   # frequency integration
   for ii in 1:size(dos)[1]-1     # y                                 x
       n += 0.5 * (dos[ii, 2] + dos[ii+1, 2]) * (dos[ii+1, 1] - dos[ii, 1]) / (2.0*pi)
   end

   #println("n = ", n)
   return n
end


function coefstrans(modelvec::Periodize.ModelVector, beta::Float64; cutoff::Float64=cutoffdefault, fout_name::String="dos.dat", maxevals::Int64=100000,
                    libintegrator::String="cubature", fctper::String="make_akw2green")
    
    dd = Dict{String, Float64}()
    (wvec, integrand_w) = calc_labk(modelvec, beta, cutoff=cutoff, maxevals=maxevals, libintegrator=libintegrator, fctper=fctper)
    dd["l11"] =  calc_l11(integrand_w, wvec, beta)
    dd["l21"] = calc_l21(integrand_w, wvec, beta)
    dd["l22"] = calc_l22(integrand_w, wvec, beta)
    dd["sigmadc"] = beta*dd["l11"]
    dd["seebeck"] = -beta*dd["l21"]/dd["l11"]
    dd["n"] = calc_n(modelvec, beta, fout_name=fout_name, maxevals=maxevals)

    fout = "coefstrans" * fctper * ".json"
    if isfile(fout)
        mv(fout, fout * "." *  Dates.format(now(), "YYYY-uu-HHh:MM") )
    end
    
    open(fout, "w") do f
            write(f, JSON.json(dd))
    end
    
    return dd
end




end
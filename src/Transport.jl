module Transport

using Mea.Periodize


function dfdw(beta::Float64, ww::Float64)
    return (-beta*exp(beta*ww)/(1.0 + exp(beta*ww))^2.0   )
end


function vz2_int(kk::Array{Float64 ,1}) #v_perp^2.0 integrated in z

    (kx, ky) = (kk[1], kk[2])
    result = 2.0
    #result = 2.0*(cos(kx) - cos(ky))^(4.0)
    return result
end



function calc_labk(modelvec::Periodize.ModelVector, beta::Float64) #calculate the k integral of the L_ab coeffiecient
    #returns: Array{Float64, 1} size of w_vec that will be integrated in frequency

  cutoff = 20.0
  cutoffidx=0
  for (ii, ww) in enumerate(modelvec.wvec_)
    if abs(ww*beta) < cutoff
      cutoffidx = ii
      break
    end
  end

  modelvec.wvec_ = modelvec.wvec_[cutoffidx:end-cutoffidx]
  modelvec.sEvec_c_ = modelvec.sEvec_c_[cutoffidx:end-cutoffidx, :, :]

  function make_lab_kintegrand(model::Periodize.Model)
      function lab_kintegrand(kk::Array{Float64, 1})
          akw2 = Periodize.make_akw2(model)
          return(akw2(kk)*vz2_int(kk))
      end
  end

  result = Periodize.calcintegral(modelvec, make_lab_kintegrand, maxevals=100000)

  integrand = zeros(Float64, size(result)[1])
  for ii in 1:size(result)[1]
      integrand[ii] = -dfdw(beta, result[ii, 1])*result[ii, 2]
  end

  return (integrand)

end


function calc_l11(modelvec::Periodize.ModelVector, beta::Float64)

    integrand = calc_labk(modelvec, beta)
    l11 = 0.0

    wvec = modelvec.wvec_
    assert(size(wvec) == size(integrand))

   # frequency integration
   for ii in 1:size(wvec)[1]-1                       # y                                 x
       l11 += 0.5 * (integrand[ii] + integrand[ii+1]) * (wvec[ii+1, 1] - wvec[ii, 1]) / (2.0*pi)
   end

   l11 /= beta
   #println(l11)
   return l11
end


function calc_l21(modelvec::Periodize.ModelVector, beta::Float64)

    integrand = calc_labk(modelvec, beta)
    l21 = 0.0
    wvec = modelvec.wvec_
    assert(size(wvec) == size(integrand))

    for jj in 1:size(integrand)[1]
        integrand[jj] *= wvec[jj]
    end

   # frequency integration
   for ii in 1:size(wvec)[1]-1                       # y                                 x
       l21 += 0.5 * (integrand[ii] + integrand[ii+1]) * (wvec[ii+1, 1] - wvec[ii, 1]) / (2.0*pi)
   end

   l21 /= beta
   #println(l21)
    return l21
end


function calc_l22(modelvec::Periodize.ModelVector, beta::Float64)

    integrand = calc_labk(modelvec, beta)
    l22 = 0.0
    wvec = modelvec.wvec_
    assert(size(wvec) == size(integrand))

    for jj in 1:size(integrand)[1]
        integrand[jj] *= (wvec[jj])^(2.0)
    end

   # frequency integration
   for ii in 1:size(wvec)[1]-1                       # y                                 x
       l21 += 0.5 * (integrand[ii] + integrand[ii+1]) * (wvec[ii+1, 1] - wvec[ii, 1]) / (2.0*pi)
   end

   l22 /= beta
   #println(l22)
   return l22
end


function calc_sigmadc(modelvec::Periodize.ModelVector, beta::Float64)

    sigmadc = beta*calc_l11(modelvec, beta)
    println(sigmadc)
    return sigmadc
end

# function calc_seebeck(...)
# end

end

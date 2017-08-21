module SigmaDC


using Mea.Periodize


function dfdw(beta_::Float64, ww::Float64)

    #if abs(beta_*ww) > 20.0
    #    return 0.0
    #end

    return (-beta_*exp(beta_*ww)/(1.0 + exp(beta_*ww))^2.0   )
end


function calc_sigmadc(modelvec::Periodize.ModelVector, beta_::Float64)

  cutoff = 20.0
  cutoffidx=0
  for (ii, ww) in enumerate(modelvec.wvec_)
    if abs(ww*beta_) < cutoff
      cutoffidx = ii
      #println("ii = ", ii)
      break
    end
  end

  #println("cutoffidx = ", cutoffidx)

  modelvec.wvec_ = modelvec.wvec_[cutoffidx:end-cutoffidx]
  modelvec.sEvec_c_ = modelvec.sEvec_c_[cutoffidx:end-cutoffidx,:,:]

  dos2 = Periodize.calcdos2(modelvec, maxevals=100000)

  integrand = zeros(Float64, size(dos2)[1])
  for ii in 1:size(dos2)[1]
      integrand[ii] = dfdw(beta_, dos2[ii, 1])*dos2[ii, 2]
  end

  sigmadc = 0.0

      for ii in 1:size(dos2)[1]-1
                              # y                                 x
          sigmadc += 0.5 * (integrand[ii] + integrand[ii+1]) * (dos2[ii+1, 1] - dos2[ii, 1])
      end
  sigmadc*= -2.0/(2.0*pi)
  println(sigmadc)
  return sigmadc
end


end

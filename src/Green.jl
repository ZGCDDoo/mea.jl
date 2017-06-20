#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Here every gf is a gfvec, i.e every green-function is a vector of green function.
To lighten the writing we supress the "vec". Maybe consider putting everything as in the README.txt.
"""

module Green

using Base.Test

function readgreen_c(fin_gf_to::String, zn_col::Int = 1)

    gf_to = readdlm(fin_gf_to)

    # Get the Matsubara frequencies on imaginary axis
    zn_vec = deepcopy(gf_to[:, zn_col])

    # construct complex gf_t (each gf is formed by the consecutif real part
    # (pair columns) and imaginary part (odd columns))
    gf_c = to_to_c(zn_vec, gf_to)
    return (zn_vec, gf_c)

end


function to_to_t(gf_to::Array{Float64, 2})
    assert(size(gf_to)[2] == 11)
    gf_to = deepcopy(gf_to)
    return (1.0*gf_to[:, 2:2:end] + 1.0im*gf_to[:, 3:2:end])
end


function t_to_to(zn_vec::Array{Float64, 1}, gf_t::Array{Complex{Float64}, 2})

    zn_vec = deepcopy(zn_vec) ; gf_t = deepcopy(gf_t)
    assert(size(zn_vec)[1] == size(gf_t)[1])
    assert(size(gf_t)[2] == 5)

    gf_to = zeros(Float64, (size(zn_vec)[1], 11))

    gf_to[:, 1] = zn_vec
    gf_to[:, 2:2:end] =  real(gf_t)
    gf_to[:, 3:2:end] =  imag(gf_t)

    return gf_to
end


function t_to_c(zn_vec::Array{Float64, 1}, gf_t::Array{Complex{Float64}, 2})
    """convert from tabular form to cluster form """
    zn_vec = deepcopy(zn_vec) ; gf_t = deepcopy(gf_t)
    #assert(size(zn_vec)[1] == size(gf_t)[1])
    #assert(size(gf_t) == 5)

    gf_c = zeros(Complex{Float64}, (size(zn_vec)[1], 4, 4))

    gf_c[:, 1, 1] = gf_t[:, 1] ; gf_c[:, 1, 2] = gf_t[:, 3] ; gf_c[:, 1, 3] = gf_t[:, 4] ; gf_c[:, 1, 4] = gf_t[:, 3]
    gf_c[:, 2, 1] = gf_t[:, 3] ; gf_c[:, 2, 2] = gf_t[:, 2] ; gf_c[:, 2, 3] = gf_t[:, 3] ; gf_c[:, 2, 4] = gf_t[:, 5]
    gf_c[:, 3, 1] = gf_t[:, 4] ; gf_c[:, 3, 2] = gf_t[:, 3] ; gf_c[:, 3, 3] = gf_t[:, 1] ; gf_c[:, 3, 4] = gf_t[:, 3]
    gf_c[:, 4, 1] = gf_t[:, 3] ; gf_c[:, 4, 2] = gf_t[:, 5] ; gf_c[:, 4, 3] = gf_t[:, 3] ; gf_c[:, 4, 4] = gf_t[:, 2]

    return gf_c
end


function c_to_t(zn_vec::Array{Float64, 1}, gf_c::Array{Complex{Float64}, 3})
    """ """
    zn_vec = deepcopy(zn_vec) ; gf_c = deepcopy(gf_c)
    gf_to = c_to_to(zn_vec, gf_c)
    gf_t = to_to_t(gf_to)
    return gf_t
end


function to_to_c(zn_vec::Array{Float64, 1}, gf_to::Array{Float64, 2})
    """convert from tabular-out form to cluster form"""
    zn_vec = deepcopy(zn_vec) ; gf_to = deepcopy(gf_to)

    # println("size(zn_vec) =  ", size(zn_vec))
    # println("size(gf_to) = ", size(gf_to))

    assert(size(zn_vec)[1] == size(gf_to)[1])
    assert(size(gf_to)[2] == 11)

    # the following three lines could be replaced by: gf_t = to_to_t(gf_to)
    # but I do this for tests purposes (np.testing)
    gf_t = zeros(Complex{Float64}, (size(gf_to)[1], 5))
    gf_t = 1.0*gf_to[:, 2:2:end] + 1.0im*gf_to[:, 3:2:end]
    @test isapprox(gf_t, to_to_t(gf_to))

    gf_c = zeros(Complex{Float64}, (size(zn_vec)[1], 4, 4))

    gf_c[:, 1, 1] = gf_t[:, 1] ; gf_c[:, 1, 2] = gf_t[:, 3] ; gf_c[:, 1, 3] = gf_t[:, 4] ; gf_c[:, 1, 4] = gf_t[:, 3]
    gf_c[:, 2, 1] = gf_t[:, 3] ; gf_c[:, 2, 2] = gf_t[:, 2] ; gf_c[:, 2, 3] = gf_t[:, 3] ; gf_c[:, 2, 4] = gf_t[:, 5]
    gf_c[:, 3, 1] = gf_t[:, 4] ; gf_c[:, 3, 2] = gf_t[:, 3] ; gf_c[:, 3, 3] = gf_t[:, 1] ; gf_c[:, 3, 4] = gf_t[:, 3]
    gf_c[:, 4, 1] = gf_t[:, 3] ; gf_c[:, 4, 2] = gf_t[:, 5] ; gf_c[:, 4, 3] = gf_t[:, 3] ; gf_c[:, 4, 4] = gf_t[:, 2]

    @test isapprox(gf_c, t_to_c(zn_vec, gf_t))
    #println("gf_c = ")
    #println(gf_c)
    #println("gf_c = ")
    #println(t_to_c(zn_vec, gf_t))
    return gf_c
end


function c_to_to(zn_vec::Array{Float64}, gf_c::Array{Complex{Float64}, 3})
    """convert from cluster form to tabular-out form """
    zn_vec = deepcopy(zn_vec) ; gf_c = deepcopy(gf_c)
    assert(size(zn_vec)[1] == size(gf_c)[1])
    assert(size(gf_c)[2] == size(gf_c)[2]) ; assert(size(gf_c)[2] ==  size(gf_c)[3], 4 )

    gf_to = zeros(Float64, (size(zn_vec)[1], 11))


    gf_to[:, 1] = zn_vec
    gf_to[:, 2] = gf_c[:, 1, 1].real ; gf_to[:, 3] = gf_c[:, 1, 1].imag
    gf_to[:, 4] = gf_c[:, 2, 2].real ; gf_to[:, 5] = gf_c[:, 2, 2].imag
    gf_to[:, 6] = gf_c[:, 1, 2].real ; gf_to[:, 7] = gf_c[:, 1, 2].imag
    gf_to[:, 8] = gf_c[:, 1, 3].real ; gf_to[:, 9] = gf_c[:, 1, 3].imag
    gf_to[:, 10] = gf_c[:, 2, 4].real ; gf_to[:, 11] = gf_c[:, 2, 4].imag

    return gf_to
end

end

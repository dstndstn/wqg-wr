using GraphPlot
using LightGraphs
using SparseArrays
using Plots

## can be replaced by LightGraphs.SimpleGraphs.grid
export gridGraph
function gridGraph(lgs,m)
    sk = path_graph(lgs[length(lgs)])
    for i in length(lgs)-1:-1:1
       sk = cartesian_product(sk , path_graph(lgs[i]))
    end
    return sk
end

struct Edgevv
    v1::Int
    v2::Int
end

## impose condition on Edgevv so that vertices are ordered?

export evv_order
function evv_order(evv::Edgevv) # order vertices
    if evv.v1 > evv.v2
        return Edgevv(evv.v2,evv.v1)
    else
        return evv
    end
end

## line graph dsk of sk

"""refer to dual skeleton vertex"""
struct Dskv
    v::Int
end

# [setup]

export l_evv_vertices
# given a list of vertices, return all relevant evv
function l_evv_vertices(sk,lv)
    list = Edgevv[]
    for i in lv, j in neighbors(sk,i)
        push!(list, evv_order(Edgevv(i,j)))
    end
    return unique(list)
end

using SparseArrays

export value_hom
function value_hom(sk, D, vl)
    n = nv(sk)
    value = spzeros(n,n)
    for i in edges(sk)
        value[i.src,i.dst] = vl
    end
    return value
end

export set_value
# update val with non-zero elements of var
function set_value(value::SparseMatrixCSC, var::SparseMatrixCSC)
    val = copy(value)
    for i in findall(!iszero, var)
        v1, v2 = i.I[1], i.I[2]
        val[v1,v2] = var[v1,v2]
    end
    return val
end

function set_value(sk, levv, var_vec)
# set configuration var_vec on edges levv, output a sparse matrix
    l = length(var_vec)
    if l != length(levv)
        return 0.0
    end
    n = nv(sk)
    value = spzeros(n,n)
    for j in 1:l
        i = levv[j]
        value[i.v1,i.v2] = var_vec[j]
    end
    value
end

# [corners]

using IterTools

export corners
function corners(sk, D, evv::Edgevv)
    v1,v2 = evv.v1, evv.v2
    list = []
    for v in [v1,v2]
        vo = filter!(x->x≠v, [v1,v2])[1]
        nb = collect(neighbors(sk,v))
        nb = filter!(x->x≠vo, nb)
        for i in subsets(nb, D-1)
            li = []
            for j in i
                push!(li, evv_order(Edgevv(v,j))) # order vertices
            end
            push!(list, li)
        end
    end
    return list
end

export corners_α_w
function corners_α_w(sk, D, evv::Edgevv, a, wvalue::SparseMatrixCSC)
    n = length(corners(sk, D, evv))
    v = 0
    for i in corners(sk, D, evv)
        rs = 1
        for j in i
            v1,v2 = j.v1, j.v2
            rs *= 1/wvalue[v1,v2]
        end
    v += rs
    end
    v *= a/n
    return v
end

# [amplitudes]

export vvmd_r
# d is spatial dimension, which equals D-1 where D is the spacetime dimension
function vvmd_r(x,d)
    sinc(x/pi)^(-d)
end


export campEdge_wr
function campEdge_wr(sk, D, wvalue, rvalue, evv::Edgevv, α)
    evvo = evv_order(evv)
    v1, v2 = evvo.v1, evvo.v2
    w = wvalue[v1,v2]
    r = rvalue[v1,v2]
    vd = real(vvmd_r((r+0im)^(1/2),d))
    cf = vd^(-1)*log(vd)
    return 1/(0.5+72*(α*cf*w)^2)
end

export campEdge_wr_measure
function campEdge_wr_measure(sk, D, wvalue, rvalue, ws, rs, evv::Edgevv, α, ϵ)
    # ϵ as IR regulator
    ## is the regulartor necessary?
    # ws and rs rescales w and r values
    evvo = evv_order(evv)
    v1, v2 = evvo.v1, evvo.v2
    w = wvalue[v1,v2]/ws
    r = rvalue[v1,v2]/rs + 0im
    vd = real(vvmd_r((r+0im)^(1/2),d))
    cf = vd^(-1)*log(vd)
    return real(1/((0.5+72*(α*cf*w)^2)*w+ϵ))
end

export campEdge_b
# boundary modulus (same for spacelike and timelike)
function campEdge_b(sk, D, wvalue, rvalue, evv::Edgevv, α)
    evvo = evv_order(evv)
    v1, v2 = evvo.v1, evvo.v2
    w = wvalue[v1,v2]
    r = rvalue[v1,v2]
    vd = real(vvmd_r((r+0im)^(1/2),d))
    cf = vd^(-1)*log(vd)
    return 1/sqrt(1+(12*α*cf*w)^2)
end

# this does not distinguish boundary and interior edges
export camp_VG_wr
function camp_VG_wr(sk, D, wvalue, rvalue, a)    
    amp = 1.0
    for i in edges(sk)
        evv = Edgevv(i.src, i.dst)
        α = corners_α_w(sk, D, evv, a, wvalue)
        amp *= campEdge_wr(sk, D, wvalue, rvalue, evv, α)
    end
    return amp
end

function camp_VG_wr(sk, D, wvalue, rvalue, a, levv_i, levv_b)
    # levv is the list of interior edges
    amp = 1.0
    for i in levv_i
        α = corners_α_w(sk, D, i, a, wvalue)
        amp *= campEdge_wr(sk, D, wvalue, rvalue, i, α)
    end
    # levv_b is the list of boundary edges
    for j in levv_b
        α = corners_α_w(sk, D, j, a, wvalue)
        amp *= campEdge_b(sk, D, wvalue, rvalue, j, α)
    end
    return amp
end

# this does not distinguish boundary and interior edges
export camp_VG_wr_measure
function camp_VG_wr_measure(sk, D, wvalue, rvalue, ws, rs, a, ϵ)    
    amp = 1.0
    for i in edges(sk)
        evv = Edgevv(i.src, i.dst)
        α = corners_α_w(sk, D, evv, a, wvalue)
        amp *= campEdge_wr_measure(sk, D, wvalue, rvalue, ws, rs, evv, α, ϵ)
    end
    return amp
end

# left out factor d from Jacobian (s,rho)->(w,r)
function camp_VG_wr_measure(sk, D, wvalue, rvalue, ws, rs, a, ϵ, levv_i, levv_b)    
    amp = 1.0
    for i in levv_i
        α = corners_α_w(sk, D, i, a, wvalue)
        amp *= campEdge_wr_measure(sk, D, wvalue, rvalue, ws, rs, i, α, ϵ)
    end
    for j in levv_b
        α = corners_α_w(sk, D, j, a, wvalue)
        amp *= campEdge_b(sk, D, wvalue, rvalue, j, α)
    end
    return amp
end

# [volume]

export volume_edge
function volume_edge(sk, D, evv::Edgevv, b, wvalue, rvalue)
    # b is a parameter because the volume is defined up to this constant
    d = D-1
    v1, v2 = evv.v1, evv.v2 # assuming v1 < v2
    w = wvalue[v1, v2]
    r = rvalue[v1, v2]
    vd = real(vvmd_r((r+0im)^(1/2),d))
    return corners_α_w(sk, D, evv::Edgevv, b, wvalue) / (abs(w)*vd)
end

export volume_wr
function volume_wr(sk, D, b, wvalue, rvalue)
    vol = 0.0
    for i in edges(sk)
        evv = Edgevv(i.src,i.dst)
        vol += volume_edge(sk, D, evv::Edgevv, b, wvalue, rvalue)
    end
    return vol
end

export volume_constant_var
# volume for homogeneous configurations w, r with boundary wfix, rfix
function volume_constant_var(sk,D,b,levv,wfix,rfix,w,r)
    #levv = l_evv_vertices(sk,lv)
    le = length(levv)
    wvar_vec = fill(w,le)
    wval = set_value(wfix::SparseMatrixCSC, set_value(sk, levv, wvar_vec)) 
    rvar_vec = fill(r,le)
    rval = set_value(rfix::SparseMatrixCSC, set_value(sk, levv, rvar_vec))
    return volume_wr(sk, D, b, wval, rval)
end

using Roots

export solve_w_volfix
# solve for w given total volume vol and r
# homogeneous interior configuration, boundary wfix, rfix
function solve_w_volfix(sk,D,b,levv,wfix,rfix,vol,r,wi,wf)
# wi, wf - range for solver
    f(w) = volume_constant_var(sk,D,b,levv,wfix,rfix,w,r) - vol
    find_zero(f, (wi, wf))
end

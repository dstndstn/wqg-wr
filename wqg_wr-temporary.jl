# coefficient of the w_evv^(-1) term in the total volume
function co_volume_edge_w(sk, D, evv::Edgevv, b, wvalue, rvalue)
    d = D-1
    v = 0
    for i in corners(sk, D, evv::Edgevv)
        wprod = 1
        f = 0
        for j in i
            v1,v2 = j.v1, j.v2
            wprod *= wvalue[v1,v2] 
            rj = rvalue[v1,v2]
            vd = real( vvmd_r( (rj+0im)^(1/2), d) )
            lcorner = length(corners(sk, D, j))
            f += ( lcorner * vd )^(-1)
        end
        f *= 1/wprod
        v += f
    end
    v *= b
    
    r = rvalue[evv.v1,evv.v2]
    vd = real( vvmd_r( (r+0im)^(1/2), d) )
    α_l = corners_α_w(sk, D, evv, 1, wvalue)
    v += α_l * vd^(-1)
end

# Input r move r_j_n. Output w move w_l_n that preserves the total volume

# j is the edge of proposed r move
# l is the edge for proposed w move
function w_new(sk, D, evv_j, evv_l, wvalue, rvalue, r_j_n)
    jv1,jv2 = evv_j.v1, evv_j.v2
    lv1,lv2 = evv_l.v1, evv_l.v2
    # make r move
    rvalue_n = copy(rvalue)
    rvalue_n[jv1,jv2] = r_j_n
    #println("r_j_n=",r_j_n)
    # change of total volume by r move    
    vol_n = volume_edge(sk, D, evv_j::Edgevv, 1, wvalue, rvalue_n)
    vol_o = volume_edge(sk, D, evv_j::Edgevv, 1, wvalue, rvalue)
    var_vol_r = vol_n - vol_o
    #println("var_vol_r=",var_vol_r)
    # alternative: var_vol_r = α_j * w_j^(-1) * (vd_j_n^(-1) - vd_j^(-1))
    co = co_volume_edge_w(sk, D, evv_l::Edgevv, 1, wvalue, rvalue_n)
    #println("co=",co)
    # new w_l:
    w_l = wvalue[lv1,lv2]
    return ( w_l^(-1) - var_vol_r / co )^(-1)
end

## data analysis

# volume at the step-th step
function data_volume(step) 
    pars = chain[step,:]
    lp = length(pars)
    W = pars[1:Int(lp/2)]
    R = pars[1+Int(lp/2):lp]
    #wvar = set_value(W)
    #rvar = set_value(R)
    #wvalue = set_value(wfix::SparseMatrixCSC, wvar::SparseMatrixCSC)
    #rvalue = set_value(rfix::SparseMatrixCSC, rvar::SparseMatrixCSC)    
    wvalue = set_value(wfix::SparseMatrixCSC, set_value(sk, levv, W))
    rvalue = set_value(rfix::SparseMatrixCSC, set_value(sk, levv, R))    
    return volume_wr(sk, D, b, wvalue, rvalue)
end

# volumes from step1 to step2
function l_volume(step1, step2)  
    list = []
    for j in step1: step2
        push!(list, data_volume(j))
    end
    return list
end

function tstvolume(x)
    return data_volume(Int(x))
end

# display mcmc data of "chain" on edge j
function mcmcdata1(sk, j, levv, chain, ndim, v1i, v1e, v1n, v2i, v2e, v2n) 
# j-th edge; 
# 2d histogram parameters: vi - initial value for v, ve - ending value for v, vn - numbers of slots for v
    evv = levv[j]
    vtx1, vtx2 = evv.v1, evv.v2
    v1 = chain[:,j]
    v2 = chain[:,j+Int(ndim/2)]
    hgram1(j)=histogram(v1, xlabel="v1", ylabel="Number of MCMC samples",title="edge no.=$j, Edgevv=($vtx1, $vtx2)")
    hgram2(j)=histogram(v2, xlabel="v2", ylabel="Number of MCMC samples")
    display(hgram1(j))
    display(hgram2(j))
    display(histogram2d(v1, v2, xlabel="v1", ylabel="v2", bins=(range(v1i, v1e, length=v1n), range(v2i, v2e, length=v2n))))
end

function plotchain1(j,chain,sk,levv)
    plot()
    evv = levv[j]
    vtx1, vtx2 = evv.v1, evv.v2
    n = length(chain[:,j])
    plot!(chain[1:n, j],label="",title="1st variable, edge no.=$j, Edgevv=($vtx1, $vtx2)")#title="j=$j")
    plot!()
end

function plotchain2(j,chain,sk,lv)
    plot()
    evv = levv[j]
    vtx1, vtx2 = evv.v1, evv.v2
    ndof=length(levv)
    n = length(chain[:,j+ndof])
    plot!(chain[1:n, j+ndof],label="",title="2st variable, edge no.=$j, Edgevv=($vtx1, $vtx2)")#title="j=$j")
    plot!()
end
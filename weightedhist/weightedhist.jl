module weightedhist
export whist

import Base.sturges

function whist!{HT}(h::AbstractArray{HT}, v::AbstractVector, edg::AbstractVector; init::Bool=true, weights::AbstractVector 
= ones(HT,length(v)), density::Bool=false)

    length(weights) == length(v) || error("length(weights) must equal length(v)")
    n = length(edg) - 1
    length(h) == n || error("length(h) must equal length(edg) - 1.")
    
    if init
        fill!(h, zero(HT))
    end
    
    for j=1:length(v)
        i = searchsortedfirst(edg, v[j])-1
        if 1 <= i <= n
            h[i] += weights[j]
        end
    end
    
    if density == true
        binwidths = [(edg[i+1] - edg[i]) for i in 1:n]
        h = h/sum(h)./binwidths
    end
    
    edg, h
end

whist(v::AbstractVector, edg::AbstractVector; weights::AbstractVector = ones(Int,length(v)), density::Bool=false) = 
whist!(Array(eltype(weights), length(edg)-1), v, edg; weights=weights, density=density)
whist(v::AbstractVector, n::Integer; weights::AbstractVector = ones(Int,length(v)), density::Bool=false) = whist(v, 
whistrange(v,n); weights=weights, density=density)
whist(v::AbstractVector; weights::AbstractVector = ones(Int,length(v)), density::Bool=false) = whist(v, 
sturges(length(v)); weights=weights, density=density)


end # module

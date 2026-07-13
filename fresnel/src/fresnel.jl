module fresnel
using FastGaussQuadrature

greet() = print("Name: Ori Gonen")

"""
    in this work we solve fresnel (cosine) integral
    We divide the calculation, using taylor expansion for certain values and gauss quadrature for others
"""

const X_cutoff = 1
const TAYLOR_DEG = 15 # picked this value from trial and error to see if it passes the tests


"""
    y = fresnel_C(x)

computes the fresnel C integral for value x
for negative x, computes -fresnel_C(-x)
divides x inputs for 0<x<X_cutoff and X_cutoff<=x
"""
function fresnel_C(x)
    # we know that C(x) = -C(-x)
    if x<0
        return -fresnel_C(-x)
    end 
    if x==0
        return 0.0
    end

    # now we assume x > 0
    """
    strategy:
    for small x (x < X_cutoff): it suffices to use taylor expansion
    for big x (x >= X_cutoff) : we use the auxillary functions
    I choose the constant X_cutoff to be 1, seems like a fair choice
    """
    if x < X_cutoff
        return taylor_C(x)
    else
        return aux_C(x)
    end
end



"""
    y=taylor_C(x)

computes fresnel C integral with taylor expansion approximation (for small x)

C(x) = sum_{k=0}^{inf} frac{ (-1)^k (π/2)^{2k} x^{4k+1} }{ (2k)! (4k+1) }

for efficiency, we'll calculate the summands recursively
    s_k = s_{k-1} * -(pi/2)^2*x^4*(4k-3)/( 2k*(2k-1)*(4k+1) )
"""
function taylor_C(x; taylor_deg=TAYLOR_DEG)
    # calculate the terms recursively
    result = x
    
    # enough to calculate only once
    minus_pi_over_2_squared_x4 = -(pi/2)^2*x^4
    s = x

    for k in 1:taylor_deg
        s *= minus_pi_over_2_squared_x4 * (4k-3)/( 2k*(2k-1)*(4k+1) )
        result += s
    end

    return result
end


"""
    y = aux_C(x)

for big x, compute the fresnel C integral with the auxillary functions
"""
function aux_C(x)
    return 0.5 + frensel_f(x)*sin(pi*x^2/2) - frensel_g(x)*cos(pi*x^2/2)
end


"""
We calculate gauss quadrature only once instead of every time we call the function
    instead of interval [-1,1], we map it then to [0, W]
as discussed in the lab, I can do it with gausslegendre from FastGaussQuadrature
"""
const N_quad = 40 # picked this value from trial and error to see if it passes the tests
const W = 6 # picked this value from trial and error to see if it passes the tests
const unshifted_nodes, unshifted_weights = gausslegendre(N_quad) # compute for [-1,1]
# map to [0,W]
const nodes = W .* (unshifted_nodes .+ 1) ./ 2
const weights = W/2 .* unshifted_weights
# would be shorter to compute the function values of those nodes in advance
f_nodes_exp = exp.(-pi .* (nodes.^2) ./ 2) # for f
g_nodes_exp = (nodes.^2) .* f_nodes_exp # for g


"""
    f(z) = sqrt(2)/(pi*z) * integral e^(-pi w^2 / 2) / (1+ (w/z)^4) dw

    compute the f auxillary function
"""
function frensel_f(z)
    int=0.0
    for i=1:length(nodes)
        int += weights[i] * f_nodes_exp[i] / (1 + (nodes[i]/z)^4 )
    end

    return sqrt(2)/(pi*z) * int
end


"""
    g(z) = sqrt(2)/(pi*z^3) * integral w^2 e^(-pi w^2 / 2) / (1+ (w/z)^4) dw
    
    compute the g auxillary function
"""
function frensel_g(z)
    int=0.0
    for i=1:length(nodes)
        int += weights[i] * g_nodes_exp[i] / (1 + (nodes[i]/z)^4 )
    end

    return sqrt(2)/(pi*z^3) * int
end



export fresnel_C
end # module fresnel

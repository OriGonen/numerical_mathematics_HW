module HW01
using Plots

greet() = print("Name: Ori Gonen\n")

# ---------------------------------
#   cubic spline
# ---------------------------------

# struct for the spline
struct Zlepek
    x::Vector{Float64} # list of x_i
    y::Vector{Float64} # list of y_i
    coefficients::Vector{Vector{Float64}} # list of coefficients (a_i,b_i,c_i,d_i) for each interval
end

# interpolate
"""
    Z = interpoliraj(x, y)
Computes the cubic spline of points (x_i, y_i).
Returns a corresponding Zlepek object.
"""
function interpoliraj(x, y)
    """
    Theoretical computation:

        S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
        S_i'(x) = b_i + 2*c_i*(x-x_i) + 3*d_i*(x-x_i)^2
        S_i''(x) = 2*c_i + 6*d_i*(x-x_i)

        
        S''(x1) = 0
        => 2*c_1+*6d_1*(x_1-x_1)=0
        => c_1 = 0
        ------------ 1 equation
        S''(x_n) = 0
        => 2*c_n+6*d_n*( x_n - x_{n-1} )=0
        => c_n+3*d_n*( x_n - x_{n-1} )=0
        ------------ 1 equation

        S_i(x_i) = y_i
        => a_i = y_i
        -------------- n-1 equations | i=1,...,n-1
        S_i(x_{i+1}) = y_{i+1}
        let h_i = x_{i+1} - x_i
        => a_i + b_i*h_i + c_i*h_i^2 + d_i*h_i^3 = y_{i+1}
        --------------- n-1 equations | i=1,...,n-1

        S_i'(x_{i+1}) = S_{i+1}'(x_{i+1})
        => b_i + 2*c_i*h_i + 3*d_i*h_i^2 = b_{i+1}
        => b_i + 2*c_i*h_i + 3*d_i*h_i^2 - b_{i+1} = 0
        --------------- n-2 equations | i=2,...,n-1

        S_i''(x_{i+1}) = S_{i+1}''(x_{i+1})
        => 2*c_i + 6*d_i*h_i = 2*c_{i+1}
        => c_i + 3*d_i*h_i - c_{i+1} = 0
        --------------- n-2 equations | i=2,...,n-1

        total: 1+1+(n-1)+(n-1)+(n-2)+(n-2)=4(n-1)
        => 4(n-1) equations. just like the number of variables!
        
    """
    
    # make x and y sorted
    perm = sortperm(x)
    x = x[perm]
    y = y[perm] # make sure each y value stays along the corresponding x value

    n = length(x)
    # there are n-1 intervals and each one has 4 coefficients, so total 4(n-1) variables
    
    # transformation for indices
    idx(i) = 4*(i-1)
    """
    now:
    a_i = coeffs_vec[idx(i)+1]
    b_i = coeffs_vec[idx(i)+2]
    c_i = coeffs_vec[idx(i)+3]
    d_i = coeffs_vec[idx(i)+4]
    """

    # create the matrix (A) and column vector (b) as we described in the beginning
    m = 4*(n-1)
    A = zeros(m,m)
    b=zeros(m)
    
    row=1 # row counter
    
    # c_1 = 0
    A[row, 3+idx(1)] = 1
    row+=1
    
    # c_n+3*( x_n - x_{n-1} ) * d_n = 0
    h = x[n]-x[n-1]
    A[row, idx(n-1)+3] = 1
    A[row, idx(n-1)+4] = 3*h
    row+=1

    # a_i = y_i     |  i=1,...,n-1
    for i=1:(n-1)
        k=idx(i)
        A[row, k+1] = 1
        b[row] = y[i]
        row+=1
    end

    # a_i + b_i*h_i + c_i*h_i^2 + d_i*h_i^3 = y_{i+1} | i=1,...,n-1
    for i=1:(n-1)
        k=idx(i)
        h = x[i+1]-x[i]
        A[row, k+1] = 1
        A[row,k+2] = h
        A[row,k+3] = h^2
        A[row,k+4] = h^3
        b[row] = y[i+1]
        row+=1
    end    

    # b_i + 2*c_i*h_i + 3*d_i*h_i^2 - b_{i+1} = 0  | i=1,...,n-2
    for i=1:(n-2)
        k=idx(i)
        h = x[i+1]-x[i]
        A[row,k+2] = 1
        A[row,k+3] = 2*h
        A[row,k+4] = 3*h^2
        A[row, idx(i+1)+2] = -1
        row+=1
    end


    # c_i + 3*d_i*h_i - c_{i+1} = 0 | i=1,...,n-2
    for i=1:(n-2)
        k=idx(i)
        h = x[i+1]-x[i]
        A[row,k+3] = 1
        A[row,k+4] = 3*h
        A[row, idx(i+1)+3] = -1
        row+=1
    end

    # solve the system to get coefficients
    coeffs_vec = A \ b
    # change the vec into list of lists
    # such that     coeffs[k] = [a_k,b_k,c_k,d_k]
    coeffs = Vector{Vector{Float64}}(undef, n-1)
    for i=1:(n-1)
        coeffs[i] = coeffs_vec[idx(i)+1:idx(i)+4]
    end

    # create the spline instance
    spline = Zlepek(x,y, coeffs)

    return spline 

end


"""
    y = vrednost(Z, x; i)
Given a spline `Z` and point `x`, it evaluates the y level of point `x` according the spline.

has optional `i` input, to be the index of the interval `x` is located in. used to speed up
computation for the `plot_spline` function.
if `i` is `nothing` (default value), the function will find the appropriate index by itself
"""
function vrednost(Z::Zlepek, x::Float64 ; i=nothing)
    # by construction the points are already sorted by x value
    
    # if i is not given, find it
    if i === nothing
        i=length(Z.x)-1
        while x<Z.x[i]
            i = i-1
        end
    end

    # we construct the cubic polyinomial number i
    a,b,c,d = Z.coefficients[i]
    x_shift = x-Z.x[i]
    return a + b*x_shift + c*x_shift^2 + d*x_shift^3
end

""" 
    plot_spline(Z)

Plots the spline `Z`. Uses funciton `vrednost` with specified index for evaluation of points.
"""
function plot_spline(Z::Zlepek)
    x = Z.x
    coeffs = Z.coefficients

    p = Plots.plot()

    for i in 1:length(coeffs)

        x_left = x[i]
        x_right = x[i+1]

        xs = range(x_left, x_right, length=100)
        
        # using our evaluation function
        ys = [
            vrednost(Z, xi ; i=i) for xi in xs
        ]

        color = isodd(i) ? :red : :blue

        Plots.plot!(p, xs, ys, color=color, label=nothing)
    end
    scatter!(p, x, Z.y, label=nothing)

    return p
end


export Zlepek, interpoliraj, vrednost, plot_spline
end # module HW01

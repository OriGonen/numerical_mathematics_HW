module HW01

greet() = print("Name: Ori Gonen")

# ---------------------------------
#   cubic spline
# ---------------------------------

# struct for the spline
struct Zlepek
    points # list of (x_i, f_i)
    coefficients # list of coefficients (a_i,b_i,c_i,d_i) for each interval
end

# interpolate
# TODO: write description
function interpoliraj(x, y)
    """
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
    n = length(x)
    
    # there are n-1 intervals and each one has 4 coefficients, so total 4(n-1) variables
    
    # transformation for indices
    idx(i) = 4*(i-1)
    """
    now:
    a_i = u[idx(i)+1]
    b_i = u[idx(i)+2]
    c_i = u[idx(i)+3]
    d_i = u[idx(i)+4]
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
    A[row, idx(n)+3] = 1
    A[row, idx(n)+4] = 3*h
    row+=1


    for i=1:(n-1)
        h = x[i+1]-x[i]
        k=idx(i)

    end
        
    
end


end # module HW01

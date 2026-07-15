module limitcycle

greet() = print("Name: Ori Gonen\n")

"""

In this homework we will find the period of a limit cycle, for the differential equation

x''(t) - 4(1-x^2) x'(t) + x = 0

target precision: accurate to 10 decimal places


Steps:
1. reduce to first order by: u1 = x     u2 = x'
                u1' = u2
                u2' = 4(1-u1^2) u2 - u1

2. solve the equation (RK4)


3. the solution will converge into the limit cycle. detect period by checking sign alteration (- to +)

4. computing multiple times until consecutive estimates agree to 1e-10

"""






"""
    du = diff_eq(u)

returns the derivative of the vector u in accordance to the differential equation
"""
function diff_eq(u)
    return [u[2], 4(1-u[1]^2) * u[2] - u[1]]
end


"""
    u_next = RK4_step(u,h)

performs RK4 with stepsize h on vector u
"""
function RK4_step(u, h)
    k1 = diff_eq(u)
    k2 = diff_eq(u .+ (h/2).*k1)
    k3 = diff_eq(u .+ (h/2).*k2)
    k4 = diff_eq(u .+ h.*k3)

    return u .+ h .* (k1/6 + k2/3 + k3/3 + k4/6)
end



"""
    t = precise_crossing(u_init, interval)

uses bisection method to get the time of crossing (velocity changing sign - to +)
"""
function precise_crossing(u_init, interval; iters = 40)
    a, b = interval
    t_init = a
    for i=1:iters
        # find with bisection
        mid = (a+b)/2
        u_mid = RK4_step(u_init, mid-t_init) # get the new u approximation. the initial interval should be already small enough for it to be ok going with one step
        
        # we check u[2] changing sign from - to +
        if u_mid[2] < 0
            a = mid # we did not reach the crossing yet
        else
            b = mid # we overshot the crossing
        end
    end
    
    # approximation for the crossing time: middle of our resulting interval
    return (a+b)/2
end



"""
    period = estimate_period(u0, h)

compute the period of the limit cycle. u0 is the initial guess, h is step size
"""
function estimate_period(u0, h; t_max = 40)
    prev_period_est = -1.0
    prev_crossing = -1.0
    u = u0
    t = 0.0


    while t < t_max

        u_next = RK4_step(u,h)
        t_next = t + h

        # check if there was a crossing
        if u[2] < 0 && u_next[2] > 0
            currect_crossing = precise_crossing(u, [t,t_next])
            # check if the new period estimate is close enough to the previous one
            current_period_est = currect_crossing - prev_crossing

            if isapprox(prev_period_est, current_period_est ; atol=1e-10)
                return current_period_est
            end

            # update 
            prev_period_est = current_period_est
            prev_crossing = currect_crossing
        end
        # update
        t = t_next
        u = u_next
        
    end

    error("could not find sufficiently accurate period")
end





export estimate_period, RK4_step, diff_eq




end # module limitcycle

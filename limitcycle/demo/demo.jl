using limitcycle

limitcycle.greet()



function find_limit_cycle_period(u0, h0 ; max_iters=20)
    prev_estimate = -2.0 # initialize to impossible values
    current_estimate = -1.0 # initialize to impossible values
    h = h0
    iter = 1

    # we will halve h until we get a sufficiently close answer
    # we want to make sure that the estimates from h and h/2 are close enough
    while iter < max_iters
        current_estimate = estimate_period(u0, h)
        if isapprox(current_estimate, prev_estimate; atol=1e-10)
            println("period estimate:")
            println(current_estimate)
            println("estimated with step size ", h)
            return current_estimate
        else
            prev_estimate = current_estimate
            h /= 2
        end
    end

    return current_estimate
end

period_est = find_limit_cycle_period([1.0, 0.0], 1e-3)
# I got more than 10 decimal places of accuracy
# period estimation: 10.203523691042903




using Plots

function integrate_trajectory(u0, h, T)
    n = Int(round(T / h))
    t = zeros(n + 1)
    traj = Matrix{Float64}(undef, n + 1, 2)
    traj[1, :] .= u0

    u = u0
    for i in 1:n
        u = RK4_step(u, h)
        traj[i+1, :] .= u
        t[i+1] = t[i] + h
    end

    return t, traj
end

h = 1e-3
T_max = 40.0
t, traj = integrate_trajectory([1.0, 0.0], h, T_max)
x = traj[:, 1]
xprime = traj[:, 2]

# (x, x') converging to the limit cycle
plot(x, xprime;
    xlabel = "x", ylabel = "x'",
    title  = "Convergence to the Limit Cycle",
    legend = false,
    linewidth = 1)
scatter!([x[1]], [xprime[1]]; color = :green, markersize = 6)
savefig("xprime_vs_x.png")

# x(t) over the last few periods, once settled onto the limit cycle
last3periods = t .>= T_max - 3 * period_est
plot(t[last3periods], x[last3periods];
    xlabel = "t", ylabel = "x(t)",
    title  = "x(t) on the Limit Cycle",
    legend = false,
    linewidth = 2)
savefig("x_vs_t.png")



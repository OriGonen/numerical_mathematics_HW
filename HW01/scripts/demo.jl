using HW01
using Plots
using Random
Random.seed!(1)

HW01.greet()

# function to print the examples in a plot
function demo_examples()
    examples = [

        ("Mixed function",
            -100:20:100,
            x -> -5*log(x^2 + 1) + 4*sin(x),
            true
        ),

        ("Random",
            0:1:20,
            x -> x + 20*randn(),
            false # no "real" function to compare to
        ),

        ("Sin with uneven spacing",
            [0.0, 1.0, 2.0, 10.0, 50.0, 100.0],
            x -> sin(x),
            true
        ),

        ("High frequency Sin (undersampled)",
            0:1/4:1,
            x -> sin(10*pi*x),
            true
        ),

        ("|x| with uneven spacing",
            [-5, -4, -3.5, -0.5, 1.5, 3, 3.3, 5],
            x -> abs(x),
            true
        ),

        ("Gaussian",
            -10:2:10,
            x -> exp(-10*x^2),
            true
        )
    ]

    plots = []

    for (name, xs, f, has_true_function) in examples
        ys = [f(x) for x in xs]
        Z = interpoliraj(xs, ys)

        p = plot_spline(Z)

        if has_true_function
            x_real = range(minimum(xs), maximum(xs), length=1000)
            y_real = [f(x) for x in x_real]

            plot!(p, x_real, y_real,
                label="True function",
                linewidth=2)
        end

        scatter!(p, xs, ys, label="Data points")
        title!(p, name)

        push!(plots, p)
    end

    plot(plots..., layout=(3,2), size=(1000,800))
end



function undersampling_demo()
    f(x) = sin(10*pi*x)

    samplings = [
        ("0.02", 0:0.02:1),
        ("0.05", 0:0.05:1),
        ("0.10", 0:0.10:1),
        ("0.25", 0:0.25:1)
    ]

    plots = []

    for (labeltxt, xs) in samplings
        ys = [f(x) for x in xs]
        Z = interpoliraj(xs, ys)

        p = plot_spline(Z)

        x_real = range(0, 1, length=2000)
        plot!(p, x_real, [f(x) for x in x_real],
            label="sin(10*pi*x)",
            linewidth=2)

        scatter!(p, xs, ys, label="Samples")
        title!(p, "interval=$labeltxt")

        push!(plots, p)
    end

    plot(plots..., layout=(2,2), size=(1000,700))
end


demo_examples()
savefig("demo.png")

undersampling_demo()
savefig("undersampling_demo.png")


using HW01
using Plots
using Random
Random.seed!(1)


HW01.greet()

function demo_examples()
    examples = [

        ("Mixed function",
            -100:20:100,
            x -> -5*log(x^2 + 1) + 4*sin(x)
        ),

        ("Random",
            0:1:20,
            x -> x + 20*randn()
        ),

        ("Sin with uneven spacing",
            [0.0, 1.0, 2.0, 10.0, 50.0, 100.0],
            x -> sin(x)
        ),

        ("High frequency Sin",
            0:1/4:1,
            x -> sin(10*pi*x)
        ),

        ("|x| with uneven spacing",
            [-5, -4, -3.5, -0.5, 1.5, 3, 3.3, 5],
            x -> abs(x)
        ),

        ("Gaussian",
            -10:2:10,
            x -> exp(-10*x^2)
        )
    ]

    plots = []

    for (name, xs, f) in examples
        ys = [f(x) for x in xs]
        Z = interpoliraj(xs, ys)

        p = plot_spline(Z)
        title!(p, name)

        push!(plots, p)
    end

    plot(plots..., layout=(3,2), size=(900,800))
end

demo_examples()


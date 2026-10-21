using fresnel
using Plots


fresnel.greet()


x = -10:0.001:10
plot(x, fresnel_C.(x);
    xlabel = "x", ylabel = "C(x)",
    title  = "Fresnel Cosine Integral",
    legend = false,
    linewidth = 2)
# assymptotes at -0.5 and 0.5
hline!([0.5, -0.5]; linestyle = :dash, color = :gray)

savefig("Fresnel Cosine Integral")


# zoom in on X_cutoff
x_zoom = 0.9:0.0001:1.1
plot(x_zoom, fresnel_C.(x_zoom);
    xlabel="x", ylabel="C(x)",
    title="Zoom near X_cutoff",
    legend=false,
    linewidth=2)
vline!([1]; linestyle = :dash, color = :gray)
savefig("Zoom in near X_cutoff")
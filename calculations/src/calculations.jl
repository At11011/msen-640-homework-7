module calculations

using CairoMakie, ForwardDiff

function problem_2_a()
    R = 8.314
    T = 600
    G(x) = @. 8400 * x * (1 - x) + R * T * (x * log(x) + (1 - x) * log(1 - x))
    dG(x) = @. ForwardDiff.derivative(G, x)

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="a_b", ylabel="G")

    x = range(0, 1, 100)
    g = G(x)
    dg = dG(x)

    aB = @. exp((g + (1 - x) * dg) / (R * T))

    lines!(ax, x, aB)

    save("../assets/fig_1.png", fig)

    fig
end

function problem_2_b()
    R = 8.314
    T = 600
    G(x) = @. 8400 * x * (1 - x) + R * T * (x * log(x) + (1 - x) * log(1 - x))
    dG(x) = @. ForwardDiff.derivative(G, x)

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="a_b", ylabel="G")

    x = range(0, 1, 100)
    g = G(x)
    dg = dG(x)

    Gfus = -1200

    aB = @. exp((g + (1 - x) * dg) / (R * T)) * exp(Gfus / (R * T))

    lines!(ax, x, aB)

    save("../assets/fig_2.png", fig)

    fig
end

function main()
    problem_2_b()
end

end



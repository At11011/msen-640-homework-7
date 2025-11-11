module calculations

using CairoMakie, ForwardDiff, NLsolve

function problem_1()
    # constants and parameters
    R = 8.314462618
    Tm1 = 1550.0   # K, component 1
    Tm2 = 1200.0   # K, component 2
    dS1 = 8.0      # J/mol/K
    dS2 = 11.0     # J/mol/K
    omega_alpha = 19500.0  # J/mol
    omega_liq = 1500.0   # J/mol

    deltaG_mix(omega, x, T) = omega * x * (1 - x) + R * T * (x * log(x) + (1 - x) * log(1 - x))
    function mu_B_of_G(omega, x, T)
        x = clamp(x, 1e-9, 1 - 1e-9)
        (omega * (1 - 2 * x) + R * T * log(x / (1 - x)) + 0.0)
    end
    G_alpha(x, T) = deltaG_mix(omega_alpha, x, T)
    function G_liq(x, T)
        dG1 = dS1 * (Tm1 - T)
        dG2 = dS2 * (Tm2 - T)
        return (1 - x) * dG1 + x * dG2 + deltaG_mix(omega_liq, x, T)
    end
    function muB_alpha(x, T)
        x = clamp(x, 1e-9, 1 - 1e-9)
        omega_alpha * (1 - 2 * x) + R * T * log(x / (1 - x))
    end
    function muB_liq(x, T)
        dG1 = dS1 * (Tm1 - T)
        dG2 = dS2 * (Tm2 - T)
        return (dG2 - dG1) + omega_liq * (1 - 2 * x) + R * T * log(x / (1 - x))
    end
    muA_from_G_muB(G, x, muB) = G - x * muB


    function find_tie_line(T)
        function equations!(F, vars)
            xa, xl = vars
            # constrain in (0,1)
            xa = clamp(xa, 1e-9, 1 - 1e-9)
            xl = clamp(xl, 1e-9, 1 - 1e-9)

            muB_a = muB_alpha(xa, T)
            muB_l = muB_liq(xl, T)
            Ga = G_alpha(xa, T)
            Gl = G_liq(xl, T)
            muA_a = muA_from_G_muB(Ga, xa, muB_a)
            muA_l = muA_from_G_muB(Gl, xl, muB_l)

            F[1] = muB_a - muB_l
            F[2] = muA_a - muA_l
        end
        guesses = [
            [0.05, 0.95], [0.1, 0.9], [0.2, 0.8], [0.4, 0.6], [0.5, 0.7]
        ]

        for g in guesses
            try
                sol = nlsolve(equations!, g, autodiff=:forward)
                if converged(sol)
                    xa, xl = sol.zero
                    if 0.0 < xa < 1.0 && 0.0 < xl < 1.0 && abs(xa - xl) > 1e-6
                        return xa, xl
                    end
                end
            catch
                continue
            end
        end

        return nothing
    end

    Ts = [800, 1000, 1200, 1300, 1450, 1500, 1550, 1650]
    xs = range(0, 1, 100)
    for T in Ts
        Ga = G_alpha.(xs, T)
        Gl = G_liq.(xs, T)
        tie = find_tie_line(T)

        fig = Figure()
        ax = Axis(fig[1, 1], xlabel="Mole fraction x_2", ylabel="G (kJ/mol)", title="G vs T at $T K")
        alpha = lines!(ax, xs, Ga / 1000.0)
        liq = lines!(ax, xs, Gl / 1000.0)

        if tie !== nothing
            xa, xl = tie
            # draw tie points
            Ga_xa = G_alpha(xa, T)
            Gl_xl = G_liq(xl, T)
            scatter!(ax, [xa], [Ga_xa / 1000.0], color=:black, marker=:rect, markersize=10)
            scatter!(ax, [xl], [Gl_xl / 1000.0], color=:black, marker=:rect, markersize=10)

            # draw common tangent line (slope = muB at those points)
            slope = muB_alpha(xa, T)
            # tangent line passing through (xa, Ga)
            xline = range(0, 1, 3)
            tangent = (slope .* (xline .- xa) .+ Ga_xa)
            t = lines!(ax, xline, tangent / 1000.0, color=:black, linestyle=:dash, alpha=0.7)

            println("T=$T K: tie-line found x_alpha=$(round(xa, digits=5)), x_liquid=$(round(xl, digits=5))")
            Legend(fig[1, 2], [alpha, liq, t], ["Gᵅ (kJ/mol)", "Gˡ (kJ/mol)", "Tie line"])
        else
            println("T=$T K: no tie-line found (single phase or solver failed).")
            Legend(fig[1, 2], [alpha, liq], ["Gᵅ (kJ/mol)", "Gˡ (kJ/mol)"])
        end

        save("../assets/G_vs_x_T_$T.png", fig)
    end


    Ts = range(1200, 1540, 20)
    xas = []
    xls = []
    for T in Ts
        tie = find_tie_line(T)
        if tie !== nothing
            xa, xl = tie
            append!(xas, xa)
            append!(xls, xl)
        end
    end
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="T", ylabel="x_2")
    a = scatter!(ax, Ts, xas)
    l = scatter!(ax, Ts, xls)
    Legend(fig[1, 2], [a, l], ["Solid phase boundary", "Liquid phase boundary"])
    save("../assets/phase_diagram.png", fig)
    fig

end

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

function problem_3()
    a = 10500
    b = 3e-4
    c = 8e-10
    R = 8.314
    T = 298
    n_a = 35
    n_b = 65
    X_A = n_a / (n_a + n_b)
    X_B = n_b / (n_a + n_b)
    P = 101325

    GE = a * (1 - b * T) * (1 - c * P) * X_A * X_B
    GID = R * T * (X_A * log(X_A) + X_B * log(X_B))
    println((GE + GID) * (n_a + n_b))
end

function main()
    problem_3()
end

end



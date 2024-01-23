### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 884380c2-2a81-11ee-35db-2bc32422997c
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ aedd096b-8103-4f07-9cbd-827088beab6f
using DifferentialEquations,
    MethodOfLines, ModelingToolkit, DomainSets, Plots, NonlinearSolve

# ╔═╡ a5e96864-549e-49ec-8d40-9ee3bd8a5756
md"# Partial Differential Equations"

# ╔═╡ 2ae423f2-daf3-4ebe-ac41-d9347be427d4
begin
    @parameters t, x, y
    @variables T(..), u(..), q(..)
end;

# ╔═╡ db1e018b-03a8-4746-8040-37848b084f77
md"## Example 1

$k\frac{\partial^2 T}{\partial x^2} - \rho c_p \frac{\partial T}{\partial t} = 0$

- `T` : temperature
- `ρ` : density
- `cₚ`: key heat conductivity component
- `k` : thermal conductivity
"

# ╔═╡ 5c3b9561-0a96-40f6-865f-a606082c64de
begin
    Dt = Differential(t)
    Dx = Differential(x)
    Dtt = Dt^2
    Dxx = Dx^2


    k = 666
    L = 0.25
    ρ = 1000
    cₚ = 500
    # q = 5e6

    # eq_1 = k * Dxx(T(t, x)) ~ ρ * cₚ * Dt(T(t, x)) + q
    eq_1 = [k * Dxx(T(t, x)) ~ ρ * cₚ * Dt(T(t, x))] #+ k * Dx(T(t, x))]


    domain_1 = [x ∈ Interval(0.0, L), t ∈ Interval(0.0, 5.0)]

    bcs_1 = [
        # Inital condition
        T(0, x) ~ 0.0,
        Dt(T(t, 0)) ~ 680,
        T(t, L) ~ 0.0,
        # q(0) ~ 0
    ]

    @named wall = PDESystem(eq_1, bcs_1, domain_1, [t, x], [T(t, x)])

    wall_discretization = MOLFiniteDifference([x => 0.01], t)
    prob_1 = discretize(wall, wall_discretization)
    wall_model = solve(prob_1; saveat = 0.5)
    heat_arcoss_wall = wall_model[T(t, x)]
end;

# ╔═╡ 9f2252c9-72d5-47ba-82b9-7e3927e0f048
plot(
    wall_model[x],
    heat_arcoss_wall[end, :];
    label = nothing,
    xlabel = "Distance (m)",
    ylabel = "Temperature (°C)",
    ylim = (0, 700),
)

# ╔═╡ 28cad6ea-064c-4774-851d-69ff94a6461c
wall_model[T(t, x)]

# ╔═╡ 86091c22-d19e-49b2-8574-8b59326c5bb8
heat_arcoss_wall[end, :]

# ╔═╡ aed303d8-faff-4906-b5cc-3b6b1dc45ebb
# gif(anim_1; fps=10)


# ╔═╡ b7c252fc-3b81-48b0-b34a-8407cca31976
md"## Example 2
$\frac{\partial^2 T}{\partial x^2} = - \frac{\partial^2 T}{\partial y^2}$
"

# ╔═╡ 6f4ec480-1bd3-4ae9-9e36-b8444a46a06b
# sol = solve(prob_2)

# ╔═╡ 3ba4a33f-f549-4d22-b2fb-2fbe0f07ec1b
md"## Example 3
$\frac{\partial^2 u(x,t)}{\partial t^2} = \frac{K L^2}{M} \frac{\partial^2 u(x,t)}{\partial x^2}$
"

# ╔═╡ 8599e14e-1fae-4dcc-9eea-2a0bf3dac11f
begin
    Dyy = Differential(y)^2

    eq = Dxx(u(x, y)) + Dyy(u(x, y)) ~ 0

    temperature_left = 100
    temperature_right = 75
    temperature_top = 250
    temperature_bottom = 300

    width = 1
    height = 1

    bcs = [
        u(0.0, y) ~ temperature_left,
        u(width, y) ~ temperature_right,
        u(x, 0) ~ temperature_bottom,
        u(x, height) ~ temperature_top,
    ]

    # Space and time domains
    domains = [x ∈ Interval(0.0, width), y ∈ Interval(0.0, height)]

    @named pdesys = PDESystem([eq], bcs, domains, [x, y], [u(x, y)])

    dx = 0.1
    dy = 0.1

    # Note that we pass in `nothing` for the time variable `t` here since we
    # are creating a stationary problem without a dependence on time, only space.
    discretization = MOLFiniteDifference([x => dx, y => dy], nothing, approx_order = 2)

    prob = discretize(pdesys, discretization)
    sol = NonlinearSolve.solve(prob, NewtonRaphson())

    u_sol = sol[u(x, y)]


    heatmap(
        sol[x],
        sol[y],
        u_sol,
        xlabel = "x values",
        ylabel = "y values",
        title = "Steady State Heat Equation",
    )
end

# ╔═╡ e638dfe8-ae17-4fc6-b9ff-8fa92c4e1404
plot(sol[x], sol[y], u_sol)

# ╔═╡ Cell order:
# ╟─a5e96864-549e-49ec-8d40-9ee3bd8a5756
# ╠═884380c2-2a81-11ee-35db-2bc32422997c
# ╠═aedd096b-8103-4f07-9cbd-827088beab6f
# ╠═2ae423f2-daf3-4ebe-ac41-d9347be427d4
# ╟─db1e018b-03a8-4746-8040-37848b084f77
# ╠═5c3b9561-0a96-40f6-865f-a606082c64de
# ╠═9f2252c9-72d5-47ba-82b9-7e3927e0f048
# ╠═28cad6ea-064c-4774-851d-69ff94a6461c
# ╠═86091c22-d19e-49b2-8574-8b59326c5bb8
# ╠═aed303d8-faff-4906-b5cc-3b6b1dc45ebb
# ╟─b7c252fc-3b81-48b0-b34a-8407cca31976
# ╠═6f4ec480-1bd3-4ae9-9e36-b8444a46a06b
# ╟─3ba4a33f-f549-4d22-b2fb-2fbe0f07ec1b
# ╠═8599e14e-1fae-4dcc-9eea-2a0bf3dac11f
# ╠═e638dfe8-ae17-4fc6-b9ff-8fa92c4e1404

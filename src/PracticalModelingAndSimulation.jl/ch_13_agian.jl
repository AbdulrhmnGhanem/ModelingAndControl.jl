### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 0438e470-2d88-11ee-38bf-bd9d69e276e2
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ 83826d5a-2be6-453f-93b8-e6bb3d57414a
using DifferentialEquations, MethodOfLines, ModelingToolkit, DomainSets, Plots, NonlinearSolve


# ╔═╡ ca632c73-3b60-43b1-b061-ef8825c570ce
md"[The Numerical Method of Lines](https://reference.wolfram.com/language/tutorial/NDSolveMethodOfLines.html)"

# ╔═╡ 5255ad70-30dc-4be6-9c68-e03a5b87493e
begin
	@parameters t, x
	@variables T(..), q(..)
	
	Dt = Differential(t)
	Dx = Differential(x)
	DDx = Differential(x)^2
	
	k = 666
	L = 0.25
	ρ = 1000
	cₚ = 500
	q_rate = 5e6
	
	eq_1 = [
			Dt(q(t)) ~ q_rate,
			k * DDx(T(t, x)) ~ ρ * cₚ * Dt(T(t, x)) + q(t)
	]
	
	domain_1 = [
		x ∈ Interval(0.0, L),
		t ∈ Interval(0.0, 5.0)
	]
	
	bcs_1 = [
		# Inital condition
		T(0, x) ~ 0.0,
		T(t,0) ~ 0.0,
		T(t, L) ~ 0.0,
		# q(t) ~ 1.0,
	]
	
	@named sys_1 = PDESystem(eq_1, bcs_1, domain_1, [t, x], [T(t, x), q(t)])
	
	d = MOLFiniteDifference([x => 0.01], t)
	prob_1 = discretize(sys_1, d)
	sol_1 = solve(prob_1; saveat=0.5)
end;

# ╔═╡ c4e906e8-7624-4103-8707-1ff19c13753d
@gif for (i, t_desc) in enumerate(sol_1[t])
	plot(sol_1[x], sol_1[T(t, x)][i,:];
		ylim=(0, 700),
		label=nothing,
		xlabel="distance (m)",
		ylabel="temperature (°C)",
		title="t = $t_desc"
	)
end fps=10

# ╔═╡ 3ffb3191-4e9e-4994-b44b-9fca571f1157
plot(sol_1[x], sol_1[T(t, x)][end,:];
	ylim=(0, 700),
	label=nothing,
)

# ╔═╡ Cell order:
# ╠═0438e470-2d88-11ee-38bf-bd9d69e276e2
# ╠═ca632c73-3b60-43b1-b061-ef8825c570ce
# ╠═83826d5a-2be6-453f-93b8-e6bb3d57414a
# ╠═5255ad70-30dc-4be6-9c68-e03a5b87493e
# ╠═c4e906e8-7624-4103-8707-1ff19c13753d
# ╠═3ffb3191-4e9e-4994-b44b-9fca571f1157

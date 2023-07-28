### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2cf8e56c-0ec1-11ee-14e3-3bc363d524f0
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate("..")
end

# ╔═╡ 2ec1db18-6f3b-4f29-9420-e210abfe4568
using Plots, DifferentialEquations, PlutoUI, ModelingToolkit, LinearAlgebra

# ╔═╡ 0fa70b6e-2557-4206-8cb8-476dec8f15f5
md"# Higher-Order and Coupled ODEs"

# ╔═╡ 06816ad6-1080-498d-bc2f-1054539bef32
md"## Fourth-Order ODE Problem

${d^4y \over dt^4} + 3{d^3y \over dt^3} - sin(t){dy \over dt} + 8y = t^2$

${y(0)} = 50, \ dy(0) = 40, \ d^2y(0) = 30, \ d^3y(0) = 20$
"

# ╔═╡ c47e3c5e-b3b8-49f6-b5f7-1e4cea434dc3
begin
	function forth_order_ode(y, p, t)
		y1, y2, y3, y4 = y

		y = y1
		dy1 = y2
	    dy2 = y3
	    dy3 = y4
	    
		dy4 = t^2 - 3dy3 + sin(t) * dy1 - 8y
	    
		return [dy1, dy2, dy3, dy4]
	end
	
	f0 = [50, 40, 30, 20]
	tspan_1 =  (0, 2π)
	prob_1 = ODEProblem(forth_order_ode, f0, tspan_1)
	sol_1 = solve(prob_1)
end;

# ╔═╡ da7da613-c325-411d-9024-d1f20b47ad59
md"Clip y axis $(@bind clip CheckBox(default=false))"

# ╔═╡ 93c37c1f-5d3e-4fd1-b3ce-8b01aa6ca873
begin
	ylim = clip ? (-700, 300) : (Inf, Inf)
	
	plot(sol_1;
		idxs=[4],
		labels="d³y",
		ylim,
	)
	
	plot!(sol_1;
		idxs=[3],
		labels="d²y",
		ylim,
	)
	
	plot!(sol_1;
		idxs=[2],
		labels="dy",
		ylim,
	)
	
	plot!(sol_1;
		idxs=[1],
		labels="y",
		ylim,
	)
end

# ╔═╡ 0231b155-e59c-448f-b174-352a1a737076
md"## Robertson Problem"

# ╔═╡ 011c2eaa-57fd-4923-b829-ffdb521507f3
begin
	@variables t, y₁(t), y₂(t), y₃(t), y₄(t), y₅(t), y₆(t)
	D = Differential(t)
	eqs_2 = [
		D(y₁) ~ -0.04y₁ + (10^4)y₂ * y₃,
		D(y₂) ~ 0.04y₁ - (10^4)y₂ * y₃ - 3 * (10^7)y₂,
		D(y₃) ~ 3 * (10^7)y₂
	]
	@named sys_2 = ODESystem(eqs_2)
	
	f0_2 = [y₁ => 1.0, y₂ => 0.0, y₃ => 0.0]
	tspan_2 = (0, 10^2)
	
	prob_2 = ODEProblem(sys_2, f0_2, tspan_2)
	sys_2
end

# ╔═╡ 11c71279-bdba-4cc5-9346-e435221512a7
sol_2 = solve(prob_2);

# ╔═╡ cb6248a0-d655-4bfa-90bb-de1baa139515
plot(sol_2; title="Robertson Problem")

# ╔═╡ aa1eb81f-932d-4f0a-afce-5a3ca224f87c
md"## Akzo-Nobel Problem
$M{dy \over dt} = akzo\_nobel(y)$
${dy \over dt} = M^{-1}akzo\_nobel(y)$

$M = M^{-1} =$
$\begin{bmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
\end{bmatrix}$
"

# ╔═╡ 888b0642-3a9e-4241-923c-7ef18d74636e
begin
	const M = diagm([1, 1, 1, 1, 1, 0])
	const k₁   = 18.7
	const k₂   = 0.58
	const k₃   = 0.09
	const k₄   = 0.42
	const K    = 34.4
	const klA  = 3.3
	const Kₛ   = 115.83
	const pCO₂ = 0.9
	const H = 737
	
	function akzo_nobel(dy, y, p, t)
		y₁, y₂, y₃, y₄, y₅, y₆ = y
		
		r₁ = k₁ * y₁^4 * y₂^0.5
		r₂ = k₂ * y₃ * y₄
		r₃ = (k₂ / K) * y₁ * y₅
		r₄ = k₃ * y₁ * y₄^2
		r₅ = k₄ * y₆^2 * y₂^(0.5)

		Fᵢₙ = klA * ((pCO₂ / H) - y₂)

		
	    dy[1] = -2r₁ + 2r₂ - r₃ + r₄
		dy[2] = -0.5r₁ - r₄ - 0.5r₅ + Fᵢₙ
		dy[3] = r₁ - r₂ + r₃
		dy[4] = -r₂ + r₃ - 2r₄
		dy[5] = r₂ - r₃  + r₅
		dy[6] = Kₛ * y₁ * y₄ - y₆
		dy = M * dy
	end
	
	f0_3 = [0.437, 0.00123, 0, 0.007, 0, Kₛ * 0.437 * 0.007]
	tspan_3 = (0, 180)
	prob_3 = ODEProblem(akzo_nobel, f0_3, tspan_3)
	sol_3 = solve(prob_3)
end;

# ╔═╡ 600904fb-2471-4faa-85fa-cf261c1fe5f3
plot(sol_3; title="Akzo-Nobel")

# ╔═╡ 03cd0a96-eb3c-423f-9fc4-83e1f5bbfaa3
md"## HIRES Problem"

# ╔═╡ 0cd6e800-4b31-4245-8d09-eb8fe532dc87
begin
	function HIRES(du, u, p, t)
		u₁, u₂, u₃, u₄, u₅, u₆, u₇, u₈ = u

		
		du[1] = -1.7u₁ + 0.43u₂ + 8.32u₃ + 0.0007
		du[2] = 1.71u₁ - 8.75u₂
		du[3] = -10.03u₃ + 0.43u₄ + 0.035u₅
		du[4] = 8.32u₂ + 1.71u₃ - 1.12u₄
		du[5] = -1.75u₅ + 0.43u₆ + 0.43u₇
		du[6] = -280 * u₆ * u₈ + 0.69u₄ + 1.71u₅ - 0.43u₆ + 0.69u₇
		du[7] = 280 * u₆ * u₈ - 1.81u₇
		du[8] = -280 * u₆ * u₈ + 1.81u₇
	end

	f0_4 = [1, 0, 0, 0, 0, 0, 0, 0.0057]
	tspan_4 = (0, 300)
	prob_4 = ODEProblem(HIRES, f0_4, tspan_4)
	sol_4 = solve(prob_4)
end;

# ╔═╡ 6603dd5f-e7a2-463c-aadb-99548d79a568
plot(sol_4; title="HIRES")

# ╔═╡ c423edee-85ef-4642-8e58-822c7956ad73
plot(sol_4[1:50]; title="HIRES")

# ╔═╡ Cell order:
# ╟─0fa70b6e-2557-4206-8cb8-476dec8f15f5
# ╠═2cf8e56c-0ec1-11ee-14e3-3bc363d524f0
# ╠═2ec1db18-6f3b-4f29-9420-e210abfe4568
# ╟─06816ad6-1080-498d-bc2f-1054539bef32
# ╠═c47e3c5e-b3b8-49f6-b5f7-1e4cea434dc3
# ╟─da7da613-c325-411d-9024-d1f20b47ad59
# ╟─93c37c1f-5d3e-4fd1-b3ce-8b01aa6ca873
# ╟─0231b155-e59c-448f-b174-352a1a737076
# ╠═011c2eaa-57fd-4923-b829-ffdb521507f3
# ╠═11c71279-bdba-4cc5-9346-e435221512a7
# ╟─cb6248a0-d655-4bfa-90bb-de1baa139515
# ╟─aa1eb81f-932d-4f0a-afce-5a3ca224f87c
# ╠═888b0642-3a9e-4241-923c-7ef18d74636e
# ╟─600904fb-2471-4faa-85fa-cf261c1fe5f3
# ╟─03cd0a96-eb3c-423f-9fc4-83e1f5bbfaa3
# ╠═0cd6e800-4b31-4245-8d09-eb8fe532dc87
# ╟─6603dd5f-e7a2-463c-aadb-99548d79a568
# ╟─c423edee-85ef-4642-8e58-822c7956ad73

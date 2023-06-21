### A Pluto.jl notebook ###
# v0.19.25

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

# ╔═╡ 1bee3436-edd1-11ed-2fbd-873b2824b63d
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate("..")
end

# ╔═╡ 98ca666d-fca1-4e8f-950d-c66a8cb81db3
using DifferentialEquations, Plots, PlutoUI

# ╔═╡ 28a9d730-4017-4786-bc60-8fb4981b4ed1
# ╠═╡ show_logs = false
begin
	using SymPy
	import_from(sympy)
end

# ╔═╡ 32cccc2d-9eef-43c5-a441-fd112c2272e1
md"# Analytical Solutions for ODEs

!!! warning 
	Julia doesn't encourage analytcial solutions, so I am going to use `DifferentialEquations.jl` which takes a numerical approach just to get used to its API.
"

# ╔═╡ f62e7973-f1ef-47bb-b29c-b838dec186a6
md"## Example 8 - second order ODE
${d^2u \over dt^2} + 100u = 2.5sin(10t)$
$u(0)=0$
${du(0) \over dt}=0$
"

# ╔═╡ 5e48716c-f25f-42f0-a88f-3e32add73f72
function ode_8(ddu, du, u, p, t)
	@. ddu = 2.5 * sin(10 * t) - 100 * u
end

# ╔═╡ ab98d373-8dbe-43ff-8457-0ab7cf970248
begin
	u0 = [0.0]
	du0 = [0.0]
	tspan = (0.0, 10.0)
end

# ╔═╡ 5c9a55da-7b48-45c8-9f9c-d6b4d79d237d
prob = SecondOrderODEProblem(ode_8, du0, u0, tspan, []);

# ╔═╡ c6dbbfd0-d46e-4f60-9287-7ce90a6a1469
sol = solve(prob, DPRKN6(); reltol=1e-12);

# ╔═╡ 8327be9e-04fd-4cf3-bfab-f326e285474c
md"### Analytical solution"

# ╔═╡ 8db1e5f5-2591-4c69-a5e0-f9521dce75ac
begin
	u(t) = 
sin(10t) * (0.00625 - 0.00625 * cos(20t)) + (0.00625 * sin(20t) - 0.125t) * cos(10t) 
#  analytical solution using wolfram alpha https://www.wolframalpha.com/input?i=u%27%27+%3D+2.5sin%2810t%29+-+100u%28t%29%2C+u%280%29%3D0%2C+u%27%280%29%3D0

	t = range(tspan..., 1000)
	us = u.(t)
end;

# ╔═╡ 4b201856-e4af-4cd0-bb57-1af254dc26f1
begin
	plot(sol;
		idx =[2,1],
		title = "Simple Harmonic Oscillator",
		xaxis = "Time",
		yaxis = "Elongation",
		label = ["du" "u"]
	)
	plot!(t, us;
		label="u (analytical)",
		linestyle=:dash
	)
end

# ╔═╡ f799db33-a4bf-4ec8-b12d-47a2005b6418
begin
	u_vals = [sol(t)[1] for t in collect(sol.t)]
	du_vals = [sol(t)[2] for t in collect(sol.t)]
	plot(u_vals, du_vals, xlabel="u", ylabel="du", legend=false)
end

# ╔═╡ 360d55ea-5da7-47bf-99bc-2dc44ad9b0af
md"## Example 9 - System of ODEs
$\begin{align}
	y'_1 &= y_2 \\
	y_2 &= -y_1 -0.125y_2 \\
	y_1(0) &= 1 \\
	y_2(0) &= 0
\end{align}$"

# ╔═╡ a434b1fa-a4c9-4d36-baef-5d34d9abc6ac
begin
	function ode_9(du, u, p, t)
		y₁, y₂ = u
		du[1] = y₂
		du[2] = -y₁ - 0.125y₂
	end
	
	u0_9 = [1.0, 0.0]
	prob_9 = ODEProblem(ode_9, u0_9, tspan)
	sol_9 = solve(prob_9, Tsit5())
end;

# ╔═╡ 28b2bb6a-d138-45a3-84d1-e7d5c7d0d712
begin
	y₁(t) = 
(cos(t * sqrt(255) / 16) / ℯ^(t / 16)) + ((sqrt(255) * sin(t * sqrt(255) / 16)) / (255 * ℯ^(t / 16)))
y₂(t) = -((16 * sqrt(255) * sin(t * sqrt(255) / 16) / (255 * ℯ^(t/16))))
	analytcial_ys₁ = y₁.(t)
	analytcial_ys₂ = y₂.(t)
end;

# ╔═╡ a207cab4-cd88-41c2-8de6-5bdb5860bd9e
begin
	plot(sol_9;
		idx =[2,1],
		xaxis = "Time",
		label = ["y₁" "y₂"]
	)
	
	plot!(t, [analytcial_ys₁, analytcial_ys₂];
		ls=:dash,
		lw=2,
		label = ["y₁ (analytical)" "y₂ (analytical)"]
	)
end

# ╔═╡ ad0d1102-a7aa-4942-8b45-76e695662663
md"## Example 10 - second order ODE
$\begin{align}
	2ÿ + 3 ẏ^3 - |y|cos(100t) &= 2 \\
	y(0) &= 0 \\
	ẏ(0) &= 0
\end{align}$"

# ╔═╡ cbb54be0-d5bc-4755-88b7-ef7e299777be
begin
	function ode_10(ddu, du, u, p, t)
		@. ddu = (abs(u) / 2) * cos(100 * t) - (3/2) * du^3 + 1
	end
	u0_10 = [0.0]
	du0_10 = [0.0]
	prob_10 = SecondOrderODEProblem(ode_10, du0_10, u0_10, tspan, []);
end

# ╔═╡ 86948555-aa8e-43e6-91d2-74268397acc5
sol_10 = solve(prob_10, DPRKN6(); reltol=1e-12);

# ╔═╡ 133981dd-d36e-41de-a07c-0dd2ee1634cc
plot(sol_10, idx=[2,1];
	label=["du" "u"]
)

# ╔═╡ 3d425689-f115-46bc-bea1-903bb64a2100
md"## Example 11 - MuPad
$\begin{align}
ÿ+ẏ &= sin(t) \\
y(0) &= 1 \\
ẏ(0) &= 2
\end{align}$
"

# ╔═╡ 7d4f90cc-0010-4e18-9b90-1b7b86cf0c94
begin
	function ode_11(ddu, du, u, p, t)
		@. ddu = sin(t) - du
	end
	u0_11 = [1.0]
	du0_11 = [2.0]
	tspan_11 = (-1, 13)
	prob_11 = SecondOrderODEProblem(ode_11, du0_11, u0_11, tspan_11, []);
end

# ╔═╡ 86a26f4b-2b10-4228-b6e4-f7f893867980
sol_11 = solve(prob_11, DPRKN6(); reltol=1e-12);

# ╔═╡ 1713685d-8b71-49a8-946a-fc64a2031319
plot(sol_11, idx=[2,1];
	label=["du" "u"],
)

# ╔═╡ ee381a04-9899-42b4-aee7-1e791312102f
md"## Example 12 - Unstable
$\begin{align}
ÿ + 3ẏ - |y|y &= 2\\
y(0) &= 0 \\
ẏ(0) &= 0
\end{align}$
"

# ╔═╡ 934924fd-60f0-47f4-92cd-8993da4c9043
begin
	function ode_12(ddu, du, u, p, t)
		@. ddu = abs(u) * u - 3 * du + 2
	end
	u0_12 = [0.0]
	du0_12 = [0.0]
	prob_12 = SecondOrderODEProblem(ode_12, du0_12, u0_12, tspan, []);
end

# ╔═╡ ba1e3724-ff59-4e80-941d-3af5c08744a3
sol_12 = solve(prob_12, DPRKN6())

# ╔═╡ a235fadf-7f3b-40bc-b3f0-9e566385ca45
md" ## Example 13- Parameters
$\begin{align}
	Aÿ + Bẏ + Cy &= D\\
	y(0) &= a \\
	ẏ(0) &= b
\end{align}$
"

# ╔═╡ f52d33d9-9cb6-4242-bbf9-c19c2431fbe4
md"""
`a = ` $(@bind a NumberField(0:0.1:10, default=1.0))

`b = ` $(@bind b NumberField(0:0.1:10, default=0.0))

`A = ` $(@bind A NumberField(0:0.1:10, default=1.3))

`B = ` $(@bind B NumberField(0:0.5:10, default=1.0))

`C = ` $(@bind C NumberField(0:0.5:10, default=2.0))

`D = ` $(@bind D NumberField(0:0.5:10, default=3.0))
"""

# ╔═╡ 2fd7c463-0cd5-4c73-9482-44705b9abdca
begin
	function ode_13(ddu, du, u, p, t)
		A, B, C, D = p
		@. ddu = (D - B * du - C * u) / A 
	end
	params = [A, B, C, D]
	u0_13 = [a]
	du0_13 = [b]
	prob_13 = SecondOrderODEProblem(ode_13, du0_13, u0_13, tspan, params);
end

# ╔═╡ f88a6192-6a39-40e9-8f16-2b3a62dd9a2a
sol_13 = solve(prob_13, DPRKN6(); reltol=1e-12);

# ╔═╡ 346cc205-1ded-4823-8ff2-8d412a8598b4
plot(sol_13.t, sol_13[2,:];
	label="u",
	title="a=$a, b=$b, A=$A, B=$B, C=$C,$D"
)

# ╔═╡ cf6eb8a6-f751-4858-9f8a-e0d7b20f6dee
md"## Example 15 - Laplace
$\begin{align}
x(t) &= sin(2t)\\
Ł\{x(t)\} &= \frac{2}{(s^2 + 4)}
\end{align}$

!!! warning 
	We're using `SymPy` here because appearently no one cares about analytcial laplace transforms in the Julia ecosystem.
"

# ╔═╡ b9232698-67b6-41f6-9116-ba9d8d2dc45f
begin
	@vars s tₛ postive=true
	f_1 = sin(2tₛ)
	F_1 = laplace_transform(f_1, tₛ, s; noconds=True)
end

# ╔═╡ f7f2c9e0-1bb5-4534-9626-ae94ba5d424f
md"## Example 17 - Laplace
$\begin{align}
y(x) &= ax^3 + b\\
Ł\{y(t)\} &= \frac{6a}{s^4} + \frac{b}{s}
\end{align}$
"

# ╔═╡ 4cd7c0cb-584e-4b84-8253-c102cdec02dc
begin
	@vars xₛ aₛ bₛ
	f_2 = aₛ * xₛ^3 + bₛ
	F_2 = laplace_transform(f_2, xₛ, s; noconds=True)
end

# ╔═╡ 02b69c62-e7e7-4498-bbd4-74bef6f11552
md"## Example 18 - Solve ODE with Laplace
$\begin{align}
ẏ + 2y &= 0 \\
y(0) &= 0.5
\end{align}$
"

# ╔═╡ 3fb87ccb-a9c1-4e52-9d7e-849f86dbd39a
begin
	@vars yₛ
	
	y_18 = SymFunction("y")(tₛ)
	ode_18 = Eq(diff(y_18, tₛ) + 2 * y_18, 0)
	lhs_18 = laplace_transform(ode_18.lhs, tₛ, s)[1]
	rhs_18 = laplace_transform(ode_18.rhs, tₛ, s)[1]
	ODE_18 = Eq(lhs_18, rhs_18)
	Y_18 = solve(ODE_18, laplace_transform(y_18, tₛ, s))[1][1]
	# set initial conditions: we can assume heaviside is 1; y(0) 0.5 > 0
	Y_18 = Y_18.subs(y_18.subs(tₛ, 0), 0.5)
	y_18 = inverse_laplace_transform(Y_18, s, tₛ).subs(Heaviside(tₛ), 1)
end

# ╔═╡ 77136f81-ea83-467c-ba91-e1f454ed7a48
md"Something is wrong with the book plot? [wolfram plot](https://www.wolframalpha.com/input?i2d=true&i=0.5*exp%5C%2840%29-2.0*t%5C%2841%29)"

# ╔═╡ 3fde40b2-b83b-4b67-9fb1-0ad469529696
md"## Example 19 - Solve second order ODE with Laplace
$\begin{align}
ÿ + ẏ &= sin(t) \\
y(0) &= 1 \\
ẏ(0) &= 2
\end{align}$
"

# ╔═╡ a479b28d-ef82-413c-a0a6-81bae6001780
begin
	y_19 = SymFunction("y")(tₛ)
	ode_19 = Eq(diff(y_19, tₛ, 2) + diff(y_19, tₛ), sin(tₛ))
	lhs = laplace_transform(ode_19.lhs, tₛ, s)[1]
	rhs = laplace_transform(ode_19.rhs, tₛ, s)[1]
	ODE_19 = Eq(lhs, rhs)
	Y_19 = solve(ODE_19, laplace_transform(y_19, tₛ, s)[1])[1]
	Y_19 = Y_19.subs(y_19.subs(tₛ, 0), 1).subs(diff(y_19, tₛ), 2)
	Y_19 = simplify(Y_19)
	y_19 = inverse_laplace_transform(Y_19, s, tₛ).subs(Heaviside(tₛ), 1)
end

# ╔═╡ 5e309c0d-e33c-42c0-b9fc-ae0f90948119
plot(y_19, 0, 14;
	title="ÿ + ẏ = sin(t)",
	xlabel="t",
	ylabel="y(t)",
	legend=false
)

# ╔═╡ Cell order:
# ╠═32cccc2d-9eef-43c5-a441-fd112c2272e1
# ╠═1bee3436-edd1-11ed-2fbd-873b2824b63d
# ╠═98ca666d-fca1-4e8f-950d-c66a8cb81db3
# ╟─f62e7973-f1ef-47bb-b29c-b838dec186a6
# ╠═5e48716c-f25f-42f0-a88f-3e32add73f72
# ╠═ab98d373-8dbe-43ff-8457-0ab7cf970248
# ╠═5c9a55da-7b48-45c8-9f9c-d6b4d79d237d
# ╠═c6dbbfd0-d46e-4f60-9287-7ce90a6a1469
# ╟─8327be9e-04fd-4cf3-bfab-f326e285474c
# ╠═8db1e5f5-2591-4c69-a5e0-f9521dce75ac
# ╠═4b201856-e4af-4cd0-bb57-1af254dc26f1
# ╠═f799db33-a4bf-4ec8-b12d-47a2005b6418
# ╟─360d55ea-5da7-47bf-99bc-2dc44ad9b0af
# ╠═a434b1fa-a4c9-4d36-baef-5d34d9abc6ac
# ╠═28b2bb6a-d138-45a3-84d1-e7d5c7d0d712
# ╠═a207cab4-cd88-41c2-8de6-5bdb5860bd9e
# ╟─ad0d1102-a7aa-4942-8b45-76e695662663
# ╠═cbb54be0-d5bc-4755-88b7-ef7e299777be
# ╠═86948555-aa8e-43e6-91d2-74268397acc5
# ╠═133981dd-d36e-41de-a07c-0dd2ee1634cc
# ╟─3d425689-f115-46bc-bea1-903bb64a2100
# ╠═7d4f90cc-0010-4e18-9b90-1b7b86cf0c94
# ╠═86a26f4b-2b10-4228-b6e4-f7f893867980
# ╠═1713685d-8b71-49a8-946a-fc64a2031319
# ╠═ee381a04-9899-42b4-aee7-1e791312102f
# ╠═934924fd-60f0-47f4-92cd-8993da4c9043
# ╠═ba1e3724-ff59-4e80-941d-3af5c08744a3
# ╟─a235fadf-7f3b-40bc-b3f0-9e566385ca45
# ╟─f52d33d9-9cb6-4242-bbf9-c19c2431fbe4
# ╠═2fd7c463-0cd5-4c73-9482-44705b9abdca
# ╠═f88a6192-6a39-40e9-8f16-2b3a62dd9a2a
# ╠═346cc205-1ded-4823-8ff2-8d412a8598b4
# ╟─cf6eb8a6-f751-4858-9f8a-e0d7b20f6dee
# ╠═28a9d730-4017-4786-bc60-8fb4981b4ed1
# ╠═b9232698-67b6-41f6-9116-ba9d8d2dc45f
# ╟─f7f2c9e0-1bb5-4534-9626-ae94ba5d424f
# ╠═4cd7c0cb-584e-4b84-8253-c102cdec02dc
# ╟─02b69c62-e7e7-4498-bbd4-74bef6f11552
# ╠═3fb87ccb-a9c1-4e52-9d7e-849f86dbd39a
# ╟─77136f81-ea83-467c-ba91-e1f454ed7a48
# ╟─3fde40b2-b83b-4b67-9fb1-0ad469529696
# ╠═a479b28d-ef82-413c-a0a6-81bae6001780
# ╠═5e309c0d-e33c-42c0-b9fc-ae0f90948119

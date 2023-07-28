### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 61fbae56-aae4-4070-9634-7c9c980f73e4
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate("..")
end

# ╔═╡ 5560fc2e-b983-4112-9bf9-f83fd9744f25
using Plots, DifferentialEquations, LinearAlgebra, ModelingToolkit

# ╔═╡ 4304e832-0fa3-11ee-0676-672cf9383cb0
md"# Implicit ODEs | Differential Algebraic Equations (DAE)"

# ╔═╡ c5ed3b1a-70cc-4e89-8a7c-2ddc971ae44a
md"## Example 1

$t^2ẏ + tẏ^2 + 2cos(t) = {1 \over y}$
$y(0) = 0.5$
"

# ╔═╡ 8cffc1f8-9405-49d3-ad0b-e7d851f132c5
begin
	function dae_1(residual, dy, y, p, t)
		@. residual = (1 / y) - t^2 * dy - t * dy^2 - 2cos(t)
	end

	tspan_1 = (0, 20)
	dy0_1 = [0.0]
	y0_1 = [0.5]
	prob_1 = DAEProblem(dae_1, dy0_1, y0_1, tspan_1; differential_vars=[true])
	sol_1 = solve(prob_1)
end;

# ╔═╡ 76188789-645e-434c-916d-6d5ca8fb3f96
plot(sol_1)

# ╔═╡ 598ad6e8-1741-4ee4-8158-56f9c7aeb750
md"## Example 2
$ÿ + y = 0$
$y(0) = 0, \ ẏ(0) =1$
"

# ╔═╡ 30e1d81b-e97e-4db5-a6ed-ceafc2d87cbf
begin
	function dae_2(residual, dy, y, p, t)
		residual[1] = dy[1] - y[2]
		residual[2] = dy[2] + y[1]
	end

	tspan_2 = (0, 7)
	dy0_2 = [0.0, 1.0]
	y0_2 = [1.0, 0.0]
	prob_2 = DAEProblem(dae_2, dy0_2, y0_2, tspan_2; differential_vars=[true, true])
	sol_2 = solve(prob_2)
end;

# ╔═╡ a35fdbf8-4255-46b2-a7d0-2dc74f803b33
plot(sol_2)

# ╔═╡ 4c1c71c2-3706-48c2-a815-07becba3fe12
md"## Example 4
$e^{ẍ} + ẍ + x = 0$
$x(0) = 1, \ ẋ(0) = 0$

let

$x_1 = x, \ ẋ_1 = x_2$

then 

$\begin{cases}
	ẋ_1 - x_2 = 0 \\
	e^{ẋ_2} + ẋ_2 + x_1 = 0
\end{cases}$
"

# ╔═╡ 1d2d9eea-c1bc-46f3-961b-5287177cee33
begin
	function dae_4(residual, dy, y, p, t)
		residual[1] = dy[1] - y[2]
		residual[2] = exp(dy[2]) + dy[2] + y[1]
	end

	tspan_4 = (0, 13)
	dy0_4 = [0.0, 0.0]
	y0_4 = [1.0, 0.0]
	prob_4 = DAEProblem(dae_4, dy0_4, y0_4, tspan_4; differential_vars=[true, true])
	sol_4 = solve(prob_4)
end;

# ╔═╡ d3f73626-5796-4d73-8aeb-3e8ed089d3fd
plot(sol_4)

# ╔═╡ 595097e9-be47-49e7-97ba-69c682198a23
md"## Example 6
$e^{0.01t} + ẏe^{y} = 2$
$y(0) = 2$
"

# ╔═╡ f094039c-52b1-4aff-82fa-da7230bdf494
begin
	function dae_6(residual, dy, y, p, t)
		@. residual = exp(0.01t) + dy * exp(y) - 2
	end

	tspan_6 = (0, 129+10e-2)
	dy0_6 = [0.5]
	y0_6 = [2.0]
	prob_6 = DAEProblem(dae_6, dy0_6, y0_6, tspan_6; differential_vars=[true])
	sol_6 = solve(prob_6)
end;

# ╔═╡ 9836a7e2-6a2c-4d1b-b1b7-5589504757aa
plot(sol_6)

# ╔═╡ Cell order:
# ╟─4304e832-0fa3-11ee-0676-672cf9383cb0
# ╠═61fbae56-aae4-4070-9634-7c9c980f73e4
# ╠═5560fc2e-b983-4112-9bf9-f83fd9744f25
# ╟─c5ed3b1a-70cc-4e89-8a7c-2ddc971ae44a
# ╠═8cffc1f8-9405-49d3-ad0b-e7d851f132c5
# ╠═76188789-645e-434c-916d-6d5ca8fb3f96
# ╟─598ad6e8-1741-4ee4-8158-56f9c7aeb750
# ╠═30e1d81b-e97e-4db5-a6ed-ceafc2d87cbf
# ╟─a35fdbf8-4255-46b2-a7d0-2dc74f803b33
# ╟─4c1c71c2-3706-48c2-a815-07becba3fe12
# ╠═1d2d9eea-c1bc-46f3-961b-5287177cee33
# ╟─d3f73626-5796-4d73-8aeb-3e8ed089d3fd
# ╠═595097e9-be47-49e7-97ba-69c682198a23
# ╠═f094039c-52b1-4aff-82fa-da7230bdf494
# ╟─9836a7e2-6a2c-4d1b-b1b7-5589504757aa

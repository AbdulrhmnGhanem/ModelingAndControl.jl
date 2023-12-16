### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 287a6854-0c8b-11ee-189b-f505088c3bad
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ f1e1bf1c-c5cd-48d9-8536-42097dada190
using Plots, DifferentialEquations, BenchmarkTools

# ╔═╡ 50c4c999-8c36-41a2-a631-aa238822cb14
md"# Stiff ODEs


Stiff ODEs (ordinary differential equations) are a type of differential equation for which certain numerical methods for solving the equation are numerically unstable. The phenomenon of stiffness arises when the differential system has widely varying time scales, such as when there are very fast and very slow reactions in a chemical reaction mechanism .

Stiff ODEs pose a challenge for numerical solvers because they require very small time steps to be stable, which can make the computation time-consuming and impractical.

To solve stiff ODEs, specialized numerical methods called stiff solvers are used. These solvers typically do more work per step but can take much larger steps, making them more efficient for stiff problems.

In summary, stiff ODEs are differential equations that require specialized numerical methods to solve due to widely varying time scales in the differential system.
"

# ╔═╡ dd37a837-c889-41b4-a4bd-c831d6750b73
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ 4a364f0f-b6f8-450a-9419-300e3c356118
function tspan_to_range(tspan, h; collect=false)
	r = tspan[1]:h:tspan[2]
	if collect
		return Base.collect(r)
	end
	return r
end

# ╔═╡ 0575d574-c655-4c24-aafd-ad4555310477
# ╠═╡ show_logs = false
begin
	first_order = ingredients("ch_02.jl")
	second_order = ingredients( "ch_03.jl")
end;

# ╔═╡ f3dc814b-caae-49a1-be98-6ff76cf7b470
md"### Example 1

${dy \over dt} + 100 y = 0$
$y(0) = 1$

We are going to solve it analytcially, with forward euler, and with Runge-Kutta.
"

# ╔═╡ d139cc0e-617c-472d-816b-35df4e868953
begin
	df_1(u, t) = -100t

	f0_1 = 1
	h_1_a =  0.015
	h_1_b = h_1_a / 2
	tspan_1 = (0, 0.2)

	analytical_f(t) = exp(-100t)
	sol_analytical_1 = tspan_to_range(tspan_1, h_1_b/2, collect=true), analytical_f.(tspan_to_range(tspan_1, h_1_b/2))

	sol_euler_1_a = first_order.euler(df_1, h_1_a, tspan_1, f0_1)
	sol_euler_1_b = first_order.euler(df_1, h_1_b, tspan_1, f0_1)
	sol_rk_1 = first_order.runge_kutta(df_1, h_1_b, tspan_1, f0_1)
end;

# ╔═╡ 8340e803-d98a-45f0-849b-6a4f7dbbe9a7
begin
	plot(sol_analytical_1;
		title="Analytical vs Euler vs Runge-Kutta Solutions",
		xlabel="t",
		ylabel="u(t)",
		label="Analytical"
	)
	plot!(sol_euler_1_a; label="Euler h=$h_1_a")
	plot!(sol_euler_1_b; label="Euler h=$h_1_b")
	plot!(sol_rk_1, label="Runge-Kutta h=$h_1_b")
end

# ╔═╡ d6915da7-facf-4e8a-b698-b3e282ebe076
md"## Using DifferentialEquations.jl

![Imgur](https://i.imgur.com/GObvC8H.png)"

# ╔═╡ 2ff9fc5a-392a-4964-af2e-d57ed68c9855
md"### Example 2
${du \over dt} + 1000u = 0$
$u(0) = 1$
"

# ╔═╡ 1dcac265-708c-457c-a0d3-3a3b972c57b3
begin
	function df_2(du, u, p, t)
		@. du = -1000u
	end
	f0_2 = [1]
	tspan_2 = (0,0.2)
	
	prob_2 = ODEProblem(df_2, f0_2, tspan_2)
end;

# ╔═╡ afcdb864-d705-4537-b5c3-3c5aaa966815
md"#### Non-stiff sovler"

# ╔═╡ c5fe5b66-016b-426a-b7cc-fe673d3091d6
sol_2_tist = @btime solve(prob_2, Tsit5();
		reltol=1e-3,
		abstol=1e-5,
);

# ╔═╡ 10dc5bce-7a41-4158-9460-408025db4694
md"#### Stiff solver"

# ╔═╡ 635aff41-8939-485c-906e-c28a849ae33c
sol_2_rodas = @btime solve(prob_2, Rodas4();
	reltol=1e-3,
		abstol=1e-5,
);

# ╔═╡ 7622f3e0-4070-4ee6-b14c-f26533be85aa
begin
	plot(sol_2_tist[1:16];
		title="Stiff and non-Stiff sovlers solution",
		xlabel="t",
		ylabel="u(t)",
		label="Tsit5 (non-stiff)",
	)
	plot!(sol_2_rodas[1:16];
		label="Rodas4 (stiff)",
		ls=:dash
	)
end

# ╔═╡ 60b75a0b-785c-450d-ab98-ff8998c73045
md"### Example 3"

# ╔═╡ 3d1022c5-3c7f-49a5-85d6-5907388d69a9
md"### DifferentialEquations.jl outperforms MATLAB

!!! danger
	Althogh the purpose of this example in the book is to illusterate that sometimes euler performs better than the MATLAB ODEs we can't do that; DE.jl solvers are indeed much better than Euler.
"

# ╔═╡ e655b3f1-0918-4b95-86d4-d96e41dfea39
begin
	function df_3(du, u, p, t)
		 du = cos(1001t)
	end
	
	f0_3 = [1.001]
	tspan_3 = (0, 0.75)
	prob_3 = ODEProblem(df_3, f0_3, tspan_3);
end;

# ╔═╡ 20f80f83-951d-4480-a9b8-dfc63e72ddbd
begin
	sol_3_rosenbrock = solve(prob_3, Rosenbrock23();
	)
	q = 0
	sol_3_euler = first_order.euler((u, t) -> df_3(q, u, [], t),  h_1_b, tspan_3, f0_3[1])
	
	sol_3_dp5 = solve(prob_3, DP5();
			reltol=1e-3,
			abstol=1e-5,
	)

	sol_3_vcabm = solve(prob_3, VCABM();
			reltol=1e-3,
			abstol=1e-5,
	)

	sol_3_qndf = solve(prob_3, QNDF();
			reltol=1e-9,
			abstol=1e-9,
	)
end;

# ╔═╡ a74c72ea-143d-4dbd-a004-b4c59988bc23
begin
	plot(sol_3_euler;
		title="Euler vs DifferentialEquations.jl sovlers solution",
		xlabel="t",
		ylabel="u(t)",
		label="Euler",
		xlims=(0,0.75),
		ylims=(0.9, 1.1)
	)
	plot!(sol_3_dp5, xlims=(0,0.75); label="DP5", ls=:dash)
	plot!(sol_3_vcabm, xlims=(0,0.75), label="VCABM", ls=:dashdot)
	plot!(sol_3_qndf, xlims=(0,0.75), label="QNDF", ls=:dot)
	
end

# ╔═╡ 239e7a5f-c883-4bda-bb26-7e025c4c6134
md"### Example 4

$ÿ = 100 ẏ - 10.9y$
$y(0) = 1, \ ẏ(0) = 0$
"

# ╔═╡ 53ca1315-53f4-412d-a7c6-122e5453ca80
function ddf_4(ddu, du, u, p, t)
	ddu = 100du - 10.9u
end

# ╔═╡ 7b5624e6-3d9f-48d1-bb09-af9a3a9f530c
begin
	tspan_4 = (0, 100)
	f0_4 = [1.0]
	df0_4 = [0.0]
	prob_4 = SecondOrderODEProblem(ddf_4, df0_4, f0_4, tspan_4)
end

# ╔═╡ 735b3c21-68e0-4194-8453-f29e5f34fbcf
sol_4_tist = @btime solve(prob_4, Tsit5();
	reltol=1e-6,
	abstol=1e-8,
);

# ╔═╡ 83cf1071-63f2-4ebf-8131-eb4fe3696010
md"### Example 6"

# ╔═╡ 08be0d20-3a65-4398-a71e-4cb3b93f904f
function ddf_6_1(ddu, du, u, p, t)
	@. ddu = (10^4) * abs(sin(333t)) - 100du - (10^4)u
end

# ╔═╡ 77f8c161-3eea-4162-8b30-8af08c7da4b5
begin
	tspan_6 = (0, 0.1)
	f0_6 = [0.0]
	df0_6 = [0.0]
	prob_6 = SecondOrderODEProblem(ddf_6_1, df0_6, f0_6, tspan_6)
end;

# ╔═╡ da96a526-8456-4127-863e-5c4e4a147bfd
sol_6_rodas = @btime solve(prob_6, Rodas4();
	reltol=1e-6,
	abstol=1e-8,
);

# ╔═╡ 4a913865-b84d-4d2b-9168-ddd0d93f04ab
begin
	input(t) = abs(sin(333t))
	ddf_6_2 = second_order.@second_order (du, u, t) -> 
				@. (10^4) * input(t) - 100du - (10^4)u

	h = 10e-5
	sol_6_euler = @btime second_order.euler_2nd_order(ddf_6_2, h, tspan_6, f0_6[1], df0_6[1])
	inputs = tspan_to_range(tspan_6, h; collect=true), input.(tspan_to_range(tspan_6, h))
end;

# ╔═╡ d122d886-9506-4847-b6d4-31869cc1da1c
begin
	plot(sol_6_euler;
		title="Electric rectifier circuit voltage simulation",
		xlabel="t",
		ylabel="voltage",
		label="output"
	)
	plot!(inputs;
		xlims=(0, 0.1),
		label="input",
	)
end

# ╔═╡ cee124e3-c3f9-4971-8499-bd6e7ab9976e
md"
The euler approximation isn't only faster than other ODE solvers, it is much faster than MATLAB implementation:
julia in mciro seconds while MATLAB in milliseconds."

# ╔═╡ Cell order:
# ╟─50c4c999-8c36-41a2-a631-aa238822cb14
# ╠═287a6854-0c8b-11ee-189b-f505088c3bad
# ╠═f1e1bf1c-c5cd-48d9-8536-42097dada190
# ╟─dd37a837-c889-41b4-a4bd-c831d6750b73
# ╟─4a364f0f-b6f8-450a-9419-300e3c356118
# ╠═0575d574-c655-4c24-aafd-ad4555310477
# ╟─f3dc814b-caae-49a1-be98-6ff76cf7b470
# ╠═d139cc0e-617c-472d-816b-35df4e868953
# ╟─8340e803-d98a-45f0-849b-6a4f7dbbe9a7
# ╟─d6915da7-facf-4e8a-b698-b3e282ebe076
# ╟─2ff9fc5a-392a-4964-af2e-d57ed68c9855
# ╠═1dcac265-708c-457c-a0d3-3a3b972c57b3
# ╟─afcdb864-d705-4537-b5c3-3c5aaa966815
# ╠═c5fe5b66-016b-426a-b7cc-fe673d3091d6
# ╟─10dc5bce-7a41-4158-9460-408025db4694
# ╠═635aff41-8939-485c-906e-c28a849ae33c
# ╟─7622f3e0-4070-4ee6-b14c-f26533be85aa
# ╟─60b75a0b-785c-450d-ab98-ff8998c73045
# ╟─3d1022c5-3c7f-49a5-85d6-5907388d69a9
# ╠═e655b3f1-0918-4b95-86d4-d96e41dfea39
# ╠═20f80f83-951d-4480-a9b8-dfc63e72ddbd
# ╟─a74c72ea-143d-4dbd-a004-b4c59988bc23
# ╟─239e7a5f-c883-4bda-bb26-7e025c4c6134
# ╠═53ca1315-53f4-412d-a7c6-122e5453ca80
# ╠═7b5624e6-3d9f-48d1-bb09-af9a3a9f530c
# ╠═735b3c21-68e0-4194-8453-f29e5f34fbcf
# ╟─83cf1071-63f2-4ebf-8131-eb4fe3696010
# ╠═08be0d20-3a65-4398-a71e-4cb3b93f904f
# ╠═77f8c161-3eea-4162-8b30-8af08c7da4b5
# ╠═da96a526-8456-4127-863e-5c4e4a147bfd
# ╠═4a913865-b84d-4d2b-9168-ddd0d93f04ab
# ╟─d122d886-9506-4847-b6d4-31869cc1da1c
# ╟─cee124e3-c3f9-4971-8499-bd6e7ab9976e

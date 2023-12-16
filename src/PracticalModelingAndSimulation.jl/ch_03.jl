### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ da3a4744-0ab0-11ee-0066-05339cb401a8
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ 98fd8b9f-436e-4050-a212-e47d0318bc66
using Plots

# ╔═╡ f87e4cef-cd95-48a2-b414-c1ed2456014b
md"# Numerical Methods for Second-Order ODEs"

# ╔═╡ 8da5c63d-50ea-4b4a-9b56-e45b4b497f9e
md"## Analytical solution"

# ╔═╡ e9596408-9aad-4748-89fe-2e97822ca8bb
macro second_order(f)
	"""Modify the functiuon to return both first and second derivatives"""
	quote
		function wrapped_function(du, u, t)
            return [du, $(esc(f))(du, u, t)]
		end
	end
end

# ╔═╡ ccce800e-9d82-4f39-839f-31a93f8f6176
md"## Forward euler"

# ╔═╡ 7ea65c37-4251-4258-bcd1-87d23e2675b6
function euler_2nd_order(f, h, tspan, f0, f1)
	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	y_ = [[f0, f1], Iterators.repeated([0.0, 0.0], steps_num)...]
	for i in range(1, steps_num-1)
		f_ = f(y_[i][2], y_[i][1], tspan[i])
		y_[i+1] = y_[i] + f_ * h
	end
	y_ = y_[1:end-1]
	return tspan, [el[1] for el in y_]
end

# ╔═╡ 5c082b73-cc18-4a64-8858-766a2db03889
begin
	tspan = (0, 10)
	ddf = @second_order (du, u, t) -> 4t - (4/5)du - 2u

	h_0 = 0.2
	f0 = 1
	f1 = 2

	sol_euler_1 = euler_2nd_order(ddf, h_0, tspan, f0, f1)
end;

# ╔═╡ 260254c2-2dd9-4640-a330-677ca7ae6516
begin
	plot(sol_euler_1;
		title="Forward Euler solution",
		xlabel="t",
		ylabel="u(t)",
		label="Euler"
	)
end

# ╔═╡ 038b27e2-9bcf-4a1b-ac53-a3912a19b31e
md"### Euler - Example 2"

# ╔═╡ 176c5d9a-2edd-463e-86c4-c3c03ec0febb
begin
	tspan_2 = (0, 2.5)
	ddf_2 = @second_order (du, u, t) -> -6du - 9u

	f_2_anaytical(u) = exp(-3u) * (2 + 6u)
	sol_analytical_2 = collect(tspan[1]:h_0:tspan[2]), f_2_anaytical.(collect(tspan[1]:h_0:tspan[2]))

	h_2_a = 0.1
	h_2_b = 0.025
	f0_2 = 2.0
	f1_2 = 0.0
	sol_euler_2_a = euler_2nd_order(ddf_2, h_2_a, tspan_2, f0_2, f1_2)
	sol_euler_2_b = euler_2nd_order(ddf_2, h_2_b, tspan_2, f0_2, f1_2)
end;

# ╔═╡ 7fcf9bd6-62c9-4ff7-88ef-a27a0eb0d0ed
begin
	plot(sol_analytical_2;
		title="Analytical vs Forward Euler (with different step size)\n solutions",
		xlabel="t",
		ylabel="u(t)",
		label="Analytical"
	)

	plot!(sol_euler_2_a;
		label="Euler (h=$h_2_a)",
		ls=:dot
	)
	plot!(sol_euler_2_b;
		label="Euler (h=$h_2_b)",
		ls=:dash
	)
end

# ╔═╡ 9ced3fdc-9a9d-49e9-a8fc-e60b8c4ad4c0
md"### Euler - Example 3"

# ╔═╡ 6068434b-15ea-47a3-a8a9-71a0c891b061
begin
	tspan_3 = (0, 12)
	ddf_3 = @second_order (du, u, t) -> -sin(t)

	f_3_anaytical(u) = sin(u) - u
	sol_analytical_3 = collect(tspan_3[1]:h_0:tspan_3[2]), f_3_anaytical.(collect(tspan_3[1]:h_0:tspan_3[2]))

	h_3 = 0.1
	f0_3 = 0.0
	f1_3 = 0.0
	sol_euler_3 = euler_2nd_order(ddf_3, h_3, tspan_3, f0_3, f1_3)
end;

# ╔═╡ 55f7cf6a-b02d-45aa-ba2f-681e55b98df7
begin
	plot(sol_analytical_3;
		title="Analytical vs Forward Euler",
		xlabel="t",
		ylabel="u(t)",
		label="Analytical"
	)

	plot!(sol_euler_3;
		label="Euler",
		ls=:dot
	)
end

# ╔═╡ 88feeded-7381-42f5-8d00-9e10687c7b26
md"### Euler - Example 4"

# ╔═╡ 6bfa972a-9c93-479f-853d-e268e06ea0ce
begin
	tspan_4 = (0, 20)
	ddf_4 = @second_order (du, u, t) -> sin(t) - du

	f_4_anaytical(u) = (1/2) * (-5 * exp(-u) - sin(u) - cos(u) + 8)
	sol_analytical_4 = collect(tspan_4[1]:h_0:tspan_4[2]), f_4_anaytical.(collect(tspan_4[1]:h_0:tspan_4[2]))

	h_4 = 0.1
	f0_4 = 1.0
	f1_4 = 2.0
	sol_euler_4 = euler_2nd_order(ddf_4, h_4, tspan_4, f0_4, f1_4)
end;

# ╔═╡ 8dea96c4-f83a-4899-8461-c42020c0a68d
begin
	plot(sol_analytical_4;
		title="Analytical vs Forward Euler",
		xlabel="t",
		ylabel="u(t)",
		label="Analytical"
	)

	plot!(sol_euler_4;
		label="Euler",
		ls=:dot
	)
end

# ╔═╡ 02b366a8-9afd-4cde-98eb-86f56e7e7217
md"### Euler - Example 5"

# ╔═╡ 205182a7-4f6a-4b80-83a1-a8e08fd81a9f
begin
	tspan_5 = (0, 20)
	ddf_5 = @second_order (du, u, t) -> sin(u*t) - du

	h_5_a = 0.025
	h_5_b = 0.05
	h_5_c = 0.1
	h_5_d = 0.2

	f0_5 = 1.0
	f1_5 = 2.0

	sol_euler_5_a = euler_2nd_order(ddf_5, h_5_a, tspan_5, f0_5, f1_5)
	sol_euler_5_b = euler_2nd_order(ddf_5, h_5_b, tspan_5, f0_5, f1_5)
	sol_euler_5_c = euler_2nd_order(ddf_5, h_5_c, tspan_5, f0_5, f1_5)
	sol_euler_5_d = euler_2nd_order(ddf_5, h_5_d, tspan_5, f0_5, f1_5)
end;

# ╔═╡ 6d3c2452-467c-44bc-84c2-75f09a739f63
begin
	plot(sol_euler_5_a;
		title="Different solution with different step sizes",
		xlabel="t",
		ylabel="u(t)",
		label="h=$h_5_a"
	)
	plot!(sol_euler_5_b;
		label="h=$h_5_b"
	)
	plot!(sol_euler_5_c;
		label="h=$h_5_c"
	)
	plot!(sol_euler_5_d;
		label="h=$h_5_d"
	)
end

# ╔═╡ 1f654a09-b327-4ec7-a493-60d8a10170d9
md"## Runge-Kutta"

# ╔═╡ c05d9b64-ac77-46ca-b82a-34899676d9b0
function runge_kutta_2nd_order(f, h, tspan, f0, f1)
	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	y_ = [[f0, f1], Iterators.repeated([0.0, 0.0], steps_num)...]
	for i in range(1, steps_num - 1)
		k1 = f(y_[i][2], y_[i][1], tspan[i])
		k2 = f(y_[i][2]+ h * k1[2] / 2, y_[i][1] + h * k1[1] / 2, tspan[i]+ h / 2)
		k3 = f(y_[i][2] + h * k2[2] / 2, y_[i][1] + h * k2[1] / 2, tspan[i] + h / 2)
		k4 = f(y_[i][2] + h * k3[2], y_[i][1] + h * k3[1], tspan[i] + h)
		y_[i+1] = y_[i] + h * (k1 + 2k2 + 2k3 + k4)/6
	end
	y_ = y_[1:end-1]
	return tspan, [el[1] for el in y_]
end

# ╔═╡ b79167cb-ea4f-44d8-9d37-3f7976a95194
md"### Runge-Kutta - Example 6"

# ╔═╡ 366175ef-e95e-4c49-ab8a-e8d9c00ca952
begin
	tspan_6 = (0, 2)
	ddf_6 = @second_order (du, u, t) -> 4t - (4/5)du - 2u

	h_6_a = 0.1
	h_6_b = 0.2
	h_6_c = 0.3
	h_6_d = 0.4

	f0_6 = 1.0
	f1_6 = 2.0

	sol_rk_6_a = runge_kutta_2nd_order(ddf_6, h_6_a, tspan_6, f0_6, f1_6)
	sol_rk_6_b = runge_kutta_2nd_order(ddf_6, h_6_b, tspan_6, f0_6, f1_6)
	sol_rk_6_c = runge_kutta_2nd_order(ddf_6, h_6_c, tspan_6, f0_6, f1_6)
	sol_rk_6_d = runge_kutta_2nd_order(ddf_6, h_6_d, tspan_6, f0_6, f1_6)
end;

# ╔═╡ 9143f1d7-76eb-4565-946a-4ee5953bdf60
begin
	plot(sol_rk_6_a;
		title="Different solution with different step sizes",
		xlabel="t",
		ylabel="u(t)",
		label="h=$h_6_a"
	)
	plot!(sol_rk_6_b;
		label="h=$h_6_b",
		ls=:dot
	)
	plot!(sol_rk_6_c;
		label="h=$h_6_c",
		ls=:dash
	)
	plot!(sol_rk_6_d;
		label="h=$h_6_d",
		ls=:dashdot
	)
end

# ╔═╡ 28895519-34b7-4484-a40a-af64a5d3f1f9
md"### Runge-Kutta - Example 7"

# ╔═╡ 8c08d1b1-52af-4f6b-ba74-4ac9f2a3d979
begin
	tspan_7 = (0, 3)
	ddf_7 = @second_order (du, u, t) -> -6du - 9u

	f_7_anaytical(u) = exp(-3u) * (6u + 2)
	sol_analytical_7 = collect(tspan_7[1]:h_0:tspan_7[2]), f_7_anaytical.(collect(tspan_7[1]:h_0:tspan_7[2]))

	h_7 = 0.2
	f0_7 = 2.0
	f1_7 = 0.0

	sol_rk_7 = runge_kutta_2nd_order(ddf_7, h_7, tspan_7, f0_7, f1_7)
	sol_euler_7 = euler_2nd_order(ddf_7, h_7, tspan_7, f0_7, f1_7)
end;

# ╔═╡ 6ef0c89d-8f93-4385-abfe-d2d5197e0dc3
begin
	plot(sol_rk_7;
		title="Runge-Kutta vs Euler vs Analytical solution",
		xlabel="t",
		ylabel="u(t)",
		label="Runge-Kutta"
	)

	plot!(sol_euler_7;
		ls=:dot,
		label="Euler"
	)

	plot!(sol_analytical_7;
		label="Analytical",
		ls=:dash
	)
end

# ╔═╡ 00c858d6-722f-4ef9-91dd-9073fc93fb14
md"### Runge-Kutta - Example 8"

# ╔═╡ 5c48ab36-0148-4d6e-afbe-65306e6f58a7
begin
	tspan_8 = (0, 15)
	ddf_8 = @second_order (du, u, t) -> -sin(t)

	h_8 = 0.3
	f0_8 = 0.0
	f1_8 = 0.0

	sol_rk_8 = runge_kutta_2nd_order(ddf_8, h_8, tspan_8, f0_8, f1_8)
	sol_euler_8 = euler_2nd_order(ddf_8, h_8, tspan_8, f0_8, f1_8)
end;

# ╔═╡ 3fc156b7-47d3-4bd0-acc4-82e5ebab5705
begin
	plot(sol_rk_8;
		title="Runge-Kutta vs Euler solution",
		xlabel="t",
		ylabel="u(t)",
		label="Runge-Kutta"
	)

	plot!(sol_euler_8;
		ls=:dot,
		label="Euler"
	)
end

# ╔═╡ c953a1bd-0c46-4e2a-a664-ea58ffcab1de
begin
	tspan_9 = (0, 15)
	ddf_9 = @second_order (du, u, t) -> sin(t) - du

	f_9_anaytical(u) = (1/2) * (-5 * exp(-u) - sin(u) - cos(u) + 8)
	sol_analytical_9 = collect(tspan_9[1]:0.0001:tspan_9[2]), f_9_anaytical.(collect(tspan_9[1]:0.0001:tspan_9[2]))

	h_9 = 0.25
	f0_9 = 1.0
	f1_9 = 2.0

	sol_rk_9 = runge_kutta_2nd_order(ddf_9, h_9, tspan_9, f0_9, f1_9)
	sol_euler_9 = euler_2nd_order(ddf_9, h_9, tspan_9, f0_9, f1_9)
end;

# ╔═╡ 13a787f3-976d-4bec-b041-47d2edf4d66f
begin
	plot(sol_rk_9;
		title="Runge-Kutta vs Euler vs Analytical solution",
		xlabel="t",
		ylabel="u(t)",
		label="Runge-Kutta"
	)

	plot!(sol_analytical_9;
		ls=:dash,
		label="Analytcial"
	)

	plot!(sol_euler_9;
		ls=:dot,
		label="Euler"
	)
end

# ╔═╡ 51b45585-82ac-4c5f-be87-be42ff69d311
md"## Adams-Multon"

# ╔═╡ 3c5bd97f-80cb-40a8-afb4-3b8dbb85d010
function adams_multon_2nd_order(f, h, tspan, f0, f1; steps=4)
	allowed_steps = (1, 2, 3, 4)
	if !(steps in allowed_steps)
		throw(error("Invalid step size for Adams-Multon! Allowed steps sizes are $allowed_steps."))
	end

	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	y_ = [[f0, f1], Iterators.repeated([0.0, 0.0], steps_num)...]

	lookahead = steps - 1
	for i in range(1, steps_num - lookahead)
		k1() = f(y_[i][2],y_[i][1], tspan[i])
		k2() = f(y_[i+1][2],y_[i+1][1], tspan[i+1])
		k3() = f(y_[i+2][2],y_[i+2][1], tspan[i+2])
		k4() = f(y_[i+3][2],y_[i+3][1], tspan[i+3])

		y_[i+1] = y_[i] + h * k1()
		if steps == 1
			continue
		end

		y_[i+1] = y_[i] + (h/2) * (k1() + k2())
		if steps == 2
			continue
		end

		y_[i+2] = y_[i+1] + (3h/2) * k2() - (h/2) * k1()
		y_[i+2] = y_[i+1] + h * ((5/12) * k3() + (2/3) * k2() - (1/2) * k1())

		if steps == 3
			continue
		end

		y_[i+3] = y_[i+2] + h * ((23/12) * k3() - (4/3) * k2() + (5/12) * k1())
		y_[i+3] = y_[i+2] + h * ((3/8) * k4() + (19/24) * k3() - (5/24) * k2() + (1/24) * k1())
	end

	y_ = y_[1:end-1]
	return tspan, [el[1] for el in y_]
end

# ╔═╡ cdd7d81c-d0e8-4310-a763-231319ed2eca
md"### Adams-Multon - Example 12"

# ╔═╡ 30592fc4-0779-407c-9b7f-b8737cd028dc
begin
	tspan_12 = (0, 15)
	ddf_12 = @second_order (du, u, t) -> sin(u*t) - du

	h_12 = 0.1
	f0_12 = 1.0
	f1_12 = 2.0

	sol_am_12_1 = adams_multon_2nd_order(ddf_12, h_12, tspan_12, f0_12, f1_12; steps=1);
	sol_am_12_2 = adams_multon_2nd_order(ddf_12, h_12, tspan_12, f0_12, f1_12; steps=2);
	sol_am_12_3 = adams_multon_2nd_order(ddf_12, h_12, tspan_12, f0_12, f1_12;steps=3);
	sol_am_12_4 = adams_multon_2nd_order(ddf_12, h_12, tspan_12, f0_12, f1_12; steps=4);
end;

# ╔═╡ e35a0804-f4f0-424f-a561-03e7cd0a7337
begin
	plot(sol_am_12_1;
		title="Adams-Multon solution with different step numbers\nh=$h_12",
		xlabel="t",
		ylabel="u(t)",
		label="steps=1",
		ls=:dash
	)

	plot!(sol_am_12_2;
		label="steps=2",
		ls=:dot
	)

	plot!(sol_am_12_3;
		label="steps=3",
		ls=:dashdot
	)

	plot!(sol_am_12_4;
		label="steps=4"
	)
end

# ╔═╡ Cell order:
# ╟─f87e4cef-cd95-48a2-b414-c1ed2456014b
# ╠═da3a4744-0ab0-11ee-0066-05339cb401a8
# ╠═98fd8b9f-436e-4050-a212-e47d0318bc66
# ╟─8da5c63d-50ea-4b4a-9b56-e45b4b497f9e
# ╟─e9596408-9aad-4748-89fe-2e97822ca8bb
# ╟─ccce800e-9d82-4f39-839f-31a93f8f6176
# ╠═7ea65c37-4251-4258-bcd1-87d23e2675b6
# ╠═5c082b73-cc18-4a64-8858-766a2db03889
# ╠═260254c2-2dd9-4640-a330-677ca7ae6516
# ╟─038b27e2-9bcf-4a1b-ac53-a3912a19b31e
# ╠═176c5d9a-2edd-463e-86c4-c3c03ec0febb
# ╟─7fcf9bd6-62c9-4ff7-88ef-a27a0eb0d0ed
# ╟─9ced3fdc-9a9d-49e9-a8fc-e60b8c4ad4c0
# ╠═6068434b-15ea-47a3-a8a9-71a0c891b061
# ╟─55f7cf6a-b02d-45aa-ba2f-681e55b98df7
# ╟─88feeded-7381-42f5-8d00-9e10687c7b26
# ╠═6bfa972a-9c93-479f-853d-e268e06ea0ce
# ╠═8dea96c4-f83a-4899-8461-c42020c0a68d
# ╟─02b366a8-9afd-4cde-98eb-86f56e7e7217
# ╠═205182a7-4f6a-4b80-83a1-a8e08fd81a9f
# ╟─6d3c2452-467c-44bc-84c2-75f09a739f63
# ╟─1f654a09-b327-4ec7-a493-60d8a10170d9
# ╠═c05d9b64-ac77-46ca-b82a-34899676d9b0
# ╟─b79167cb-ea4f-44d8-9d37-3f7976a95194
# ╠═366175ef-e95e-4c49-ab8a-e8d9c00ca952
# ╠═9143f1d7-76eb-4565-946a-4ee5953bdf60
# ╟─28895519-34b7-4484-a40a-af64a5d3f1f9
# ╠═8c08d1b1-52af-4f6b-ba74-4ac9f2a3d979
# ╠═6ef0c89d-8f93-4385-abfe-d2d5197e0dc3
# ╟─00c858d6-722f-4ef9-91dd-9073fc93fb14
# ╠═5c48ab36-0148-4d6e-afbe-65306e6f58a7
# ╠═3fc156b7-47d3-4bd0-acc4-82e5ebab5705
# ╠═c953a1bd-0c46-4e2a-a664-ea58ffcab1de
# ╠═13a787f3-976d-4bec-b041-47d2edf4d66f
# ╟─51b45585-82ac-4c5f-be87-be42ff69d311
# ╠═3c5bd97f-80cb-40a8-afb4-3b8dbb85d010
# ╟─cdd7d81c-d0e8-4310-a763-231319ed2eca
# ╠═30592fc4-0779-407c-9b7f-b8737cd028dc
# ╟─e35a0804-f4f0-424f-a561-03e7cd0a7337

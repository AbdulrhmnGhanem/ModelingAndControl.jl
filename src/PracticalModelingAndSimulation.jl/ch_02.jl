### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 3a12645a-ff8f-4540-a1b2-1f75f8a4e3f9
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ f82bcf92-9501-4aac-8399-2659fdbd9593
using Plots, DifferentialEquations

# ╔═╡ 4a57db60-efcf-11ed-34d4-690120ffb2c4
md"# Numerical Methods for First-Order ODEs"

# ╔═╡ 641fad95-8abb-4818-a297-a4244e261283
begin
	f(t, y) = 3 - 2y + exp(-t)
	tspan = (1, 10)
	f0_0 = 1
	h_0 = 0.1
end;

# ╔═╡ 11cee654-fe7b-4376-9a51-5d6ae99ef447
md"## Euler Method"

# ╔═╡ 3e883f70-3061-4588-a3fd-86c651c294c0
md"
### forward euler
$$\begin{cases}
    y_i = y_{i-1} + hk \\
    k = f(t_{i-1}, y_{i-1}) \\
    t_i = t_{i-1} + h|
\end{cases}$$"

# ╔═╡ 015dcf60-51dc-4d6e-ac07-0f2c847067e7
function euler(f, h, tspan, f0)
	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	f_ = zeros(steps_num)
	y_ = [f0, zeros(steps_num)...]
	for i in range(1, steps_num)
		f_[i] = f(tspan[i], y_[i])
		y_[i+1] = y_[i] + f_[i] * h
	end
	return tspan, y_[1:end-1]
end;

# ╔═╡ 02d916f4-1fa3-4fbf-b118-ef9b30cf40d6
sol_euler_1 = euler(f, h_0, tspan, f0_0);

# ╔═╡ 1131ac8d-42a3-4777-a538-6b89a805e13f
md"#### forward euler vs auto diff"

# ╔═╡ d8cc309f-1035-4b7b-9215-af465eb8fd9b
begin
	function ode_1(du, u, p, t)
		@. du = 3 - 2u + exp(-t)
	end
	prob_1 = ODEProblem(ode_1, [f0_0], tspan)
	sol_ode_1 = solve(prob_1)
end;

# ╔═╡ ecff5b0c-5474-4094-ab0b-36aac496bca6
begin
	plot(sol_ode_1;
		title="Euler vs ODE solution",
		label="ODE"
	)
	plot!(sol_euler_1;
		ls=:dash,
		label="Euler"
	)
end

# ╔═╡ 6e685a10-a662-4a9d-b9e0-102a20919ea2
md"### Improved Euler (Heun's Method)

$$y_{i+1} = y_i + (f(t_i, y_i)+f(t_i+h, y_i+hf_i)){h \over 2}$$
"

# ╔═╡ 478f44f4-6470-4abf-9258-da6908b02a99
function improved_euler(f, h, tspan, f0)
	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	y_ = [f0, zeros(steps_num)...]
	for i in range(1, steps_num)
		k1 = f(tspan[i], y_[i])
		k2 = f(tspan[i] + h, y_[i] + h * k1)
		y_[i+1] = y_[i] + (k1 + k2) * (h / 2)
	end
	return tspan, y_[1:end-1]
end

# ╔═╡ 1cfe11dc-e7e8-4fe6-9e5b-e3ddb2f09964
sol_euler_2 = improved_euler(f, h_0, tspan, f0_0);

# ╔═╡ 67c546d1-6d5e-4a3b-9730-28f5d00a2c6e
md"### Euler vs Improved Euler vs auto diff"

# ╔═╡ ad9fd0c7-4172-4d40-93bd-3e83c1689b54
begin
	plot(sol_ode_1;
		title="Euler vs Imporved Euler vs solution",
		label="ODE"
	)
	plot!(sol_euler_1;
		ls=:dashdotdot,
		label="Euler"
	)
	plot!(sol_euler_2;
		ls=:dash,
		label="Improved Euler"
	)
end

# ╔═╡ 0664d3fa-89f2-4338-8240-e0aa9c0866d8
md"### Bacward Euler
$\begin{cases}
    t_i = t_{i-1} + h \\
    k = f(t_{i}, y_{i}) \\
	y^*_{i} = y_i + hk \\
	y_i = y_{i-1} + hf(t_i, y^*)
\end{cases}$
"

# ╔═╡ 4c4acf15-99b2-41ac-b9b1-31557c6393a2
function backward_euler(f, h, tspan, f0)
	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	y_ = [f0, zeros(steps_num)...]
	for i in range(1, steps_num)
		k = f(tspan[i], y_[i])
		y_star = y_[i] + h * k
		y_[i+1] = y_[i] + h * f(tspan[i], y_star)
	end
	return tspan, y_[1:end-1]
end

# ╔═╡ 2d6cc393-2c31-4699-bb2f-14d5014b42ca
sol_euler_3 = backward_euler(f, h_0, tspan, f0_0);

# ╔═╡ 6fcabf32-e25a-494b-8871-61c1ba5fef79
begin
	plot(sol_ode_1;
		title="Bacward Euler vs Imporved Euler vs ODE solution",
		label="ODE"
	)
	plot!(sol_euler_3;
		ls=:dashdotdot,
		label="Bacward Euler"
	)
	plot!(sol_euler_2;
		ls=:dash,
		label="Improved Euler"
	)
end

# ╔═╡ 4d69f735-fa09-40b8-b568-7f0c9c959843
md"### The Ralston method
$\begin{cases}
	y_{i+1} = y_i + ({k_1 \over 3} + {2k_2 \over 3})h \\
	k_1 = f(t_i, y_i) \\
	k_2 = f(t_i + {3h \over 4}, y_i + {3hk \over 4})
\end{cases}$
"

# ╔═╡ c5b5138b-5df6-4cfb-a728-1989d5418e3c
function ralston(f, h, tspan, f0)
	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	y_ = [f0, zeros(steps_num)...]
	for i in range(1, steps_num)
		k1 = f(tspan[i], y_[i])
		k2 = f(tspan[i] + 3h/4, y_[i] + 3h*k1/4)
		y_[i+1] = y_[i] + (k1/3 + 2k2/3) * h
	end
	return tspan, y_[1:end-1]
end

# ╔═╡ e38fc5bc-4d7d-4fd8-916c-a08a5e48a557
sol_ralston = ralston(f, h_0, tspan, f0_0);

# ╔═╡ d55f07f6-a8f7-463d-ae20-cb7083229c02
begin
	plot(sol_euler_1;
		title="Ralston vs euler vs ODE solution",
		label="euler"
	)
	plot!(sol_ode_1;
		ls=:dot,
		label="ODE"
	)
	plot!(sol_ralston;
		ls=:dashdotdot,
		label="Ralston"
	)
end

# ╔═╡ efa8ee74-b4d8-47e8-8583-54622d4ce7a5
md"## Runge-Kutta

$y_{i+1} = y_i + h {(k_1 + 2k_2 + 2k_3 + k_4)\over 6}$
$\begin{align}
k_1 &= f(t_i, y_i) \\
k_2 &= f(t_i + {h \over 2}, y_i + {k_1h \over 2}) \\
k_3 &= f(t_i + {h \over 2}, {k_2h \over 2}) \\
k_4 &= f(t_i + {h \over 2}, y_i + k_3h)
\end{align}$
"

# ╔═╡ 2473152b-0fa2-494a-85a0-e9462aaeaf2f
function runge_kutta(f, h, tspan, f0)
	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	y_ = [f0, zeros(steps_num)...]
	for i in range(1, steps_num)
		k1 = f(tspan[i], y_[i])
		k2 = f(tspan[i] + h/2, y_[i] + h * k1/2)
		k3 = f(tspan[i] + h/2, y_[i] + h * k2/2)
		k4 = f(tspan[i] + h/2, y_[i]+ h * k3)
		y_[i+1] = y_[i] + (k1 + 2k2 + 2k3 + k4) * h/6
	end
	return tspan, y_[1:end-1]
end

# ╔═╡ a4fab0a5-5290-4f37-a054-33d9cefe36a7
sol_runge_kutta = runge_kutta(f, h_0, tspan, f0_0);

# ╔═╡ f7a39609-5c16-4441-af5c-7b173f54a8d6
begin
	plot(sol_ode_1;
		title="ODE vs Runge Kutta solution",
		label="ODE"
	)
	plot!(sol_runge_kutta;
		ls=:dashdotdot,
		label="Runge Kutta"
	)
end

# ╔═╡ 4ddb71d0-8fa3-442c-83da-7b7622d2a991
md"### Runge-Kutta-Gill
$y_i = y_i + {1 \over 6} (k_1 + k_4) + {1 \over 3}(bk_2 + dk_3)$

$\begin{align}
k_1 &= hf(t_i, y_i) \\
k_2 &= hf(t_i + {h \over 2}, y_i + {k_1 \over 2}) \\
k_3 &= hf(t_i + {h \over 2}, y_i + ak_1 + bk_2) \\
k_4 &= hf(t_i +h, y_i +cbk_2 + dk_3)
\end{align}$
$a = {{\sqrt{2}-1} \over 2}, \ b = {{2-\sqrt{2}} \over 2}, \ c = - {\sqrt{2} \over 2}, \ d = 1 + {{\sqrt 2}\over 2}$
"

# ╔═╡ 2d8fb7af-c561-424d-b6ca-03a9756b93d5
function runge_kutta_gill(f, h, tspan, f0)
	a = (√2 - 1) / 2
	b = (2 - √2) / 2
	c = √2/2
	d = 1 + √2/2

	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	y_ = [f0, zeros(steps_num)...]
	for i in range(1, steps_num)
		k1 = h * f(tspan[i], y_[i])
		k2 = h * f(tspan[i] + h/2, y_[i] + k1/2)
		k3 = h * f(tspan[i] + h/2, y_[i] + a * k1 + b * k2)
		k4 = h * f(tspan[i] + h, y_[i]+ c * b * k2 + d * k3)
		y_[i+1] = y_[i] + (k1 + k4)/6 + (b * k2 + d * k3)/3
	end
	return tspan, y_[1:end-1]
end

# ╔═╡ 2ce16b6a-bad7-4bca-8f8c-a5d765fa5718
sol_runge_kutta_gill = runge_kutta_gill(f, h_0, tspan, f0_0);

# ╔═╡ f2450fcd-acb2-495c-93dc-bbbeb08baf48
begin
	plot(sol_ode_1;
		title="ODE vs Runge Kutta Gill solution",
		label="ODE"
	)
	plot!(sol_runge_kutta_gill;
		ls=:dashdotdot,
		label="Runge Kutta Gill"
	)
end

# ╔═╡ e59da910-c833-4c1e-b7a5-88f71c8ea3b5
md"# Adams-Bashforth"

# ╔═╡ e9a5623f-7120-449e-abd5-96826100d28d
function adams_bashforth(f, h, tspan, f0)
	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)
	y_ = [f0, zeros(steps_num - 1)...]
	for i in range(1, steps_num - 2)
		# using euler method to get an estimate for y_[i+1]
		if i == 1
			y_[i+1] = y_[i] + h * f(tspan[i], y_[i])
		end
		y_[i+2] = y_[i+1] + (3 * h * f(tspan[i+1], y_[i+1]) - h * f(tspan[i], y_[i]))/2
	end
	return tspan[1:end-1], y_[1:end-1]
end

# ╔═╡ 59d76095-6ede-42da-b220-0054e0686a9c
sol_adams_bashforth = adams_bashforth(f, h_0, tspan, f0_0);

# ╔═╡ ac6ac5f2-fc5b-41e6-a88a-c2e4f6e64d64
begin
	plot(sol_ode_1;
		title="ODE vs Adam-Bashforth vs Euler solution",
		label="ODE"
	)
	plot!(sol_euler_1;
		ls=:dot,
		label="Euler"
	)
	plot!(sol_adams_bashforth;
		ls=:dashdotdot,
		label="Adams Bashforth"
	)
end

# ╔═╡ 3277d14e-5ae8-4364-9380-c0f4a530d4e9
md"# Milne
$y_{i+3} = y_i + {4h \over 3}(2f(t_{i+2}, y_{i+2}) - f(t_{i+1}, y_{i+1}) + 2f(t_i, y_i))$

For Milne we have to use other alogrithms to perdict $y_{2}$, and $y_{i+2}$. We are using **Forward Euler** to predict $y_{2}$, and **Adams-Bashforth** to predict $y_{i+2}$.
"

# ╔═╡ a3461d3b-fe86-4e19-ad7a-009155df4113
function milne(f, h, tspan, f0)
	tspan = collect(tspan[1]:h:tspan[2])
	steps_num = length(tspan)

	# using euler method to get an estimate for y_[2]
	f1 = f0 + h * f(tspan[1], f0)
	y_ = [f0, f1, zeros(steps_num - 2)...]
	for i in range(1, steps_num - 3)
		# using adams-bashforth to get y_[i+2]
		y_[i+2] = y_[i+1] + (3 * h * f(tspan[i+1], y_[i+1]) - h * f(tspan[i], y_[i]))/2
		y_[i+3] = y_[i] + (4h/3) * (2f(tspan[i+2], y_[i+2]) - f(tspan[i+1], y_[i+1]) + 2f(tspan[i], y_[i]))
	end
	return tspan[1:end-1], y_[1:end-1]
end

# ╔═╡ fc512ce3-d49e-4b38-a82f-9feaf2067ca8
sol_milne = milne(f, h_0, tspan, f0_0);

# ╔═╡ 456309e4-78e2-427d-bd80-50787694e971
begin
	plot(sol_ode_1;
		title="ODE vs Milne vs Euler solution",
		label="ODE"
	)
	plot!(sol_euler_1;
		ls=:dot,
		label="Euler"
	)
	plot!(sol_milne;
		ls=:dashdotdot,
		label="Milne"
	)
end

# ╔═╡ e0bb8e3e-893c-4bc8-87f5-c00b7c7fbb2e
md" # Translation between MatLab and Julia ODE solvers
![Imgur](https://i.imgur.com/PCimdOA.png)
#

![Imgur](https://i.imgur.com/GObvC8H.png)"

# ╔═╡ c11023dc-c05d-43e6-9a1e-99c3b014db4f
md"
!!! info
	The rest of the chapter is using stanadard MatLab ODE sovlers,
	so this means I will do the same using julia, I have already used it extensively in the first chapter so I skipped this section.
"

# ╔═╡ Cell order:
# ╟─4a57db60-efcf-11ed-34d4-690120ffb2c4
# ╠═3a12645a-ff8f-4540-a1b2-1f75f8a4e3f9
# ╠═f82bcf92-9501-4aac-8399-2659fdbd9593
# ╠═641fad95-8abb-4818-a297-a4244e261283
# ╟─11cee654-fe7b-4376-9a51-5d6ae99ef447
# ╠═3e883f70-3061-4588-a3fd-86c651c294c0
# ╠═015dcf60-51dc-4d6e-ac07-0f2c847067e7
# ╠═02d916f4-1fa3-4fbf-b118-ef9b30cf40d6
# ╟─1131ac8d-42a3-4777-a538-6b89a805e13f
# ╠═d8cc309f-1035-4b7b-9215-af465eb8fd9b
# ╠═ecff5b0c-5474-4094-ab0b-36aac496bca6
# ╟─6e685a10-a662-4a9d-b9e0-102a20919ea2
# ╠═478f44f4-6470-4abf-9258-da6908b02a99
# ╠═1cfe11dc-e7e8-4fe6-9e5b-e3ddb2f09964
# ╟─67c546d1-6d5e-4a3b-9730-28f5d00a2c6e
# ╠═ad9fd0c7-4172-4d40-93bd-3e83c1689b54
# ╟─0664d3fa-89f2-4338-8240-e0aa9c0866d8
# ╠═4c4acf15-99b2-41ac-b9b1-31557c6393a2
# ╠═2d6cc393-2c31-4699-bb2f-14d5014b42ca
# ╟─6fcabf32-e25a-494b-8871-61c1ba5fef79
# ╠═4d69f735-fa09-40b8-b568-7f0c9c959843
# ╠═c5b5138b-5df6-4cfb-a728-1989d5418e3c
# ╠═e38fc5bc-4d7d-4fd8-916c-a08a5e48a557
# ╠═d55f07f6-a8f7-463d-ae20-cb7083229c02
# ╟─efa8ee74-b4d8-47e8-8583-54622d4ce7a5
# ╠═2473152b-0fa2-494a-85a0-e9462aaeaf2f
# ╠═a4fab0a5-5290-4f37-a054-33d9cefe36a7
# ╠═f7a39609-5c16-4441-af5c-7b173f54a8d6
# ╟─4ddb71d0-8fa3-442c-83da-7b7622d2a991
# ╠═2d8fb7af-c561-424d-b6ca-03a9756b93d5
# ╠═2ce16b6a-bad7-4bca-8f8c-a5d765fa5718
# ╠═f2450fcd-acb2-495c-93dc-bbbeb08baf48
# ╟─e59da910-c833-4c1e-b7a5-88f71c8ea3b5
# ╠═e9a5623f-7120-449e-abd5-96826100d28d
# ╠═59d76095-6ede-42da-b220-0054e0686a9c
# ╠═ac6ac5f2-fc5b-41e6-a88a-c2e4f6e64d64
# ╟─3277d14e-5ae8-4364-9380-c0f4a530d4e9
# ╠═a3461d3b-fe86-4e19-ad7a-009155df4113
# ╠═fc512ce3-d49e-4b38-a82f-9feaf2067ca8
# ╠═456309e4-78e2-427d-bd80-50787694e971
# ╟─e0bb8e3e-893c-4bc8-87f5-c00b7c7fbb2e
# ╟─c11023dc-c05d-43e6-9a1e-99c3b014db4f

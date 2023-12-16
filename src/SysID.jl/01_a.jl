### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 5bda20ac-95a2-11ee-0f56-851af662b57a
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ 1757b3df-36f3-4bac-9465-5c126e3f553a
using Random, Distributions, LinearAlgebra, Plots, Unitful

# ╔═╡ 1beef37b-ddcc-4b28-8c32-e0363e47386b
md"# Exercise 1.a: Least squares estimation of the value of a resistor"

# ╔═╡ 5f33232b-a16e-4aa1-9993-35f79ed601a4
begin
	R₀ = 1000u"Ω"
	iₘₐₓ = 0.01u"A"
	current_distribution = Uniform(-iₘₐₓ.val,iₘₐₓ.val)
	nᵤ = Normal(0, 1)  # voltage disturbance
	num_of_measurements = [10 100 1000 10_000]
	num_repeations = 100

	R̂ = zeros((num_repeations, length(num_of_measurements),))u"Ω"

	for n in 1:length(num_of_measurements)
		for r in 1:num_repeations
			i = rand(current_distribution, num_of_measurements[n])u"A"
			u = R₀ * i + rand(nᵤ, num_of_measurements[n])u"V"
			# the `ustrip` is only necessary because unitful hasn't been proberly
			# integraded with linear solvers.
			R̂[r, n] =	ustrip(i) \ ustrip(u)u"Ω"
		end
	end
end

# ╔═╡ e2b6fffc-388e-4bf8-8cf5-643fbd3ed34f
begin
	plot(R̂;
		layout=4,
		seriestype=:scatter,
		ylims=(900, 1100),
		legend=false,
		title=["$r Measurements" for r in num_of_measurements]
	)
	plot!(fill(R₀, size(R̂)); linewidth=3)
end

# ╔═╡ Cell order:
# ╠═5bda20ac-95a2-11ee-0f56-851af662b57a
# ╠═1757b3df-36f3-4bac-9465-5c126e3f553a
# ╟─1beef37b-ddcc-4b28-8c32-e0363e47386b
# ╠═5f33232b-a16e-4aa1-9993-35f79ed601a4
# ╟─e2b6fffc-388e-4bf8-8cf5-643fbd3ed34f

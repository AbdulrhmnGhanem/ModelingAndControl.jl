### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 5bda20ac-95a2-11ee-0f56-851af662b57a
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.activate("../..")
end

# ╔═╡ 1757b3df-36f3-4bac-9465-5c126e3f553a
using Random, Distributions, LinearAlgebra, Plots, Unitful

# ╔═╡ 1beef37b-ddcc-4b28-8c32-e0363e47386b
md"# Exercise 5.a: Weighted least squares estimation"

# ╔═╡ 5f33232b-a16e-4aa1-9993-35f79ed601a4
begin
	R₀ = 1000u"Ω"
	iₘₐₓ = 0.01u"A"
	current_distribution = Uniform(-iₘₐₓ.val,iₘₐₓ.val)
	nᵤ1 = Normal(0, 1)
	nᵤ2 = Normal(0, 4)
	num_of_measurements = 100
	num_repeations = Int(10e5)
	
	R̂ = zeros((num_repeations, 2))u"Ω"
	w = [fill(nᵤ1.:σ, num_of_measurements)... fill(nᵤ2.:σ^2, num_of_measurements)...]'
	
	
	for r in 1:num_repeations
		i1 = rand(current_distribution, num_of_measurements)u"A"
		i2 = rand(current_distribution, num_of_measurements)u"A"
		u1 = R₀ * i1 + rand(nᵤ1, num_of_measurements)u"V"
		u2 = R₀ * i2 + rand(nᵤ2, num_of_measurements)u"V"
		i = [i1... i2...]'
		u = [u1... u2...]'
		
		R̂[r, 1] = (ustrip(i) \ ustrip(u))[1]u"Ω"
		R̂[r, 2] = (ustrip((u' * (i ./ w)) / (i' * (i ./ w))))[1]u"Ω"
	end
end

# ╔═╡ e2b6fffc-388e-4bf8-8cf5-643fbd3ed34f
stephist(R̂;
	normalize=:pdf,
	ylims=(0, 0.025),
	labels=reshape(["LS", "WLS"], 1, :),
	title="LS vs WLS",
)

# ╔═╡ Cell order:
# ╠═5bda20ac-95a2-11ee-0f56-851af662b57a
# ╠═1757b3df-36f3-4bac-9465-5c126e3f553a
# ╟─1beef37b-ddcc-4b28-8c32-e0363e47386b
# ╠═5f33232b-a16e-4aa1-9993-35f79ed601a4
# ╟─e2b6fffc-388e-4bf8-8cf5-643fbd3ed34f

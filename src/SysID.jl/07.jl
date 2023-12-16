### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 1aaf132a-96b3-11ee-1a47-7fe4b07f94c5
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ c05f3ace-fb4c-472b-ae6a-b7bcafb102b3
using Random, Distributions, LinearAlgebra, Plots

# ╔═╡ a6fd6c00-113c-4ff7-84ea-4930574ee6dc
md"# Exercise 7: Characterizing 2-dimensional parameter estimate"

# ╔═╡ bc1962a7-fa1a-4a1d-ad54-6e6c522e1dd3
begin
	num_of_measurements = 100
	num_repeations = Int(10e3)
	u1₀ = LinRange(-3, 3, num_of_measurements)
	u2₀ = LinRange(2, 5, num_of_measurements)
	
	a = 0.1
	nₜ = Normal(0, 1)
	y1₀ = a * u1₀
	y2₀ = a * u2₀

	ŷ1 = zeros((num_repeations, 2))
	ŷ2 = zeros((num_repeations, 2))

	for i in 1:num_repeations
		y1 = y1₀ + rand(nₜ, num_of_measurements)
		y2 = y2₀ + rand(nₜ, num_of_measurements)

		K1 = [u1₀ ones(size(u1₀))]
		K2 = [u2₀ ones(size(u2₀))]
		ŷ1[i, :] = K1 \ y1
		ŷ2[i, :] = K2 \ y2
	end
end

# ╔═╡ 83a443e5-8ab9-401a-8e83-feb9d1d45d96
md"
!!! warning
	What is modeled output? and add its visualization.
"

# ╔═╡ 38270932-1210-42fd-b2c2-721dbda5cd73
begin
	plot([ŷ1[:, 1], ŷ2[:, 1]], [ŷ1[:, 2], ŷ2[:, 2]]; 
		seriestype=:scatter,
		labels=reshape(["u₀ ∈ [-3,3]", "u₀ ∈ [2, 5]"], 1, :),
		xlabel="Slope",
		ylabel="Offset",
		ylims=(-2, 2),
		xlims=(-0.5, 0.5)
	)
end

# ╔═╡ 214f5024-2a01-4bdc-8c0d-e751d09b4fbf
md"### Covariance matrices"

# ╔═╡ 632810ed-4f25-4172-b3a7-3f43acdaa47d
cov(ŷ1), cov(ŷ2)

# ╔═╡ 90e15eb2-224b-4ef2-b976-3a76289986bb
md"### Correlation matrices"

# ╔═╡ f13e8e55-0d89-4dfa-b4b9-6567bfb935e4
cor(ŷ1), cor(ŷ2)

# ╔═╡ Cell order:
# ╠═1aaf132a-96b3-11ee-1a47-7fe4b07f94c5
# ╠═c05f3ace-fb4c-472b-ae6a-b7bcafb102b3
# ╟─a6fd6c00-113c-4ff7-84ea-4930574ee6dc
# ╠═bc1962a7-fa1a-4a1d-ad54-6e6c522e1dd3
# ╟─83a443e5-8ab9-401a-8e83-feb9d1d45d96
# ╟─38270932-1210-42fd-b2c2-721dbda5cd73
# ╟─214f5024-2a01-4bdc-8c0d-e751d09b4fbf
# ╟─632810ed-4f25-4172-b3a7-3f43acdaa47d
# ╟─90e15eb2-224b-4ef2-b976-3a76289986bb
# ╟─f13e8e55-0d89-4dfa-b4b9-6567bfb935e4

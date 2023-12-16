### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ ad8895e8-979b-11ee-3913-5d95eb77891a
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ 41f31b91-41b9-4b87-a496-776daff525b6
using Random, Distributions, LinearAlgebra, Plots, Unitful

# ╔═╡ f4adbb3b-a04f-4bcc-87fb-cbccb16e83c6
md"# Exercise 10: Influence of the number of the parameters on the model uncertainty

!!! purpose
	To show that when we have a prior knowledge about the model (e.g., in this exercise at `t=0`, `y=0`)
	increasing the number of the parameters, increases the uncertainty.
    But doesn't necessarily introduce systematic error.
"

# ╔═╡ 67aff6ab-1860-46e4-b2d4-61f0aab6b427
begin
	a₀ = 1
	N = 1000
	nₜ = Normal(0, 1)
	num_of_repeations = Int(10e5)

	t = LinRange(0, 1, N)
	y₀ = a₀ * t  # the true system

	ŷ₁ = zeros((num_of_repeations, 2))
	ŷ₂ = zeros(num_of_repeations)

	Threads.@threads for i in 1:num_of_repeations
		y = y₀ + rand(nₜ, N)

		k1 = [t ones(size(t))]
		k2 = t

		@inbounds ŷ₁[i, :] = k1 \ y
		@inbounds ŷ₂[i] = k2 \ y
	end

	means = round.(mean.([ŷ₁[:, 1], ŷ₂]); digits=4)
	stds = round.(std.([ŷ₁[:, 1], ŷ₂]); digits=4)
end;

# ╔═╡ 6f276e54-df15-465e-a1b4-b86d39b32bb7
md"
|    |Two-Parameter Model|One-Parameter Model|
|----|-------------------|-------------------|
| μ  |      $(means[1])  |   $(means[2])     |
| σ  |      $(stds[1])   |   $(stds[2])      |
"

# ╔═╡ 934a56b0-84ba-40fa-8509-d7bb0f4afbae
stephist([ŷ₁[:,1] ŷ₂];
	xlabel="Slope",
	labels=reshape(["Two-Parameter Model (at+b)", "One-Parameter Model (at)"], 1, :),
	seriestype=:scatter,
	normalize=:pdf,
)

# ╔═╡ Cell order:
# ╟─ad8895e8-979b-11ee-3913-5d95eb77891a
# ╠═41f31b91-41b9-4b87-a496-776daff525b6
# ╟─f4adbb3b-a04f-4bcc-87fb-cbccb16e83c6
# ╠═67aff6ab-1860-46e4-b2d4-61f0aab6b427
# ╟─6f276e54-df15-465e-a1b4-b86d39b32bb7
# ╟─934a56b0-84ba-40fa-8509-d7bb0f4afbae

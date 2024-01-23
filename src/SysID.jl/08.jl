### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ e29a8f8a-9736-11ee-0308-19e029629fa5
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ 61158423-326c-472b-9de6-9ad5319ae1ed
using Random, Distributions, LinearAlgebra, Plots, Unitful

# ╔═╡ f03d2f62-761e-4863-b4b4-b2506862e357
md"# Exercise 8: Dependence of the optimal cost function on the distribution of disturbing noise

!!! purpose
	The purpose of this exercise to show that for each noise distribution, there corresponds an optimal choice of the cost function.

	Normal → LS

	Laplace → LAV
"

# ╔═╡ f8cfa388-a043-4a33-a462-3cdd4d0e6df2
md"
!!! warning
	In the book the Laplacian noise is defined with `μ=0` and `σ=1` while in the MATLAB code it's defined with `μ=0` and `b=1`, these are two different distributions.

	A laplace distribution with `μ=0` and `σ=1` is equivalent to `μ=0` and `b=1/√2`.
	So the normal and laplace distributions in this exercise don't have the same standard deviation.
"

# ╔═╡ 2ef4ec2c-d649-44f2-9a40-5664b51f2ebe
begin
    R₀ = 1000u"Ω"
    i₀ = 0.01u"A"
    ĩ = Uniform(0, ustrip(i₀))
    Normal_nᵤ = Normal(0, 1)
    # for the laplace distribution with σ=1 use this ↓ instead
    # Laplace_nᵤ = Laplace(0, 1/√2)
    Laplace_nᵤ = Laplace(0, 1)
    num_of_measurements = 100
    num_repeations = Int(10e4)

    R̂1 = zeros((num_repeations, 2))u"Ω"
    R̂2 = zeros((num_repeations, 2))u"Ω"

    Vₗₐᵥ(R, u, i) = sum(abs.(u - R * i)) / num_of_measurements
    Ω = R₀ * (0.9:0.001:1.1)

    Threads.@threads for r = 1:num_repeations
        i = rand(ĩ, num_of_measurements)u"A"

        uₙ = R₀ * i + rand(Normal_nᵤ, num_of_measurements)u"V"
        R̂1[r, 1] = ustrip(i) \ ustrip(uₙ)u"Ω"  # LS
        R̂1[r, 2] = argmin(r -> Vₗₐᵥ(r, uₙ, i), Ω)  # LAV

        uₗ = R₀ * i + rand(Laplace_nᵤ, num_of_measurements)u"V"
        R̂2[r, 1] = ustrip(i) \ ustrip(uₗ)u"Ω"  # LS
        R̂2[r, 2] = argmin(r -> Vₗₐᵥ(r, uₗ, i), Ω)  # LAV
    end
end

# ╔═╡ b5dfea80-4ae9-470a-aa59-6391753227a3
stephist(
    [R̂1 R̂2];
    layout = (2, 1),
    normalize = :pdf,
    ylims = (0, 0.04),
    labels = reshape(["LS", "LS", "LAV", "LAV"], 1, :),
    titles = reshape(["Guassian Noise", "Laplacian Noise"], 1, :),
)

# ╔═╡ Cell order:
# ╠═e29a8f8a-9736-11ee-0308-19e029629fa5
# ╠═61158423-326c-472b-9de6-9ad5319ae1ed
# ╟─f03d2f62-761e-4863-b4b4-b2506862e357
# ╟─f8cfa388-a043-4a33-a462-3cdd4d0e6df2
# ╠═2ef4ec2c-d649-44f2-9a40-5664b51f2ebe
# ╟─b5dfea80-4ae9-470a-aa59-6391753227a3

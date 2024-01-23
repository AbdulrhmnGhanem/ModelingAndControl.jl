### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 057dc498-9761-11ee-2942-a58010fc0fe8
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ e6323fa8-7bde-4beb-bd77-972215bd2d13
using Random, Distributions, LinearAlgebra, Plots, Unitful

# ╔═╡ 4e035c3b-ab6e-4352-a04c-8301052310cd
md"# Exercise 9: Identification in the presence of outiliers

!!! purpose
	Like the previous exercise, the purpose of this exercise to show that for each calibration method (dealing with outliers), there corresponds an optimal choice of the cost function.

	Mean calibration → LS

	Median calibration → LAV
"

# ╔═╡ 1e55a00d-2f3d-494c-9613-1196b744c719
begin
    R₀ = 1000u"Ω"
    i₀ = 0.01u"A"
    ĩ = Uniform(0, ustrip(i₀))
    Χ² = Chisq(1)  # voltage noise distribution
    num_of_measurements = 100
    num_repeations = Int(10e4)


    χ² = rand(Χ², num_repeations * 100)
    Χ²_μ, Χ̃² = mean(χ²)u"V", median(χ²)u"V"

    Vₗₐᵥ(R, u, i) = sum(abs.(u - R * i)) / num_of_measurements
    Ω = R₀ * (0.8:0.001:1.1)


    R̂1 = zeros((num_repeations, 2))u"Ω"
    R̂2 = zeros((num_repeations, 2))u"Ω"

    Threads.@threads for r = 1:num_repeations
        i = rand(ĩ, num_of_measurements)u"A"

        uₙ = R₀ * i + rand(Χ², num_of_measurements)u"V"
        uₙ_calibrated_with_mean = uₙ .- Χ²_μ
        uₙ_calibrated_with_median = uₙ .- Χ̃²

        R̂1[r, 1] = ustrip(i) \ ustrip(uₙ_calibrated_with_mean)u"Ω"  # LS
        R̂1[r, 2] = argmin(r -> Vₗₐᵥ(r, uₙ_calibrated_with_mean, i), Ω)  # LAV

        R̂2[r, 1] = ustrip(i) \ ustrip(uₙ_calibrated_with_median)u"Ω"  # LS
        R̂2[r, 2] = argmin(r -> Vₗₐᵥ(r, uₙ_calibrated_with_median, i), Ω)  # LAV
    end
end

# ╔═╡ 72683c53-ee97-4b4a-a65d-dfcd6b140459
stephist(
    [R̂1 R̂2];
    layout = (2, 1),
    normalize = :pdf,
    ylims = (0, 0.04),
    xlims = (800, 1200),
    labels = reshape(["LS", "LS", "LAV", "LAV"], 1, :),
    titles = reshape(["Calibration with mean", "Calibration with median"], 1, :),
)

# ╔═╡ Cell order:
# ╟─057dc498-9761-11ee-2942-a58010fc0fe8
# ╠═e6323fa8-7bde-4beb-bd77-972215bd2d13
# ╟─4e035c3b-ab6e-4352-a04c-8301052310cd
# ╠═1e55a00d-2f3d-494c-9613-1196b744c719
# ╟─72683c53-ee97-4b4a-a65d-dfcd6b140459

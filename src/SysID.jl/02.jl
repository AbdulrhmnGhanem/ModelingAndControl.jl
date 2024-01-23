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
md"# Exercise 2: Study of the asymptotic distribution of an estimate"

# ╔═╡ 5f33232b-a16e-4aa1-9993-35f79ed601a4
begin
    R₀ = 1000u"Ω"
    i₀ = 0.01u"A"
    σᵤ = 0.2
    nᵤ1 = Normal(0, σᵤ)  # voltage disturbance
    nᵤ2 = Uniform(-√3σᵤ, √3σᵤ)
    num_of_measurements = [1 2 4 8]
    num_repeations = Int(10e5)

    R̂1 = zeros((num_repeations, length(num_of_measurements)))u"Ω"
    R̂2 = zeros((num_repeations, length(num_of_measurements)))u"Ω"

    for n = 1:length(num_of_measurements)
        for r = 1:num_repeations
            i = fill(i₀, num_of_measurements[n])
            u1 = R₀ * i + rand(nᵤ1, num_of_measurements[n])u"V"
            u2 = R₀ * i + rand(nᵤ2, num_of_measurements[n])u"V"
            R̂1[r, n] = ustrip(i) \ ustrip(u1)u"Ω"
            R̂2[r, n] = ustrip(i) \ ustrip(u2)u"Ω"
        end
    end
end

# ╔═╡ 3411d10f-26b4-498e-a099-fedaaa1d60bd
md"$Std(nᵤ2) = \sqrt{\frac{(b-a)^2}{12}}$
$\sqrt{\frac{(0.2\sqrt{3}+0.2\sqrt{3})^2}{12}} = 0.2 = \sigma$
"

# ╔═╡ ffa6f39c-92c3-4318-9fd1-dfeb9c4c1dac
stephist(
    [R̂1 R̂2];
    layout = 4,
    normalize = :pdf,
    ylims = (0, 0.06),
    fillalpha = 0.0,
    labels = reshape([repeat(["uniform"], 4)..., repeat(["normal"], 4)...], 1, :),
    title = ["$r Measurements" for r in num_of_measurements],
)

# ╔═╡ Cell order:
# ╠═5bda20ac-95a2-11ee-0f56-851af662b57a
# ╠═1757b3df-36f3-4bac-9465-5c126e3f553a
# ╟─1beef37b-ddcc-4b28-8c32-e0363e47386b
# ╠═5f33232b-a16e-4aa1-9993-35f79ed601a4
# ╟─3411d10f-26b4-498e-a099-fedaaa1d60bd
# ╟─ffa6f39c-92c3-4318-9fd1-dfeb9c4c1dac

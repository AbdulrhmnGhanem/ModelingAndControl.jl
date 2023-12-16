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
md"# Exercise 4: Importance of the choice of the independent variable or input"

# ╔═╡ 5f33232b-a16e-4aa1-9993-35f79ed601a4
begin
	R₀ = 1000u"Ω"
	i₀ = 0.01u"A"
	σᵤ = 1
	nᵤ = Normal(0, σᵤ)  # voltage disturbance
	num_of_measurements = 100
	num_repeations = Int(10e5)
	
	R̂ = zeros((num_repeations, 2))u"Ω"
	
	for r in 1:num_repeations
		i = fill(i₀, num_of_measurements)
		u = R₀ * i + rand(nᵤ, num_of_measurements)u"V"
		R̂[r,1] =	ustrip(i) \ ustrip(u)u"Ω"
		R̂[r, 2] = 1/ ((u' * i) / (u' * u))
	end
end

# ╔═╡ 3411d10f-26b4-498e-a099-fedaaa1d60bd
md"
!!! note
	The signal with the highest SNR should be used as independent variable in order to reduce the systematic error.
	
	The bias is proportional to the inverse of SNR ($\frac{noise \ power}{signal \ power})$
"

# ╔═╡ ffa6f39c-92c3-4318-9fd1-dfeb9c4c1dac
stephist(R̂;
	normalize=:pdf,
	ylims=(0, 0.05),
	fillalpha=0.0,
	labels=reshape(["regressor = i", "regressor = v"], 1, :),
	title="Choice of regressor variable",
)

# ╔═╡ Cell order:
# ╠═5bda20ac-95a2-11ee-0f56-851af662b57a
# ╠═1757b3df-36f3-4bac-9465-5c126e3f553a
# ╟─1beef37b-ddcc-4b28-8c32-e0363e47386b
# ╠═5f33232b-a16e-4aa1-9993-35f79ed601a4
# ╟─3411d10f-26b4-498e-a099-fedaaa1d60bd
# ╟─ffa6f39c-92c3-4318-9fd1-dfeb9c4c1dac

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
using Random, Distributions, LinearAlgebra, Plots, Unitful, Statistics

# ╔═╡ 1beef37b-ddcc-4b28-8c32-e0363e47386b
md"# Impact of noise on the regressor (input) measurements"

# ╔═╡ 5f33232b-a16e-4aa1-9993-35f79ed601a4
begin
	R₀ = 1000u"Ω"
	i₀ = 0.01u"A"
	σᵤ = 1
	nᵤ1 = Normal(0, σᵤ)  # voltage disturbance
	i_distribution = Uniform(-10e-3, 10e-3)
	nᵢs = [Normal(0, σᵢ) for σᵢ in (0, 0.5e-3, 1e-3)]
	num_of_measurements = 100
	num_repeations = Int(10e5)
	
	R̂1 = zeros((num_repeations, length(nᵢs)))u"Ω"
	R̂2 = zeros((num_repeations, length(nᵢs)))u"Ω"
	
	for nᵢ in 1:length(nᵢs)
		for r in 1:num_repeations
			i = rand(i_distribution, num_of_measurements)u"A"
			noisy_i = i + rand(nᵢs[nᵢ], num_of_measurements)u"A"
			u = R₀ * i + rand(nᵤ1, num_of_measurements)u"V"
			R̂1[r, nᵢ] =	ustrip(i) \ ustrip(u)u"Ω"
			R̂2[r, nᵢ] =	ustrip(noisy_i) \ ustrip(u)u"Ω"
		end
	end
end

# ╔═╡ b2fd135f-e925-4e27-81fc-f483be7f5126
begin
	s = [
		["σ = $(n.:σ) Ω" for n in nᵢs],
		["σ = $(round(typeof(1u"Ω"), s;))" for s in std(R̂2;dims=1)],
		["μ = $(round(typeof(1u"Ω"), m))" for m in mean(R̂2;dims=1)],	
	]
	
	for i in 1:3
		@info ("when `I` noise has $(s[1][i]) →  $(s[2][i]), $(s[3][i])")
	end
end

# ╔═╡ ffa6f39c-92c3-4318-9fd1-dfeb9c4c1dac
stephist([R̂1 R̂2];
	layout=(3, 1),
	normalize=:pdf,
	ylims=(0, 0.025),
	fillalpha=0.0,
	labels=reshape([repeat(["noise on V"], 3)..., repeat(["noise on V and I"], 3)...], 1, :),
)

# ╔═╡ Cell order:
# ╠═5bda20ac-95a2-11ee-0f56-851af662b57a
# ╠═1757b3df-36f3-4bac-9465-5c126e3f553a
# ╟─1beef37b-ddcc-4b28-8c32-e0363e47386b
# ╠═5f33232b-a16e-4aa1-9993-35f79ed601a4
# ╠═b2fd135f-e925-4e27-81fc-f483be7f5126
# ╟─ffa6f39c-92c3-4318-9fd1-dfeb9c4c1dac

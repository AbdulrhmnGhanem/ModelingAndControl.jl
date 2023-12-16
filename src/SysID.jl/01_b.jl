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
md"# Exercise 1.b: Analysis of the standard deviation"

# ╔═╡ 5f33232b-a16e-4aa1-9993-35f79ed601a4
begin
	R₀ = 1000u"Ω"
	i₀ = 0.01u"A"
	nᵤ = Normal(0, 1)  # voltage disturbance
	num_of_measurements = [10 100 1000 10_000]
	num_repeations = 1000
	
	R̂ = zeros((num_repeations, length(num_of_measurements),))u"Ω"
	
	for n in 1:length(num_of_measurements)
		for r in 1:num_repeations
			i = fill(i₀, num_of_measurements[n])
			u = R₀ * i + rand(nᵤ, num_of_measurements[n])u"V"
			R̂[r, n] =	ustrip(i) \ ustrip(u)u"Ω"
		end
	end
end

# ╔═╡ 51c027b5-504b-4fe3-be44-05528dc900f6
begin
	emperical_std = std(R̂; dims=[1])
	theoritical_std_f(N) = (nᵤ.:σ)u"V" / (sqrt(N) * i₀)
	theoritical_std = theoritical_std_f.(1:10e4)
end;

# ╔═╡ b339d553-9821-4ce1-9211-29f25c3c4732
md"### The standard deviation is proportional to $-\frac{1}{\sqrt N}$"

# ╔═╡ 54e3a4e9-d05d-462c-83e3-135dbf344198
begin
	plot(num_of_measurements, emperical_std;
		seriestype=:scatter, 
		xaxis=:log,
		yaxis=:log,
		legend=false,
		title="Std vs number of experiments",
		xlabel="No. experiments",
	)
	plot!(theoritical_std, ylabel="σ",)
end

# ╔═╡ Cell order:
# ╠═5bda20ac-95a2-11ee-0f56-851af662b57a
# ╠═1757b3df-36f3-4bac-9465-5c126e3f553a
# ╟─1beef37b-ddcc-4b28-8c32-e0363e47386b
# ╠═5f33232b-a16e-4aa1-9993-35f79ed601a4
# ╠═51c027b5-504b-4fe3-be44-05528dc900f6
# ╟─b339d553-9821-4ce1-9211-29f25c3c4732
# ╟─54e3a4e9-d05d-462c-83e3-135dbf344198

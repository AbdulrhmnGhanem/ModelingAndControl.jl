### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ dc2d1942-accd-11ee-34f8-ef2f689c65b2
# ╠═╡ show_logs = false
begin 
 # If you are running this notebook as a stannalone notebook disable this cell.
 import Pkg 
 Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 6ce6c76b-c24b-4a6b-8297-f62cdce1191e
using Plots, FFTW, DSP, Distributions

# ╔═╡ 40338bfb-3c0d-4990-b84b-0ca6d42373c0
md"## Exercise 44: Analysis of the noise errors on FRF measurements

!!! warning
	I haven't finish this exercise yet!
"

# ╔═╡ 0345624c-b9f8-4a13-9371-f08f428523ce
function cheby1(n, r, wp)
    digitalfilter(Lowpass(wp), Chebyshev1(n, r))
end

# ╔═╡ 7c4713b0-e397-4f06-9abf-77e73480a671
function random_phase_mutli_sine(lines, period, N)
	u_all = zeros(period, N ÷ period)
	for r in 1:N÷period
		U = zeros(Complex{Float64}, period, 1)
		U[lines] = exp.(2π * im * rand(length(lines), 1))
		u = 2real(ifft(U))
		u_all[:, r] = u / std(u)
	end
	vec(reshape(u_all, N, 1))
end

# ╔═╡ 3d1a9815-e856-4c2d-a7ab-2ff2d5daadd2
begin
	period = 1024
	fₛ = 128
	t = LinRange(0, (period - 1) / fₛ, period)
	Lines = 1:period÷2
	m = 16
	N_experiments = 10
	A = 0.33
	σ = 0.02
	Nₜᵣₐₙₛ = period

	filter_order = 2
	filter_resonance = 10
	filter_cutoff = 0.1fₛ

	h = cheby1(filter_order, filter_resonance, 2filter_cutoff / fₛ)
	H, w = freqresp(h)


	GDiff = zeros(ComplexF64, length(Lines)-1, N_experiments)
	G₀Diff = similar(GDiff)
	GSin = zeros(ComplexF64, length(Lines), N_experiments)
	CUY = zeros(ComplexF64, length(Lines)-1, 1)
	for i in 1:N_experiments
		N = m * period
		# Excite the system with random input
		u = A * randn(N + Nₜᵣₐₙₛ)
		y₀ = filt(h, u)
		# Add noise to the output
		yₙ = σ * randn(size(y₀))
		y = y₀ + yₙ
		# Remove transient
		deleteat!(u, 1:Nₜᵣₐₙₛ) 
		deleteat!(y₀, 1:Nₜᵣₐₙₛ)
		deleteat!(y, 1:Nₜᵣₐₙₛ)
		# Average over `m` realization
		u = reshape(u, period, m)
		y = reshape(y, period, m)
		y₀ = reshape(y₀, period, m)
	
		# FRF
		U = fft(u, 1) / (period/2)
		U = U[Lines, :]
		U = diff(U, dims=1)
	
		Y = fft(y, 1) / (period/2)
		Y = Y[Lines, :]
		Y = diff(Y, dims=1)
	
		Y₀ = fft(y₀, 1) / (period/2)
		Y₀ = Y₀[Lines, :]
		Y₀ = diff(Y₀, dims=1)
	
		UU = mean(abs2.(U[:, 1:m]), dims=2)
		YY = mean(abs2.(Y[:, 1:m]), dims=2)
		YU = mean(Y[:, 1:m] .* conj.(U[:, 1:m]), dims=2)
		global CUY = abs2.(YU)  ./ (abs.(UU) .* abs.(YY))
		Y₀U = mean(Y₀[:, 1:m] .* conj.(U[:, 1:m]), dims=2)
		GDiff[:, i] = YU ./ UU
		G₀Diff[:, i] = Y₀U ./ UU
	end

	Var_GDiff = abs.(GDiff[:,end]) .^ 2 .* (1 .- CUY) ./ CUY ./ m
	σ_G₀Diff = sqrt.(Var_GDiff)
	μ_GDiff = mean(GDiff, dims=2)
	σ_GDiff = mapslices(x -> std(x, dims=1), GDiff, dims=2)
end;

# ╔═╡ a381b387-1ff5-484a-ad28-a2df6e529aea
begin
	mapped_f =  w * (fₛ / (2 * π))
	db(s) = amp2db.(abs.(s))
	p1 = scatter(mapped_f, 
		db.([H, μ_GDiff[1:2:end] - H[1:end-1], σ_GDiff[1:end-1] / sqrt(N_experiments)]);
		ylims=(-60, 1),
		labels=reshape(["G₀Diff", "|GDiff - G₀Diff|", "measured σ"], 1, :),
		title="Noise excitation (averaged)",
		ylabel="Amplitude (dB)",
	)
	p2 = scatter(mapped_f, db.([H, σ_GDiff, σ_G₀Diff]);
		ylims=(-60, 1),
		labels=reshape(["G₀Diff", "measured σ", "theoretical σ"], 1, :),
		title="Noise excitation (single realization)",
		ylabel="Amplitude (dB)",
		xlabel="f (Hz)"
	)

	plot(p1, p2; layout=(2,1))
end

# ╔═╡ Cell order:
# ╠═dc2d1942-accd-11ee-34f8-ef2f689c65b2
# ╠═6ce6c76b-c24b-4a6b-8297-f62cdce1191e
# ╟─40338bfb-3c0d-4990-b84b-0ca6d42373c0
# ╠═0345624c-b9f8-4a13-9371-f08f428523ce
# ╠═7c4713b0-e397-4f06-9abf-77e73480a671
# ╠═3d1a9815-e856-4c2d-a7ab-2ff2d5daadd2
# ╟─a381b387-1ff5-484a-ad28-a2df6e529aea

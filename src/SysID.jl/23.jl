### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ a66ab85a-a2a8-11ee-20fa-2bdb7af3d2fd
# ╠═╡ show_logs = false
begin 
 # If you are running this notebook as a stannalone notebook disable this cell.
 import Pkg 
 Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 3e1fc4d8-064d-49f6-8e69-23af74e32bec
using Plots, DSP, Distributions, FFTW

# ╔═╡ 6c2536b2-b96e-4cc5-bc52-559e7165ecc4
md"
# Exercise 23: Generation on multisine with reduced crest factor using random phase generation
!!! purpose
	It's pretty much the title!
"

# ╔═╡ 6ed5acec-c91c-48d4-b3d3-7a5a7fb9780b
begin
	N₂₅ = 256
	N₄₀₀ = 4096
	Ω = 1:10_000

	function random_phase_mutli_sine(N)
		F = length(0:0.1N)
		ϕ_rand = rand(Uniform(0, 2π), N)
		U = zeros(Complex{Float64}, N, 1)
		U[2:F+1] = exp.(2π * im * ϕ_rand[2:F+1])
		u = 2real(ifft(U))
		u / std(u)
	end

	function opt(N)
		u_with_smallest_crest = zeros(N, 1)
		highest_peak_in_u_with_smallest_crest = Inf
		
		for n in Ω
			u = random_phase_mutli_sine(N)
			peak = maximum(abs.(u))
			if  peak < highest_peak_in_u_with_smallest_crest
				# we are only minizing the peaks (negative and positive) because the rms to all signals is `1` due to normalization: `u = u / std(u)`.
		        u_with_smallest_crest = u
		        highest_peak_in_u_with_smallest_crest = peak
			end
		end
		u_with_smallest_crest
	end
end;

# ╔═╡ 33049560-8a16-4322-9156-e0cef97bf17a
begin
	u₂₅ = opt(N₂₅)
	u₄₀₀ = opt(N₄₀₀)
	
	md"### Signal with best crest factor exiting N frequencies"
end

# ╔═╡ 1b955b77-c344-4cd3-9b25-2dda0feb6d01
begin
	p1 = plot(u₂₅;
		xlabel="Time (ms)",
		ylabel="Amplitude",
		ylims=(-4, 4),
		legend=false,
		title="25 exited frequencies"
	)
	
	p2 = histogram(u₂₅;
		xlabel="Amplitude",
		ylabel="N",
		xlims=(-4, 4),
		bins=N₂₅÷2,
		legend=false,
	)
	
	plot(p1, p2; layout=(2, 1))
end

# ╔═╡ 21d49ca0-2c9e-427d-9151-d552e04b812c
begin
	p3 = plot(u₄₀₀;
		xlabel="Time (ms)",
		ylabel="Amplitude",
		ylims=(-4, 4),
		legend=false,
		title="400 exited frequencies"
	)
	
	p4 = histogram(u₄₀₀;
		xlabel="Amplitude",
		ylabel="N",
		xlims=(-4, 4),
		bins=:sqrt,
		legend=false,
	)
	
	plot(p3, p4; layout=(2, 1))
end

# ╔═╡ Cell order:
# ╠═a66ab85a-a2a8-11ee-20fa-2bdb7af3d2fd
# ╠═3e1fc4d8-064d-49f6-8e69-23af74e32bec
# ╟─6c2536b2-b96e-4cc5-bc52-559e7165ecc4
# ╠═6ed5acec-c91c-48d4-b3d3-7a5a7fb9780b
# ╠═33049560-8a16-4322-9156-e0cef97bf17a
# ╟─1b955b77-c344-4cd3-9b25-2dda0feb6d01
# ╟─21d49ca0-2c9e-427d-9151-d552e04b812c

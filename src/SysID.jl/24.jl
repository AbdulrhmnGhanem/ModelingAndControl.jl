### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 4fcf1762-a2be-11ee-1fc3-bbf06040eb20
# ╠═╡ show_logs = false
begin 
 # If you are running this notebook as a stannalone notebook disable this cell.
 import Pkg 
 Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 4ccfbbb9-8038-4a24-8a46-899ee42f329d
using Plots, FFTW, Distributions, DSP, Optimization, OptimizationBBO

# ╔═╡ f5659074-8999-4c7b-abc6-1c6bb91b0bf8
# ╠═╡ disabled = true
#=╠═╡
begin
	function generate_multisine(ϕ, N)
		F = length(0:0.1N)
	    U = zeros(Complex{Float64}, N, 1)
	    U[2:F+1] = exp.(im * ϕ)
	    u = 2real(ifft(U))
	    u / std(u)
	end

	function crest_factor(u)
	    return maximum(u) / rms(u)
	end

	function objective_function(ϕ, N)
	    u = generate_multisine(ϕ, N)
	    return crest_factor(u)
	end

	N₄₀₀ = 4096
	ϕ₀_length₄₀₀ = 410
	ϕ₀₄₀₀ = zeros(ϕ₀_length₄₀₀)
	
	f = OptimizationFunction(objective_function)
	prob = Optimization.OptimizationProblem(f, ϕ₀₄₀₀, N₄₀₀;
		lb=-π*ones(ϕ₀_length₄₀₀, 1),
		ub=π*ones(N₄₀₀, 1),
	)
	sol = solve(prob, BBO_separable_nes();
		maxiters = 100000,
    	maxtime = 1000,
	)
end;
  ╠═╡ =#

# ╔═╡ 481f60c0-3fed-468c-8959-b34accf93fc8
#=╠═╡
plot(generate_multisine(sol, N₄₀₀))
  ╠═╡ =#

# ╔═╡ 538654c2-e0b3-43d3-b7a2-e80cdbbaac06
#=╠═╡
crest_factor(generate_multisine(sol, N₄₀₀))
  ╠═╡ =#

# ╔═╡ 3bb40383-374e-4345-9e70-c78f5166c578
#=╠═╡
histogram(generate_multisine(sol, N₄₀₀);
		xlabel="Amplitude",
		ylabel="N",
		xlims=(-4, 4),
		bins=N₄₀₀÷16,
		legend=false,
	)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═4fcf1762-a2be-11ee-1fc3-bbf06040eb20
# ╠═4ccfbbb9-8038-4a24-8a46-899ee42f329d
# ╠═f5659074-8999-4c7b-abc6-1c6bb91b0bf8
# ╠═481f60c0-3fed-468c-8959-b34accf93fc8
# ╠═538654c2-e0b3-43d3-b7a2-e80cdbbaac06
# ╠═3bb40383-374e-4345-9e70-c78f5166c578

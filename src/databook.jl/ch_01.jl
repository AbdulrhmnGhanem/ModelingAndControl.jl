### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 31f2b14a-d9f2-11ee-01d6-4dd9050e0a83
# ╠═╡ show_logs = false
begin 
 # If you are running this notebook as a stannalone notebook disable this cell.
 import Pkg 
 Pkg.activate(joinpath("..", ".."))
 data_path = joinpath("..", "..", "books", "databook", "DATA")
end;

# ╔═╡ d7a6d9ba-066d-4214-aabc-13c9124389ab
using FileIO, LinearAlgebra, ImageShow, ImageCore, Plots

# ╔═╡ 1c78d592-8b55-4171-b116-9afb4e7f9d25
md"# Chapter 1 - SVD"

# ╔═╡ 59db64bf-fb86-4183-b217-17d23abe8944
md"## Exercise 1.1
Load the image dog.jpg and compute the full SVD. Choose a rank r < m and
confirm that the matrix **U'∗U** is the _r × r_ identity matrix. Now confirm that **UU'∗** is not the identity matrix. Compute the norm of the error between **UU∗** and the n × n identity matrix as the rank r varies from 1 to n and plot the error."

# ╔═╡ 40eb744e-4150-442b-b8cf-9ee4bdc38542
begin
	A = load(joinpath(data_path, "dog.jpg"))
	X = channelview(Gray.(A))
	U, S, V = svd(X; full=true)
	r = 100
	Ũ, S̃, Ṽ = U[:, 1:r], S[1:r], V[:, 1:r]
	M = Ũ * Diagonal(S̃) * Ṽ'
	# the matrix U'∗U is the r × r identity matrix
	@assert (Ũ' * Ũ  ≈ I) && size(Ũ' * Ũ) == (r, r)
	@assert !(Ũ * Ũ' ≈ I)
end;

# ╔═╡ ca108010-a5c2-429f-a884-ccb7240f9a10
simshow(A), simshow(M)

# ╔═╡ dfd82d31-2523-4bf6-9b58-51c48632170b
begin
	m, n = size(A)
	errors_norm = zeros(n)
	ranks = 1:n
	Threads.@threads for r in ranks
		@inbounds Ũ = @view U[:, 1:r]
		@inbounds errors_norm[r] = norm((Ũ * Ũ') - I, 2) / norm(Ũ)
	end
	plot(errors_norm;
		xlabel="rank",
		ylabel="||e||",
	)
end

# ╔═╡ 627d5352-c43a-4afe-8fb1-3d234e8f594e
md"## Exercise 1.2
Load the image dog.jpg and compute the economy SVD. Compute the relative reconstruction error of the truncated SVD in the Frobenius norm as a function of
the rank _r_. Square this error to compute the fraction of missing variance as a function of _r_. You may also decide to plot 1 minus the error or missing variance to visualize the amount of norm or variance captured at a given rank r. Plot these quantities along with the cumulative sum of singular values as a function of _r_. Find  the rank _r_ where the reconstruction captures 99% of the total variance. Compare this with the rank _r_ where the reconstruction captures 99% in the Frobenius norm and with the rank _r_ that captures 99% of the cumulative sum of singular values."

# ╔═╡ a569ad57-9c2a-4f3e-9059-1b039340a9d4
begin
	Û, Ŝ, V̂ = svd(X; full=false)
	function reconstrunction_error(r)
		U, S, V = @views Û[:, 1:r], Ŝ[1:r], V̂[:, 1:r]
		S = Diagonal(S)
		reconstructed = U * S * V'
		norm(X - reconstructed, 2)
	end
	errors = zeros(n)
	Threads.@threads for r in 1:n
		errors[r] = reconstrunction_error(r) 
	end
	singular_values = Ŝ[1:n]
	
	variance = errors .^ 2
	energy(v) = cumsum(v) / sum(v)
	threshold(energy) = energy > 0.90
end

# ╔═╡ 70622515-7b18-4110-a751-60dc4a710f44
begin
	p1 = plot(1:n, errors;
		ylabel="error",
		title="Frobenius norm",
	)
	vline!([findfirst(threshold, energy(errors))])

	p2 = plot(1:n, variance;
		ylabel="error",
		title="Variance",
	)
	vline!([findfirst(threshold, energy(variance))])

	p3 = plot(1:n, singular_values;
		xlabel="rank",
		ylabel="error",
		title="Singular Values",
	)
	vline!([findfirst(threshold, energy(singular_values))])

	plot(p1, p2, p3; layout=(3,1), size=(600, 600))
end

# ╔═╡ Cell order:
# ╠═31f2b14a-d9f2-11ee-01d6-4dd9050e0a83
# ╠═d7a6d9ba-066d-4214-aabc-13c9124389ab
# ╟─1c78d592-8b55-4171-b116-9afb4e7f9d25
# ╟─59db64bf-fb86-4183-b217-17d23abe8944
# ╠═40eb744e-4150-442b-b8cf-9ee4bdc38542
# ╠═ca108010-a5c2-429f-a884-ccb7240f9a10
# ╠═dfd82d31-2523-4bf6-9b58-51c48632170b
# ╟─627d5352-c43a-4afe-8fb1-3d234e8f594e
# ╠═a569ad57-9c2a-4f3e-9059-1b039340a9d4
# ╠═70622515-7b18-4110-a751-60dc4a710f44

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
	U, S, V = svd(X)
	r = 100
	Û, Ŝ, V̂ = U[:, 1:r], S[1:r], V[:, 1:r]
	M = Û * Diagonal(Ŝ) * V̂'
	# the matrix U'∗U is the r × r identity matrix
	@assert (Û' * Û  ≈ I) && size(Û' * Û) == (r, r)
	@assert !(Û * Û' ≈ I)
end;

# ╔═╡ ca108010-a5c2-429f-a884-ccb7240f9a10
simshow(A), simshow(M)

# ╔═╡ dfd82d31-2523-4bf6-9b58-51c48632170b
begin
	n = size(A)[2]
	norms = zeros(n)
	Threads.@threads for r in 1:n
		@inbounds Û = @view U[:, 1:r]
		@inbounds norms[r] = norm((Û * Û') - I, 2) / norm(Û)
	end
	plot(norms;
		xlabel="rank",
		ylabel="||e||",
	)
end

# ╔═╡ Cell order:
# ╠═31f2b14a-d9f2-11ee-01d6-4dd9050e0a83
# ╠═d7a6d9ba-066d-4214-aabc-13c9124389ab
# ╟─1c78d592-8b55-4171-b116-9afb4e7f9d25
# ╟─59db64bf-fb86-4183-b217-17d23abe8944
# ╠═40eb744e-4150-442b-b8cf-9ee4bdc38542
# ╠═ca108010-a5c2-429f-a884-ccb7240f9a10
# ╠═dfd82d31-2523-4bf6-9b58-51c48632170b

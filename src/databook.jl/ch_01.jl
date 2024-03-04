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
using FileIO, LinearAlgebra, Images, ImageShow, ImageCore, Plots, MAT, StatsBase

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
    @assert (Ũ' * Ũ ≈ I) && size(Ũ' * Ũ) == (r, r)
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
        ylabel="Frobenius norm",
        yaxis=:log,
        label="||e||",
    )
    vline!([findfirst(threshold, energy(errors))]; label="90% of energy",)

    p2 = plot(1:n, variance;
        ylabel="Variance",
        yaxis=:log,
        label="variance",
    )
    vline!([findfirst(threshold, energy(variance))]; label="90% of energy",)

    p3 = plot(1:n, singular_values;
        xlabel="rank",
        ylabel="Singular Values",
        yaxis=:log,
        label="signular value",
    )
    vline!([findfirst(threshold, energy(singular_values))]; label="90% of energy",)

    plot(p1, p2, p3; layout=(3, 1), size=(600, 600))
end

# ╔═╡ 4329f763-fb5a-4ec1-b9da-b62817da812d
md"## Exercise 1.3
Load the Yale B image database and compute the economy SVD using a standard svd command. Now compute the SVD with the method of snapshots. Compare the singular value spectra on a log plot. Compare the first 10 left singular vectors using each method (remember to reshape them into the shape of a face). Now compare a few singular vectors farther down the spectrum. Explain your findings.
"

# ╔═╡ 1ef5cbdb-d780-4cd3-9bd8-2ebcea768f90
begin
    faces_dataset = matread(joinpath(data_path, "allFaces.mat"))
    faces = faces_dataset["faces"]
    faces_m = faces_dataset["m"] |> Int
    faces_n = faces_dataset["n"] |> Int
    num_faces = faces_dataset["nfaces"]' .|> Int

    training_faces = faces[:, 1:sum(num_faces[1:36])]
    average_face = mean(training_faces; dims=2)
    a = reshape(repeat(average_face, size(training_faces)[2]), size(training_faces))
    faces_X = training_faces - a
    faces_U, faces_S, faces_Vt = svd(faces_X; full=false)
	Λ, v =  eigen(faces_X' * faces_X)
	Σ = sqrt(Diagonal(sort(abs.(Λ); rev=true)))
	snapshot_U = faces_X * v * inv(Σ)
end;

# ╔═╡ c7f9e533-95ed-4d6e-8708-956132c6ad9a
begin
	plot(faces_S;
	    xlabel="rank",
	    yaxis=:log,
	    label="builtin",
		ls=:dashdot
	)
	plot!(diag(Σ); 
		ls=:dash,
		label="snapshots"
	)
end

# ╔═╡ af4172af-0813-42de-8bd0-e5cf429dae59
md"### First 10 singular values"

# ╔═╡ b40e48db-fcd5-481a-a813-b61abc43ef0e
md"#### snapshot"

# ╔═╡ fa78d85a-23ac-4b76-8e31-cd7267aadd9f
[simshow(reshape(snapshot_U[:, i], faces_n, faces_m)) for i in 1:10]

# ╔═╡ 9674bfdd-ff14-4325-81df-e2de2ee05857
md"#### normal"

# ╔═╡ 57e3a924-0a95-4365-bd1e-0effba947df1
[simshow(reshape(faces_U[:, i], faces_n, faces_m)) for i in 1:10]

# ╔═╡ e3de1114-4f35-430a-9c51-be0a514c3088
md"### Last 10 singular values"

# ╔═╡ b77c56f1-6d81-40dc-9501-51306257af91
md"#### snapshot"

# ╔═╡ a05eb516-4fbf-4e77-a679-325374e8367d
[simshow(reshape(snapshot_U[:, i], faces_n, faces_m)) for i in size(snapshot_U)[2]-10:size(snapshot_U)[2]]

# ╔═╡ df6dab99-df36-405e-b3d1-fca2c534d12b
md"#### normal"

# ╔═╡ 00b71834-546c-4cda-a5e5-47836f5fd436
[simshow(reshape(faces_U[:, i], faces_n, faces_m)) for i in size(faces_U)[2]-10:size(faces_U)[2]]

# ╔═╡ 823f0389-26a2-4c5f-ae31-8b76ff2655c1
md"
!!! comparison
	The further we go into the specturm the clearer the singular vectors looks (the details of the face are more clear), whichi is the opposite to the normal svd apporach.
"

# ╔═╡ f1af5098-248f-4175-a647-2b5b5c6d187b
md"## Exercise 1.4
Generate a random 100×100 matrix, i.e., a matrix whose entries are sampled from a normal distribution. Compute the SVD of this matrix and plot the singular values.
Repeat this 100 times and plot the distribution of singular values in a box-and-whisker plot. Plot the mean and median singular values as a function of r. Now repeat this for different matrix sizes (e.g., 50 × 50, 200 × 200, 500 × 500, 1000 × 1000, etc.)
"

# ╔═╡ efbad08c-40e4-4b47-9433-419a33a76047
md"## Exercise 1.5
Compare the singular value distributions for a 1000 × 1000 uniformly distributed random matrix and a Gaussian random matrix of the same size. Adapt the Gavish–Donoho algorithm to filter uniform noise based on this singular value distribution. Add uniform noise to a data set (either an image or the test low-rank signal) and apply this thresholding algorithm to filter the noise. Vary the magnitude of the noise and compare the results. Is the filtering good or bad?"

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
# ╟─4329f763-fb5a-4ec1-b9da-b62817da812d
# ╠═1ef5cbdb-d780-4cd3-9bd8-2ebcea768f90
# ╟─c7f9e533-95ed-4d6e-8708-956132c6ad9a
# ╟─af4172af-0813-42de-8bd0-e5cf429dae59
# ╟─b40e48db-fcd5-481a-a813-b61abc43ef0e
# ╟─fa78d85a-23ac-4b76-8e31-cd7267aadd9f
# ╟─9674bfdd-ff14-4325-81df-e2de2ee05857
# ╟─57e3a924-0a95-4365-bd1e-0effba947df1
# ╟─e3de1114-4f35-430a-9c51-be0a514c3088
# ╟─b77c56f1-6d81-40dc-9501-51306257af91
# ╟─a05eb516-4fbf-4e77-a679-325374e8367d
# ╟─df6dab99-df36-405e-b3d1-fca2c534d12b
# ╟─00b71834-546c-4cda-a5e5-47836f5fd436
# ╟─823f0389-26a2-4c5f-ae31-8b76ff2655c1
# ╟─f1af5098-248f-4175-a647-2b5b5c6d187b
# ╟─efbad08c-40e4-4b47-9433-419a33a76047

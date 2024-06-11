### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(
                Base.UUID("6e696c72-6542-2067-7265-42206c756150"),
                "AbstractPlutoDingetjes",
            )].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 31f2b14a-d9f2-11ee-01d6-4dd9050e0a83
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
    data_path = joinpath("..", "..", "books", "DATA")
end;

# ╔═╡ d7a6d9ba-066d-4214-aabc-13c9124389ab
# ╠═╡ show_logs = false
using FileIO,
    LinearAlgebra, ImageShow, ImageCore, Plots, MAT, StatsBase, StatsPlots, Compat, PlutoUI

# ╔═╡ 1c78d592-8b55-4171-b116-9afb4e7f9d25
md"# Chapter 1 - SVD"

# ╔═╡ 59db64bf-fb86-4183-b217-17d23abe8944
md"## Exercise 1.1
Load the image dog.jpg and compute the full SVD. Choose a rank r < m and
confirm that the matrix **U'U** is the _r × r_ identity matrix. Now confirm that **UU'** is not the identity matrix. Compute the norm of the error between **UU'** and the n × n identity matrix as the rank r varies from 1 to n and plot the error."

# ╔═╡ 40eb744e-4150-442b-b8cf-9ee4bdc38542
begin
    A = load(joinpath(data_path, "dog.jpg"))
    X = channelview(Gray.(A))
    U, S, V = svd(X; full = true)
    r = 100
    Ũ, S̃, Ṽ = U[:, 1:r], S[1:r], V[:, 1:r]
    M = Ũ * Diagonal(S̃) * Ṽ'
    # the matrix U'∗U is the r × r identity matrix
    @assert (Ũ' * Ũ ≈ I) && size(Ũ' * Ũ) == (r, r)
    @assert !(Ũ * Ũ' ≈ I)
end;

# ╔═╡ dfc64480-e2fb-4569-8ef9-be06eefafe81
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
    plot(errors_norm; xlabel = "rank", ylabel = "||e||", legend = false)
end

# ╔═╡ 627d5352-c43a-4afe-8fb1-3d234e8f594e
md"## Exercise 1.2
Load the image dog.jpg and compute the economy SVD. Compute the relative reconstruction error of the truncated SVD in the Frobenius norm as a function of
the rank _r_. Square this error to compute the fraction of missing variance as a function of _r_. You may also decide to plot 1 minus the error or missing variance to visualize the amount of norm or variance captured at a given rank r. Plot these quantities along with the cumulative sum of singular values as a function of _r_. Find  the rank _r_ where the reconstruction captures 99% of the total variance. Compare this with the rank _r_ where the reconstruction captures 99% in the Frobenius norm and with the rank _r_ that captures 99% of the cumulative sum of singular values."

# ╔═╡ a569ad57-9c2a-4f3e-9059-1b039340a9d4
begin
    Û, Ŝ, V̂ = svd(X; full = false)
    function reconstrunction_error(r)
        U, S, V = @views Û[:, 1:r], Ŝ[1:r], V̂[:, 1:r]
        S = Diagonal(S)
        reconstructed = U * S * V'
        norm(X - reconstructed, 2)
    end
    errors = zeros(n)
    Threads.@threads for r = 1:n
        errors[r] = reconstrunction_error(r)
    end
    singular_values = Ŝ[1:n]

    variance = errors .^ 2
    energy(v) = cumsum(v) / sum(v)
    threshold(energy) = energy > 0.90
end

# ╔═╡ 70622515-7b18-4110-a751-60dc4a710f44
begin
    p1 = plot(1:n, errors; ylabel = "Frobenius norm", yaxis = :log, label = "||e||")
    vline!([findfirst(threshold, energy(errors))]; label = "90% of energy")

    p2 = plot(1:n, variance; ylabel = "Variance", yaxis = :log, label = "variance")
    vline!([findfirst(threshold, energy(variance))]; label = "90% of energy")

    p3 = plot(
        1:n,
        singular_values;
        xlabel = "rank",
        ylabel = "Singular Values",
        yaxis = :log,
        label = "signular value",
    )
    vline!([findfirst(threshold, energy(singular_values))]; label = "90% of energy")

    plot(p1, p2, p3; layout = (3, 1), size = (600, 600))
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
    average_face = mean(training_faces; dims = 2)
    a = reshape(repeat(average_face, size(training_faces)[2]), size(training_faces))
    faces_X = training_faces - a
    # builtin svd
    faces_U, faces_S, faces_Vt = svd(faces_X; full = false)
    # snapshots method
    Λ, v = eigen(faces_X' * faces_X)
    Σ = sqrt(Diagonal(sort(abs.(Λ); rev = true)))
    snapshot_U = faces_X * v * inv(Σ)
end;

# ╔═╡ c7f9e533-95ed-4d6e-8708-956132c6ad9a
begin
    plot(faces_S; xlabel = "rank", yaxis = :log, label = "builtin", ls = :dashdot)
    plot!(diag(Σ); ls = :dash, label = "snapshots")
end

# ╔═╡ af4172af-0813-42de-8bd0-e5cf429dae59
md"### First 10 singular values"

# ╔═╡ b40e48db-fcd5-481a-a813-b61abc43ef0e
md"#### snapshot"

# ╔═╡ fa78d85a-23ac-4b76-8e31-cd7267aadd9f
[simshow(reshape(snapshot_U[:, i], faces_n, faces_m)) for i = 1:10]

# ╔═╡ 9674bfdd-ff14-4325-81df-e2de2ee05857
md"#### normal"

# ╔═╡ 57e3a924-0a95-4365-bd1e-0effba947df1
[simshow(reshape(faces_U[:, i], faces_n, faces_m)) for i = 1:10]

# ╔═╡ e3de1114-4f35-430a-9c51-be0a514c3088
md"### Last 10 singular values"

# ╔═╡ b77c56f1-6d81-40dc-9501-51306257af91
md"#### snapshot"

# ╔═╡ a05eb516-4fbf-4e77-a679-325374e8367d
[
    simshow(reshape(snapshot_U[:, i], faces_n, faces_m)) for
    i = size(snapshot_U)[2]-10:size(snapshot_U)[2]
]

# ╔═╡ df6dab99-df36-405e-b3d1-fca2c534d12b
md"#### normal"

# ╔═╡ 00b71834-546c-4cda-a5e5-47836f5fd436
[
    simshow(reshape(faces_U[:, i], faces_n, faces_m)) for
    i = size(faces_U)[2]-10:size(faces_U)[2]
]

# ╔═╡ 823f0389-26a2-4c5f-ae31-8b76ff2655c1
md"
!!! remark
	The further we go into the specturm the clearer the singular vectors looks (the details of the face are more clear), whichi is the opposite to the normal svd apporach.
"

# ╔═╡ f1af5098-248f-4175-a647-2b5b5c6d187b
md"## Exercise 1.4
Generate a random 100×100 matrix, i.e., a matrix whose entries are sampled from a normal distribution. Compute the SVD of this matrix and plot the singular values.
Repeat this 100 times and plot the distribution of singular values in a box-and-whisker plot. Plot the mean and median singular values as a function of r. Now repeat this for different matrix sizes (e.g., 50 × 50, 200 × 200, 500 × 500, 1000 × 1000, etc.)
"

# ╔═╡ 9e92eedc-17d1-477d-b767-818abb59d247
begin
    function solve_four(s)
        singular_values = Array{Float64,2}(undef, (s, s))
        for i = 1:s
            X = randn(s, s)
            U, S, V = svd(X)
            singular_values[i, :] = S
        end
        averaged = mean(singular_values, dims = 1)
        μ = mean.(averaged)
        med = median.(averaged)
        singular_values', μ', med'
    end

    function plot_four(sol)
        p1 = boxplot(sol[1]; legend = false, title = "Box-whisker plot", ylims = (0, 50))
        p2 = plot(sol[1]; legend = false, title = "Singular values")
        p3 = plot(
            [sol[2], sol[3]];
            labels = reshape(["mean", "median"], 1, :),
            ls = reshape([:dot, :dash], 1, :),
            titlle = "",
        )
        plot(p1, p2, p3; layout = (3, 1), size = (600, 600))
    end
end

# ╔═╡ a8717234-cf67-4a12-8e29-c252ea6f05cd
md"### Matrix 50x50"

# ╔═╡ e55c3f39-cdf4-476c-b359-c4ade88aa2c1
plot_four(solve_four(50))

# ╔═╡ a81a9346-3cd9-458f-987b-e3ac01c7faba
md"### Matrix 100x100"

# ╔═╡ beba89ea-4886-4bd1-bd8d-b0b957232bdd
plot_four(solve_four(100))

# ╔═╡ c4329ed7-cb3e-4568-aaae-3f78d94534ee
md"### Matrix 200x200"

# ╔═╡ 9c2c124f-366c-42a4-a1cc-740aa14f7672
plot_four(solve_four(200))

# ╔═╡ 1bf98e6c-51fa-44c2-b9bb-c14b0bd393e7
md"### Matrix 500x500"

# ╔═╡ 7f34f4f4-1311-47c1-a343-dcfb36eb74b4
plot_four(solve_four(500))

# ╔═╡ d3ce0903-84a0-4b1f-bdb2-ee95bbac3f06
md"
!!! remark
	The size of the matrix doesn't change the statstical profile, but it increases its magnitude.  
"

# ╔═╡ efbad08c-40e4-4b47-9433-419a33a76047
md"## Exercise 1.5
Compare the singular value distributions for a 1000 × 1000 uniformly distributed random matrix and a Gaussian random matrix of the same size. Adapt the Gavish–Donoho algorithm to filter uniform noise based on this singular value distribution. Add uniform noise to a data set (either an image or the test low-rank signal) and apply this thresholding algorithm to filter the noise. Vary the magnitude of the noise and compare the results. Is the filtering good or bad?"

# ╔═╡ 21c2f9af-54c8-4e8f-aa7f-aa58a66a44e3
begin
    function solve_five(σ)
        # part 1
        n = 1000
        X1 = rand(n, n)
        X2 = randn(n, n)
        U1, S1, V1 = svd(X1)
        U2, S2, V2 = svd(X2)

        # part 2: adding noise to the dog image and filtering it
        # load the image, convert it gray scale, crop it to make it 1000x1000, add noise
        X =
            joinpath(data_path, "dog.jpg") |>
            load .|>
            Gray |>
            channelview |>
            float |>
            v -> v[1:n, 1:n]
        X_noisy_uniform = X + σ * rand(n, n)
        X_noisy_gaussian = X + σ * randn(n, n)
        U1, S1, Vt1 = svd(X_noisy_uniform)
        U2, S2, Vt2 = svd(X_noisy_gaussian)
        cutoff = (4 / √3) * σ * √n

        r = findfirst(s -> s < cutoff, S1)
        filtered_uniform = U1[:, 1:r] * Diagonal(S1[1:r]) * Vt1[:, 1:r]'
        r = findfirst(s -> s > 0.9, cumsum(S1) / sum(S1))
        truncated_uniform = U1[:, 1:r] * Diagonal(S1[1:r]) * Vt1[:, 1:r]'


        r = findfirst(s -> s < cutoff, S2)
        filtered_gaussian = U2[:, 1:r] * Diagonal(S2[1:r]) * Vt2[:, 1:r]'
        r = findfirst(s -> s > 0.9, cumsum(S2) / sum(S2))
        truncated_gaussian = U2[:, 1:r] * Diagonal(S2[1:r]) * Vt2[:, 1:r]'
        S1,
        S2,
        X_noisy_uniform,
        X_noisy_gaussian,
        filtered_uniform,
        filtered_gaussian,
        truncated_uniform,
        truncated_gaussian
    end

    function plot_five(sol)
        S1,
        S2,
        X_noisy_uniform,
        X_noisy_gaussian,
        filtered_uniform,
        filtered_gaussian,
        truncated_uniform,
        truncated_gaussian = sol
        p1 = plot(Gray.(Diagonal(S1)); axis = ([], false), title = "Uniform")

        p2 = plot(Gray.(Diagonal(S2)); axis = ([], false), title = "Gaussian")
        p3 = plot([S1, S2]; labels = reshape(["Uniform", "Gaussian"], 1, :))
        p4 = plot(Gray.(X_noisy_uniform); axis = ([], false), title = "Noisy uniform")
        p5 = plot(Gray.(filtered_uniform); axis = ([], false), title = "Filtered uniform")
        p6 = plot(
            Gray.(truncated_uniform);
            axis = ([], false),
            title = "Truncated (90%) uniform",
        )
        p7 = plot(Gray.(X_noisy_gaussian); axis = ([], false), title = "Noisy gaussian")
        p8 = plot(Gray.(filtered_gaussian); axis = ([], false), title = "Filtered gaussian")
        p9 = plot(
            Gray.(truncated_gaussian);
            axis = ([], false),
            title = "Truncated (90%) gaussian",
        )
        plot(p1, p2, p3, p4, p5, p6, p7, p8, p9; size = (1000, 1000))
    end
end

# ╔═╡ 2468b1f9-2d31-4388-9368-0b826af18e27
md"Noise: $(@bind σ5 Slider(0.2:0.01:1))"

# ╔═╡ 49858770-312c-428f-b427-305fae1a2369
plot_five(solve_five(σ5))

# ╔═╡ a452e424-2c3c-4606-abae-f841dad6b46a
md"
!!! remark
	The singular values for the uniformly distributed random noise fades slightly faster than the guaissian distributed noise. 
"

# ╔═╡ 93cced67-e8ae-4a76-a38e-9a996bdb25df
md"## Exercise 1.6
This exercise will test the concept of condition number. We will test the
accuracy of solving **Ax = b** when noise is added to **b** for matrices **A** with different condition numbers.

1. To build the two matrices, generate a random **U** ∈ $R^{100×100}$ and **V** ∈ $R^{100×100}$ and then create two **Σ** matrices: the first **Σ** will have singular values spaced logarithmically from 100 to 1, and the second **Σ** will have singular values spaced logarithmically from 100 to $10^{−6}$ . Use these matrices to create two **A** matrices, one with a condition number of 100 and the other with a condition number of 100 million. Now create a random **b** vector, solve for **x** using the two methods, and compare the results. Add a small **ϵ** to **b**, with norm $10^{−6}$ smaller than the norm of **b**. Now solve for **x** using this new **b + ϵ** and compare the results.
2. Now repeat the experiment above with many different noise **ϵ** vectors and compute the distribution of the error; plot this error as a histogram and explain the shape.
3. Repeat the above experiment comparing two **A** matrices with different singular value distributions: the first **Σ** will have values spaced linearly from 100 to 1 and the second **Σ** will have values spaced logarithmically from 100 to 1. Does anything change? Please explain why yes or no.

4. Repeat the above experiment, but now with an **A** matrix that has size 100 × 10. Explain any changes
"

# ╔═╡ bf245148-4bda-4bc2-ba8e-15fd0d416c07
begin
    function solve_six()
        U = randn(100, 100)
        V = randn(100, 100)
        b = randn(100, 1)
        ϵ = (norm(b) - 1e-6) * randn(100, 1)

        # part 1
        Σ₁ = reverse(logrange(1, 100, 100))
        Σ₂ = reverse(logrange(1e-6, 100, 100))
        A₁ = U * Diagonal(Σ₁) * V'
        A₂ = U * Diagonal(Σ₂) * V'
        x₁ = A₁ \ b
        x₂ = A₂ \ b
        xϵ₁ = A₁ \ (b + ϵ)
        xϵ₂ = A₂ \ (b + ϵ)

        # part 2
        errs = zeros(10000)
        for i = 1:10000
            ϵ = (norm(b) - 1e-6) * randn(100, 1)
            xϵᵢ = A₂ \ (b + ϵ)
            errs[i] = norm(xϵᵢ)
        end

        # part 3
        Σ₃ = reverse(LinRange(1, 100, 100))
        Σ₄ = reverse(logrange(1, 100, 100))
        A₃ = U * Diagonal(Σ₃) * V'
        A₄ = U * Diagonal(Σ₄) * V'
        x₃ = A₃ \ b
        x₄ = A₄ \ b
        xϵ₃ = A₃ \ (b + ϵ)
        xϵ₄ = A₄ \ (b + ϵ)

        # part 4
        U = randn(100, 100)
        V = randn(10, 10)
        b = randn(100, 1)
        ϵ = (norm(b) - 1e-6) * randn(100, 1)
        Σ₅ = reverse(logrange(1, 100, 10))
        Σ₆ = reverse(logrange(1e-6, 100, 10))
        A₅ = U * [Diagonal(Σ₅); zeros(90, 10)] * V'
        A₆ = U * [Diagonal(Σ₆); zeros(90, 10)] * V'
        x₅ = A₅ \ b
        x₆ = A₆ \ b
        xϵ₅ = A₅ \ (b + ϵ)
        xϵ₆ = A₆ \ (b + ϵ)

        x₁, x₂, x₃, x₄, x₅, x₆, xϵ₁, xϵ₂, xϵ₃, xϵ₄, xϵ₅, xϵ₆, errs
    end

    function plot_six(sol)
        x₁, x₂, x₃, x₄, x₅, x₆, xϵ₁, xϵ₂, xϵ₃, xϵ₄, xϵ₅, xϵ₆, errs = sol
        p1 = plot(
            [x₁, x₂, xϵ₁, xϵ₂];
            labels = reshape(["x₁", "x₂", "xϵ₁", "xϵ₂"], 1, :),
            title = "(a)",
        )
        p2 = histogram(errs; title = "(b)", label = "errors")
        p3 = plot(
            [x₃, x₄, xϵ₃, xϵ₄];
            labels = reshape(["x₃", "x₄", "xϵ₃", "xϵ₄"], 1, :),
            title = "(c)",
        )
        p4 = plot(
            [x₅, x₆, xϵ₅, xϵ₆];
            labels = reshape(["x₅", "x₆", "xϵ₅", "xϵ₆"], 1, :),
            ls = reshape([:dash, :dashdot, :dash, :dashdot], 1, :),
            title = "(d)",
        )

        plot(p1, p2, p3, p4; layout = (4, 1), size = (800, 800))
    end
end

# ╔═╡ 4476ca35-590b-448d-848e-2f4bcdcb229f


# ╔═╡ 436c13d6-1bc1-442c-bd54-e95586294bd9
plot_six(solve_six())

# ╔═╡ adb59fde-27e8-474c-a8c5-49eca0125d4f
md"## Exercise 1.7
Load the data set for fluid flow past a cylinder (you can either download this
from our book [https://databookuw.com/](https://databookuw.com/) or generate it using the IBPM code on GitHub). Each
column is a flow field that has been reshaped into a vector.

1. Compute the SVD of this data set and plot the singular value spectrum and the leading singular vectors. The **U** matrix contains eigenflow fields and the **ΣV**∗ represents the amplitudes of these eigenflows as the flow evolves in time.

2. Write a code to plot the reconstructed movie for various truncation values _r_. Compute the _r_ value needed to capture 90%, 99%, and 99.9% of the flow energy based on the singular value spectrum (recall that energy is given by the Frobenius norm squared). Plot the movies for each of these truncation values and compare the fidelity. Also compute the squared Frobenius norm of the error between the true matrix **X** and the reconstructed matrix **X̃**, where **X** is the flow field movie.

3. Fix a value r = 10 and compute the truncated SVD. Each column $w_k$ ∈ $R^{10}$ of the matrix $W = \tilde\Sigma \tilde V^*$ represents the mixture of the first 10 eigenflows in the kth column of **X**. Verify this by comparing the kth snapshot of **X** with $\tilde U w_k$

4. Now, build a linear regression model for how the amplitudes w_k evolve in time. This will be a dynamical system: 
$w_{k+1} = Aw_k$
4. Create a matrix **W** with the first 1 through m−1 columns of $\Sigma V^{*}$ and another matrix **W′** with the 2 through m columns of $\Sigma V^{*}$ . We will now try to solve for a best-fit **A** matrix so that
$W' \approx AW$
4. Compute the SVD of **W** and use this to compute the pseudo-inverse of **W** to solve for **A**. Compute the eigenvalues of **A** and plot them in the complex plane.

5. Use this **A** matrix to advance the state $w_{k} = A^{k−1}w_1$ starting from w_1 . Plot the reconstructed flow field using these predicted amplitude vectors and compare with the true values.
"

# ╔═╡ 9962945b-6dc2-4978-bf63-3a6c75e106dd


# ╔═╡ Cell order:
# ╠═31f2b14a-d9f2-11ee-01d6-4dd9050e0a83
# ╠═d7a6d9ba-066d-4214-aabc-13c9124389ab
# ╟─1c78d592-8b55-4171-b116-9afb4e7f9d25
# ╟─59db64bf-fb86-4183-b217-17d23abe8944
# ╠═40eb744e-4150-442b-b8cf-9ee4bdc38542
# ╠═dfc64480-e2fb-4569-8ef9-be06eefafe81
# ╠═dfd82d31-2523-4bf6-9b58-51c48632170b
# ╟─627d5352-c43a-4afe-8fb1-3d234e8f594e
# ╠═a569ad57-9c2a-4f3e-9059-1b039340a9d4
# ╟─70622515-7b18-4110-a751-60dc4a710f44
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
# ╟─9e92eedc-17d1-477d-b767-818abb59d247
# ╟─a8717234-cf67-4a12-8e29-c252ea6f05cd
# ╟─e55c3f39-cdf4-476c-b359-c4ade88aa2c1
# ╟─a81a9346-3cd9-458f-987b-e3ac01c7faba
# ╟─beba89ea-4886-4bd1-bd8d-b0b957232bdd
# ╟─c4329ed7-cb3e-4568-aaae-3f78d94534ee
# ╟─9c2c124f-366c-42a4-a1cc-740aa14f7672
# ╟─1bf98e6c-51fa-44c2-b9bb-c14b0bd393e7
# ╠═7f34f4f4-1311-47c1-a343-dcfb36eb74b4
# ╟─d3ce0903-84a0-4b1f-bdb2-ee95bbac3f06
# ╟─efbad08c-40e4-4b47-9433-419a33a76047
# ╠═21c2f9af-54c8-4e8f-aa7f-aa58a66a44e3
# ╟─2468b1f9-2d31-4388-9368-0b826af18e27
# ╟─49858770-312c-428f-b427-305fae1a2369
# ╟─a452e424-2c3c-4606-abae-f841dad6b46a
# ╟─93cced67-e8ae-4a76-a38e-9a996bdb25df
# ╠═bf245148-4bda-4bc2-ba8e-15fd0d416c07
# ╠═4476ca35-590b-448d-848e-2f4bcdcb229f
# ╟─436c13d6-1bc1-442c-bd54-e95586294bd9
# ╟─adb59fde-27e8-474c-a8c5-49eca0125d4f
# ╠═9962945b-6dc2-4978-bf63-3a6c75e106dd

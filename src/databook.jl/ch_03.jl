### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 03fc6b06-dd40-11ee-2447-c16ee3267487
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
    data_path = joinpath("..", "..", "books", "DATA")
end;

# ╔═╡ 72a551d1-1282-4cce-9034-2dee7fabd6d7
using FileIO, ImageShow, ImageCore, FFTW, LinearAlgebra, Compat, Plots, Distributions, Distances, JuMP, Optim, StatsPlots, Distributions, MAT, SparseArrays

# ╔═╡ 00f5bf70-27cd-47cb-af66-1963dea45fe5
md"# Chapter 3 - Sparsity and Compressed Sensing"

# ╔═╡ 2098bbc5-0f7d-4c83-85d1-88a55cb878e2
md"## Exercise 3.1

Load the image dog.jpg and convert to grayscale. We will repeat Exercise 2.1,
using the FFT to compress the image at different compression ratios. However, now, we
will compare the error versus compression ratio for the image downsampled at different
resolutions. Compare the original image (2000 × 1500) and downsampled copies of the
following resolutions: 1000 × 750, 500 × 375, 200 × 150, and 100 × 75. Plot the error
versus compression ratio for each image resolution on the same plot. Explain the observed
trends.
"

# ╔═╡ 1dd2f677-ce86-4292-940f-581a20ee258b
function solve_one(resolutions, ratios)
    A = joinpath(data_path, "dog.jpg") |> load .|> Gray |> channelview |> float

    errors = zeros(length(resolutions), length(ratios))
    for (i, resolution) in enumerate(resolutions)
        x_scale, y_scale = size(A) .÷ resolution
        A = A[1:x_scale:end, 1:y_scale:end]
        B = fft(A)
        B_abs = abs.(B)
        B_abs_sorted = sort(B_abs[:])

        for (j, ratio) in enumerate(ratios)
            threshold = B_abs_sorted[floor((1 - ratio) * length(B_abs_sorted))|>Int]
            mask = B_abs .> threshold
            coeffecients = B .* mask
            compressed = real.(ifft(coeffecients))
            errors[i, j] = norm(A - compressed)
        end
    end

    plot(ratios, errors';
        xlabel="Compression ratio",
        ylabel="error",
        labels=reshape(["$(x)x$y" for (x, y) in resolutions], 1, :),
        legendtitle="Resolution",
    )
end

# ╔═╡ b1274594-558b-47d8-8c66-b5fa7f6076df
solve_one([[2000, 1500], [1000, 750], [500, 375], [200, 150], [100, 75]], logrange(0.001, 0.9, 10))

# ╔═╡ b55f6990-88fe-4fd6-8827-c39e1adf6eb4
md"
!!! remark
	The error decreases as the resolution decreases.
"

# ╔═╡ d4c32227-7f0b-4962-ba92-78df8978cd92
md"## Exercise 3.2

This example will explore geometry and sampling probabilities in high-dimensional spaces. Consider a two-dimensional square dart board with length L = 2 on both sides and a circle of radius R = 1 in the middle. Write a program to throw 10000 darts by generating a uniform random x and y position on the square. Compute the radius for each point and compute what fraction land inside the circle (i.e., how many have radius < 1). Is this consistent with your expectation based on the area of the circle and the square? Repeat this experiment, throwing 10000 darts randomly (sampled from a uniform distribution) on an N-dimensional cube (length L = 2) with an N-dimensional sphere inside (radius R = 1), for N = 2 through N = 10. For a given N, what fraction of the points land inside the sphere. Plot this fraction versus N. Also compute the histogram of the radii of the randomly sampled points for each N and plot these. What trends do you notice in the data?
"

# ╔═╡ 902809dc-675d-4e0b-80ac-d99d88ecc27d
begin
    function solve_two()
        L = 2
        R = 1
        Ns = 2:10
        num_darts = 10_000
        dist = Uniform(0, L)

        fractions = []
        radiis = []
        for n in Ns
            darts = Vector{Matrix{Float64}}(undef, num_darts)
            for i in 1:num_darts
                darts[i] = rand(dist, 1, n)
            end
            distance(p) = euclidean(vec(p), ones(n))
            radii = map(distance, darts)
            fractions_inside = length(filter(d -> d .< R, radii)) / num_darts
            push!(radiis, radii)
            push!(fractions, fractions_inside)
        end
        area_ratio = π * R^2 / L^2
        Ns, fractions, radiis, area_ratio
    end

    function plot_two(Ns, fractions_inside, radiis)
        p1 = histogram(radiis;
            xlabel="Radius",
            ylabel="Count",
            legendtitle="N",
            labels=reshape(2:10, 1, :),
            fillalpha=0.7,
        )
        p2 = plot(Ns, fractions_inside;
            xlabel="N",
            ylabel="Fraction of darts inside",
            legends=false,
        )
        plot(p1, p2;
            layout=(2, 1),
            size=(800, 800),
        )
    end
end

# ╔═╡ 61953b78-378e-4127-86f6-3ad4afc2c4bc
begin
    Ns, fractions_inside, radiis, area_ratio = solve_two()
    md"
   The area ratio = $(round(area_ratio; digits=4))
   	
   The fraction of darts inside = $(round(fractions_inside[1]; digits=4))
   "
end

# ╔═╡ c6e556d5-27c8-41b1-a60d-c167b6bbf97d
plot_two(Ns, fractions_inside, radiis)

# ╔═╡ 78bafda2-6465-4df9-b652-6e35f956548d
begin
    function build_and_test_dftmtx()
        """A 2d discrete fourier transform matrix"""
        function dftmtx(N)
            fft(Matrix{ComplexF64}(I, N, N), 2)
        end
        x = [1, 2, 3, 0, 3]
        @assert all(dftmtx(5) * x .≈ fft(x))
        dftmtx
    end
    dftmtx = build_and_test_dftmtx()
    heatmap(real.(dftmtx(200)), c=:inferno, color=:black, title="DFT Matrix Ψ")
end

# ╔═╡ 38bbd6a1-c5b6-4e13-8ecd-f986b22bd4ec
md"## Exercise 3.3

This exercise will explore the relationship between the sparsity _K_, the signal
size _n_, and the number of samples _p_ in compressed sensing.

1. For _n_ = 1000 and _K_ = 5, create a _K_-sparse vector **s** of Fourier coefficients in a Fourier basis **Ψ**. For each _p_ from 1 to 100, create a Gaussian random sampling matrix $\textbf{C} \in \mathcal{R}^{p \times n}$ to create a measurement vector y = **CΨs**. Use compressed sensing based on this measurement to estimate $\hat s$. For each p, repeat this with at least 10 realizations of the random measurement matrix **C**. Plot the average relative error of $||\hat s - s||_2 / ||s||$ versus p; it may be helpful to visualize the errors with a box-and- whisker plot. Explain the trends. Also plot the average $\mathcal{l}_1$ and $\mathcal{l}_0$ error versus p.

2. Repeat the above experiment for K = 1 through K = 20. What changes?

3. Now repeat the above experiment for K = 5, varying the signal size using n = 100, n = 500, n = 1000, n = 2000, and n = 5000.
"

# ╔═╡ 7265beb5-13af-4afa-ac7d-604029b9ff85
function solve_three(distribution)
    n = 1000
    k = 5
    Ψ = dftmtx(n)
    s = zeros(ComplexF64, n)
    non_zero_indices = sample(1:n, k, replace=false)
    s[non_zero_indices] .= 1

    num_of_selection = 100
    num_of_realizations = 10
    step = 10
    
	errs = zeros(num_of_selection ÷ step, num_of_realizations)
	l₀= similar(errs)
	l₁ = similar(errs)
	
    for (j, p) in collect(enumerate(1:step:num_of_selection))
        Threads.@threads for i in 1:num_of_realizations
            C = rand(distribution, p, n)
            y = C * Ψ * s
            model = Model(Optim.Optimizer)
            @variable(model, s_L1[1:n])
            @objective(model, Min, sum(abs.(s_L1)))
            @constraint(model, C * Ψ * s_L1 .== y)
            optimize!(model)
            ŝ = value.(s_L1)
			l₀[j, i] = norm(ŝ, 0)
			l₁[j, i] = norm(ŝ, 1)
            errs[j, i] = norm(ŝ - s) / norm(s, 1)
        end
    end
	
	errs2 = zeros(20, num_of_realizations)
	l₀2= similar(errs2)
	l₁2 = similar(errs2)
	
	for (j, k) in collect(enumerate(1:20))
		Threads.@threads for i in 1:num_of_realizations
			s = zeros(ComplexF64, n)
    		non_zero_indices = sample(1:n, k, replace=false)
    		s[non_zero_indices] .= 1

			C = rand(distribution, 10, n)
			y = C * Ψ * s

			model = Model(Optim.Optimizer)
            @variable(model, s_L1[1:n])
            @objective(model, Min, sum(abs.(s_L1)))
            @constraint(model, C * Ψ * s_L1 .== y)
            optimize!(model)
            ŝ = value.(s_L1)
			l₀2[j, i] = norm(ŝ, 0)
			l₁2[j, i] = norm(ŝ, 1)
            errs2[j, i] = norm(ŝ - s) / norm(s, 1)
		end
	end

	ns = [100, 500, 1000, 2000, 5000]
	errs3 = zeros(length(ns), num_of_realizations)
	l₀3= similar(errs3)
	l₁3 = similar(errs3)
	
	for (j, n) in collect(enumerate(ns))
		k = 5
		Ψ = dftmtx(n)
		s = zeros(ComplexF64, n)
		non_zero_indices = sample(1:n, k, replace=false)
		s[non_zero_indices] .= 1

		Threads.@threads for i in 1:num_of_realizations
			C = rand(distribution, 10, n)
			y = C * Ψ * s

			model = Model(Optim.Optimizer)
            @variable(model, s_L1[1:n])
            @objective(model, Min, sum(abs.(s_L1)))
            @constraint(model, C * Ψ * s_L1 .== y)
            optimize!(model)
            ŝ = value.(s_L1)
			l₀3[j, i] = norm(ŝ, 0)
			l₁3[j, i] = norm(ŝ, 1)
            errs3[j, i] = norm(ŝ - s) / norm(s, 1)
		end
	end
    errs, l₀, l₁, errs2, l₀2, l₁2, errs3, l₀3, l₁3, ns
end

# ╔═╡ 0077bafb-2beb-4743-a1f1-88b54d345058
function plot_three(sol)
	errs, l₀, l₁, errs2, l₀2, l₁2, errs3, l₀3, l₁3, ns = sol
	p1 = boxplot(1:10:100, errs';
		legend=false,
		title="ẽ",
		xlabel="p",
		outliers=false,
	)
	p2 = boxplot(1:10:100, l₀';
		legend=false,
		title="L₀",
		xlabel="p",
		outliers=false,
	)
		
	p3 = boxplot(1:10:100, l₁';
		legend=false,
		title="L₁",
		xlabel="p",
		outliers=false,
	)
	
	p4 = boxplot(errs2';
		legend=false,
		title="ẽ",
		xlabel="k",
		outliers=false,
	)
	p5 = boxplot(l₀2';
		legend=false,
		title="L₀",
		xlabel="k",
		outliers=false,
	)
	p6 = boxplot(l₁2';
		legend=false,
		title="L₁",
		xlabel="k",
		outliers=false,
	)

	p7 = boxplot(ns', errs3';
		legend=false,
		xlabel="n",
		title="ẽ",
		xrotation=90,
		outliers=false,
	)
	p8 = boxplot(ns', l₀3';
		legend=false,
		title="L₀",
		xlabel="n",
		xrotation=90,
		outliers=false,
	)
	p9 = boxplot(ns', l₁3';
		legend=false,
		xlabel="n",
		title="L₁",
		xrotation=90,
		outliers=false,
	)

	plot(p1, p2, p3, p4, p5, p6, p7, p8, p9;
		layout=(3, 3),
		size=(600, 600)
	)
end

# ╔═╡ 56f955be-2dba-4080-8881-bfc66892a783
plot_three(solve_three(Normal()))

# ╔═╡ 01a25614-364d-473d-87dd-76f33174bb5e
md"## Exercise 3.4

Repeat the above exercise with a uniformly sampled random sample matrix. Also repeat with a Bernoulli random matrix and a matrix that comprises random single pixels. Plot the average relative errors for these different sampling matrices on the same plot (including the plot for Gaussian random sampling). Discuss the trends.
"

# ╔═╡ c2383726-826b-43d9-a670-9f9b84d1a292
md"### unifromly distributed sample matrix"

# ╔═╡ e8f85762-9018-4e3d-b6d1-66e555f9fb5d
plot_three(solve_three(Uniform()))

# ╔═╡ 3fcc4f34-2742-4364-a1b8-8e924423aca1
md"### Bernoulli distributed sample matrix"

# ╔═╡ 1de5daed-9ac7-4780-a276-9ff46ae396d7
plot_three(solve_three(Bernoulli()))

# ╔═╡ c2810938-988d-4558-9d61-271528ef024f
md"## Exercise 3.5
Generate a DFT matrix ψ for n = 512.  We will use this as a basis for compressed sensing, and we will compute the incoherence of this basis and different measurement matrices. For p = 16, create different random measurement matrices **C** given by Gaussian random measurements, Bernoulli random measurements, and random single- pixel measurements. For each matrix, normalize the length of each row to 1. Now, for each measurement matrix type, compute the incoherence μ(**C, ψ**). Repeat this for many random instances of each C matrix type and compare the histogram of incoherence values for each matrix type. Further, compare the histogram of each inner product $\sqrt n \langle c_k, \psi_j \rangle$ for each matrix type. Discuss any trends and the implications for compressed sensing with these measurement matrices. Are there other factors that are relevant for the sensing matrix?

##### the incoherence equation from the book
$\mu(\mathbf{C}, \mathbf{\Psi}) = \sqrt{n} \max_{j,k} \left| \langle c_k, \psi_j \rangle \right|,$
"

# ╔═╡ ab49b66b-9744-4f39-894f-e5b1ab657336
function solve_five()
	normalize_row_length(m) = mapslices(r -> r ./ norm(r), m, dims=2)
	μ(C, Ψ) = √n * maximum(abs.(map(m -> m[1]' * m[2], zip(eachrow(C), eachcol(Ψ)))))
	
	n = 512
	p = 16
	# normalize each col length in Ψ to 1
	Ψ = dftmtx(n) |> m -> mapslices(r -> r ./ norm(r), m, dims=1)
	num_of_trials = 10_000 
	incoherence = zeros(num_of_trials, 3)
	for i in 1:num_of_trials
		C₁ = randn(p, n) |> normalize_row_length # gaussian
		C₂ = rand(Bernoulli(), p, n) |> normalize_row_length  # bernoulli
		C₃ = zeros(p, n)
		C₃[2, 2] = 1  # single pixel
		
		incoherence[i,1] = μ(C₁, Ψ)
		incoherence[i,2] = μ(C₂, Ψ)
		incoherence[i,3] = μ(C₃, Ψ)
	end
	incoherence
end

# ╔═╡ c2700784-253f-462d-a0db-2257b54e2202
function plot_five(sol)
	incoherence = sol

	p1 = histogram(incoherence[:, 1];
		title="Normal",
		xlabel="μ",
		ylabel="N",
		
	)
	p2 = histogram(incoherence[:, 2];
		title="Bernoulli",
		xlabel="μ",
		ylabel="N",
	)
	p3 = histogram(incoherence[:, 3];
		title="Single pixel",
		xlabel="μ",
		ylabel="N",
	)

	plot(p1, p2, p3;
		layout=(3, 1),
		legend=false,
	)
end

# ╔═╡ 19dc20eb-1e63-429e-aef0-bc29fedfed3c
plot_five(solve_five())

# ╔═╡ 5d9b91f0-c6da-43f4-af6f-95d7df6db2de
md"
!!! remark
	- Choosing C with a normal distribution results in low coherence that has a distribtuion close to a beta distribution.
	
	- Choosing C with a bernoulli distribution lead to much higher incoherence with what looks like a mixture distribution.

	- As expected choosing one pixel is the worst (high incoherence is desired for sparse representation).
"

# ╔═╡ 63e19dbb-5845-412f-91ad-26a8f89af46c
md"## Exercise 3.6
This exercise will explore sparse representation from Section 3.6 to estimate a fluid flow field, following Callaham et al. [146]. Load the cylinder flow data set. Coarsen each flow field by a factor of 20 in each direction using imresize, and build a library of these coarsened measurements (i.e., a matrix, where each column contains these downsampled measurements). Plot a movie of the flow field in these new coordinates. Now, pick a column of the full flow field matrix and add Gaussian random noise to this field. Downsample the noisy field by a factor of 20 and use SRC to find the closest downsampled library element. Then use this column of the full flow field library as your reconstructed estimate. Try this approach with different levels of noise added to the original flow field. See how much noise is required before the method breaks. Try different approaches to creating a low-dimensional representation of the image (i.e., instead of downsampling, you can measure the flow field in a small 5 × 10 patch and use this as the low-dimensional feature for SRC).
"

# ╔═╡ 48e9967c-657b-4ed8-8b02-04a058114b6c
function plot_six()
	V = matread(joinpath(data_path, "CYLINDER_ALL.mat"))["VORTALL"]
	vortmin = -5
    vortmax = 5
	V[V .> vortmax] .= vortmax
    V[V .< vortmin] .= vortmin
	
	grid_dims = (449, 199)

	p = plot()
	t = (1:100) / 100 .* 2 .* π
    x = 49 .+ 25 .* sin.(t)
    y = 99 .+ 25 .* cos.(t)
    
	@gif for i in 1:size(V, 2)
		heatmap!(p, reshape(V[:, i], grid_dims);
			color=:viridis,
			clim=(vortmin, vortmax),
			legend=false,
			ylims=(0, 200),
			xlims=(0, 200)
		)
		plot!(p, x, y, linecolor=:yellow, linewidth=4)
	end
end

# ╔═╡ 43a50a1b-13d3-4e9f-a892-acafe76c01ae
# ╠═╡ disabled = true
#=╠═╡
plot_six()
  ╠═╡ =#

# ╔═╡ 2a95f1b3-fe43-4cc3-942d-1ea443b6515d
md"## Exercise 3.7

This exercise will explore RPCA from Section 3.7 for robust flow field anal-
ysis, following Scherl et al.

1. Load the cylinder flow data set. Compute the SVD as in Exercise 1.7 and plot the movie of the flow. Also plot the singular values and first 10 singular vectors.

2. Now, contaminate a random 10% of the entries of the data matrix with salt-and- pepper noise. The contaminated points should not be the same for each column, and the salt-and-pepper noise should be ±5η, where η is the standard deviation of the entire data set. Compute the SVD of the contaminated matrix and plot the movie of the flow along with the singular values and first 10 singular vectors.

3. Clean the contaminated data set by applying RPCA and keeping the low-rank portion **L**. Again, compute the SVD of the decontaminated matrix L and plot the movie of  the flow along with the singular values and first 10 singular vectors. Compare these with the results from the original clean and contaminated data sets.

4. Try to clean the data by applying the Gavish–Donoho threshold to the data matrix contaminated with salt-and-pepper noise. Does this work? Explain why or why not.
"

# ╔═╡ 1c4c379f-edeb-472c-8ffe-d78b58caa222


# ╔═╡ 705cd927-f900-45ff-9471-981495ad5a05
md"## Exercise 3.8

This exercise will explore the sparse sensor selection approach based on QR
from Section 3.8.

1. Load the Yale B faces data set. Randomly choose one person to omit from the data matrix and compute the SVD of the remaining data. Compute the QR sensor locations for p = 100 using the first r = 100 modes of this SVD basis $\hat U$. Use these sensor locations to reconstruct the images of the person that was left out of the matrix for the SVD. Compare the reconstruction error using these QR sensor locations with reconstruction using p = 100 randomly chosen sensors, as in Fig. 3.22.

2. Now, repeat this experiment 36 times, each time omitting a different person from the data before computing the SVD, and use the sensor locations to reconstruct the images of the omitted person. This will provide enough reconstruction errors on which to perform statistics. For each experiment, also compute the reconstruction error using 36 different configurations of p = 100 random sensors. Plot the histograms of the error for the QR and random sensors, and discuss.

3. Finally, repeat the above experiments for different sensor number p = 10 through p = 200 in increments of 10. Plot the error distributions versus p for QR and random sensor configurations. Because each value of p corresponds to many reconstruction errors, it would be best to plot this as a box-and-whisker plot or as a violin plot.
"

# ╔═╡ 7da43eff-c7c3-44b2-9276-1d48f2ad442a
function solve_eight()
	faces_dataset = matread(joinpath(data_path, "allFaces.mat"))
    faces = faces_dataset["faces"] |> m -> mapslices(r -> r ./ norm(r), m, dims=1)
    faces_m = faces_dataset["m"] |> Int
    faces_n = faces_dataset["n"] |> Int
    num_faces = faces_dataset["nfaces"]' .|> Int

	reconstruction_err(test_face, reconstructed) = norm(test_face - reconstructed, 2)

	ps = 10:10:200
	errs = zeros(length(ps), length(num_faces), 2)
	
	cursor = 1
	for (i, num) in collect(enumerate(num_faces[1:end-1]))
		cursor += num
		v1 = @view faces[:, 1:sum(num_faces[1:i-1])]
		v2 = @view faces[:, sum(num_faces[1:i])+1:end]
		training_faces = hcat(v1, v2)
		test_face = @view faces[:, cursor+1]

	    average_face = mean(training_faces; dims=2)
	    a = reshape(repeat(average_face, size(training_faces)[2]), size(training_faces))
	    X = training_faces - a
		U, S, Vt = svd(X; full=false)
		Threads.@threads for (j, p) in collect(enumerate(10:10:200))
			ψᵣ =  Ũ = @view U[:, 1:p]
			Q, R, pivot = qr(ψᵣ', ColumnNorm())
		
			C = zero(ψᵣ')
			for k in 1:p
				@inbounds C[k, pivot[k]] = 1
			end
			@inbounds Θ = C * ψᵣ
		
			@inbounds y_qr_sensor_placement = @view test_face[pivot[1:p], :]
			@inbounds y_random_sensor_placement = @view test_face[rand(1:length(test_face), p), :]
			a_qr_sensor_placement = Θ \ y_qr_sensor_placement
			a_random_sensor_placement = Θ \ y_random_sensor_placement
			@inbounds face_recon_qr_sensor_placement = ψᵣ * a_qr_sensor_placement
			@inbounds face_recon_random_sensor_placement = ψᵣ * a_random_sensor_placement
		
			errs[j, i, 1] = reconstruction_err(test_face, face_recon_qr_sensor_placement)
			errs[j, i, 2] = reconstruction_err(test_face, face_recon_random_sensor_placement)
		end
	end
	errs, ps
end

# ╔═╡ cab7eecb-6e79-47be-a58b-6f44420bb191
function plot_eight(sol)
	errs, ps = sol
	p1 = plot(;
		xlabel="p",
		ylabel="reconstruction error",
		legend=false,
		ylims=(0.1, 1.8),
		title="QR sensors"
	)
	xticks!(([5.0, 10.0, 15.0, 20], ["50", "100", "150", "200"]))

	p2 = plot(;
		xlabel="p",
		legend=false,
		ylims=(0.1, 1.8),
		title="Random sensors"
	)
	xticks!(([5.0, 10.0, 15.0, 20], ["50", "100", "150", "200"]))
	
	for i in 1:length(10:10:200)
		boxplot!(p1, errs[i,:,1];
			outliers=false,		
		)
		boxplot!(p2, errs[i, :, 2];
			outliers=false,
		)
	end
	plot(p1, p2)
end

# ╔═╡ 1be97da5-57f0-4902-94e2-dcc8f6904c84
md"!!! remark
	Increasing the number of QR sensors significantly decreases the reconstruction error. Compared to the randomly placed sensors which becomes asymptotic around 50 sensors.
"

# ╔═╡ 1b782cf5-84dd-4a7a-a45a-48d62296e840
plot_eight(solve_eight())

# ╔═╡ 6af0479e-ae0c-4c30-952f-a13702613fc0
md"## Exercise 3.9
In the exercise above, for p = 100, compare the reconstruction results using
the p = 100 QR sensors to reconstruct in the first r = 100 modes, versus using the same
sensors to reconstruct in the first r = 90 SVD modes. Is one more accurate than the
other? Compare the condition number of the 100 × 100 and 100 × 90 matrices obtained
by sampling the p rows of the r = 100 and r = 90 columns of $\hat U$ from the SVD.
"

# ╔═╡ 79a94dcb-2948-4ea0-adf2-8a8b108f7cf9
function solve_nine()
	faces_dataset = matread(joinpath(data_path, "allFaces.mat"))
    faces = faces_dataset["faces"] |> m -> mapslices(r -> r ./ norm(r), m, dims=1)
    faces_m = faces_dataset["m"] |> Int
    faces_n = faces_dataset["n"] |> Int
    num_faces = faces_dataset["nfaces"]' .|> Int

	reconstruction_err(test_face, reconstructed) = norm(test_face - reconstructed, 2)

	errs = zeros(length(num_faces), 4)
	
	cursor = 1
	for (i, num) in collect(enumerate(num_faces[1:end-1]))
		cursor += num
		v1 = @view faces[:, 1:sum(num_faces[1:i-1])]
		v2 = @view faces[:, sum(num_faces[1:i])+1:end]
		training_faces = hcat(v1, v2)
		test_face = @view faces[:, cursor+1]

	    average_face = mean(training_faces; dims=2)
	    a = reshape(repeat(average_face, size(training_faces)[2]), size(training_faces))
	    X = training_faces - a
		U, S, Vt = svd(X; full=false)
		p = 100
		r = 90
		ψᵣ =  Ũ = @view U[:, 1:p]
		Q, R, pivot = qr(ψᵣ', ColumnNorm())
		C = zero(ψᵣ')
		for k in 1:p
			@inbounds C[k, pivot[k]] = 1
		end
		@inbounds Θ = C * ψᵣ
		Θ̃ = @view Θ[:, 1:r]

		@inbounds y_qr_sensor_placement = @view test_face[pivot[1:p], :]
		@inbounds y_random_sensor_placement = @view test_face[rand(1:length(test_face), p), :]

		a_qr_sensor_placement = Θ \ y_qr_sensor_placement
		a_random_sensor_placement = Θ \ y_random_sensor_placement
		
		ã_qr_sensor_placement = Θ̃ \ y_qr_sensor_placement
		ã_random_sensor_placement = Θ̃ \ y_random_sensor_placement
		
		@inbounds face_recon_qr_sensor_placement = ψᵣ * a_qr_sensor_placement
		@inbounds face_recon_random_sensor_placement = ψᵣ * a_random_sensor_placement
		
		@inbounds face_recon_qr_sensor_placement2 = ψᵣ[:,1:r] * ã_qr_sensor_placement
		@inbounds face_recon_random_sensor_placement2 = ψᵣ[:, 1:r] * ã_random_sensor_placement
	
		errs[i, 1] = reconstruction_err(test_face, face_recon_qr_sensor_placement)
		errs[i, 2] = reconstruction_err(test_face, face_recon_random_sensor_placement)

		errs[i, 3] = reconstruction_err(test_face, face_recon_qr_sensor_placement2)
		errs[i, 4] = reconstruction_err(test_face, face_recon_random_sensor_placement2)
	end
	errs
end

# ╔═╡ 35ab1fcd-f12c-411e-9229-04a3533f4a1f
function plot_nine(sol)
	groupedbar(sol;
	bar_position = :stack,
	xlabel="Face index",
	ylabel="reconstruction error",
	labels=reshape(["QR (r = 100)", "Random (r = 100)", "QR (r = 90)", "Random (r = 90)"], 1, :)
)
end

# ╔═╡ 47e0accf-5b04-4a7d-af11-0b54b6b44bc7
plot_nine(solve_nine())

# ╔═╡ c83b207c-d532-4dbc-befc-4c14cb5996dd
md"!!! remark
	Re-constructing less modes than the sensor leads to a better error as expected (see page 126).
"

# ╔═╡ Cell order:
# ╠═03fc6b06-dd40-11ee-2447-c16ee3267487
# ╠═72a551d1-1282-4cce-9034-2dee7fabd6d7
# ╟─00f5bf70-27cd-47cb-af66-1963dea45fe5
# ╟─2098bbc5-0f7d-4c83-85d1-88a55cb878e2
# ╠═1dd2f677-ce86-4292-940f-581a20ee258b
# ╟─b1274594-558b-47d8-8c66-b5fa7f6076df
# ╟─b55f6990-88fe-4fd6-8827-c39e1adf6eb4
# ╟─d4c32227-7f0b-4962-ba92-78df8978cd92
# ╠═902809dc-675d-4e0b-80ac-d99d88ecc27d
# ╟─61953b78-378e-4127-86f6-3ad4afc2c4bc
# ╟─c6e556d5-27c8-41b1-a60d-c167b6bbf97d
# ╠═78bafda2-6465-4df9-b652-6e35f956548d
# ╟─38bbd6a1-c5b6-4e13-8ecd-f986b22bd4ec
# ╠═7265beb5-13af-4afa-ac7d-604029b9ff85
# ╠═0077bafb-2beb-4743-a1f1-88b54d345058
# ╠═56f955be-2dba-4080-8881-bfc66892a783
# ╟─01a25614-364d-473d-87dd-76f33174bb5e
# ╟─c2383726-826b-43d9-a670-9f9b84d1a292
# ╠═e8f85762-9018-4e3d-b6d1-66e555f9fb5d
# ╟─3fcc4f34-2742-4364-a1b8-8e924423aca1
# ╠═1de5daed-9ac7-4780-a276-9ff46ae396d7
# ╟─c2810938-988d-4558-9d61-271528ef024f
# ╠═ab49b66b-9744-4f39-894f-e5b1ab657336
# ╠═c2700784-253f-462d-a0db-2257b54e2202
# ╟─19dc20eb-1e63-429e-aef0-bc29fedfed3c
# ╟─5d9b91f0-c6da-43f4-af6f-95d7df6db2de
# ╟─63e19dbb-5845-412f-91ad-26a8f89af46c
# ╠═48e9967c-657b-4ed8-8b02-04a058114b6c
# ╠═43a50a1b-13d3-4e9f-a892-acafe76c01ae
# ╟─2a95f1b3-fe43-4cc3-942d-1ea443b6515d
# ╠═1c4c379f-edeb-472c-8ffe-d78b58caa222
# ╟─705cd927-f900-45ff-9471-981495ad5a05
# ╠═7da43eff-c7c3-44b2-9276-1d48f2ad442a
# ╠═cab7eecb-6e79-47be-a58b-6f44420bb191
# ╟─1be97da5-57f0-4902-94e2-dcc8f6904c84
# ╟─1b782cf5-84dd-4a7a-a45a-48d62296e840
# ╟─6af0479e-ae0c-4c30-952f-a13702613fc0
# ╠═79a94dcb-2948-4ea0-adf2-8a8b108f7cf9
# ╠═35ab1fcd-f12c-411e-9229-04a3533f4a1f
# ╟─47e0accf-5b04-4a7d-af11-0b54b6b44bc7
# ╟─c83b207c-d532-4dbc-befc-4c14cb5996dd

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
    data_path = joinpath("..", "..", "books", "databook", "DATA")
end;

# ╔═╡ 72a551d1-1282-4cce-9034-2dee7fabd6d7
using FileIO, LinearAlgebra, Plots, PlutoUI

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


# ╔═╡ d4c32227-7f0b-4962-ba92-78df8978cd92
md"## Exercise 3.2

This example will explore geometry and sampling probabilities in high-dimensional spaces. Consider a two-dimensional square dart board with length L = 2 on both sides and a circle of radius R = 1 in the middle. Write a program to throw 10000 darts by generating a uniform random x and y position on the square. Compute the radius for each point and compute what fraction land inside the circle (i.e., how many have radius < 1). Is this consistent with your expectation based on the area of the circle and the square? Repeat this experiment, throwing 10 000 darts randomly (sampled from a uniform distribution) on an N-dimensional cube (length L = 2) with an N-dimensional sphere inside
(radius R = 1), for N = 2 through N = 10. For a given N, what fraction of the points land inside the sphere. Plot this fraction versus N. Also compute the histogram of the radii of the randomly sampled points for each N and plot these. What trends do you notice in the data?
"

# ╔═╡ 902809dc-675d-4e0b-80ac-d99d88ecc27d


# ╔═╡ 38bbd6a1-c5b6-4e13-8ecd-f986b22bd4ec
md"## Exercise 3.3

This exercise will explore the relationship between the sparsity K, the signal
size n, and the number of samples p in compressed sensing.

1. For n = 1000 and _K_ = 5, create a _K_-sparse vector **s** of Fourier coefficients in a Fourier basis Ψ. For each p from 1 to 100, create a Gaussian random sampling matrix C ∈ Rp×n to create a measurement vector y = **CΨs**. Use compressed sensing based on this measurement to estimate $\hat s$. For each p, repeat this with at least 10 realizations of the random measurement matrix C. Plot the average relative error of $||\hat s - s||_2 / ||s||$ versus p; it may be helpful to visualize the errors with a box-and- whisker plot. Explain the trends. Also plot the average $l_1$ and $l_0$ error versus p.

2. Repeat the above experiment for K = 1 through K = 20. What changes?

3. Now repeat the above experiment for K = 5, varying the signal size using n = 100, n = 500, n = 1000, n = 2000, and n = 5000.
"

# ╔═╡ 7265beb5-13af-4afa-ac7d-604029b9ff85


# ╔═╡ 01a25614-364d-473d-87dd-76f33174bb5e
md"## Exercise 3.4

Repeat the above exercise with a uniformly sampled random sample matrix. Also repeat with a Bernoulli random matrix and a matrix that comprises random single pixels. Plot the average relative errors for these different sampling matrices on the same plot (including the plot for Gaussian random sampling). Discuss the trends.
"

# ╔═╡ 77d4e06b-5dfd-4404-b1a5-195c5ba34b35


# ╔═╡ c2810938-988d-4558-9d61-271528ef024f
md"## Exercise 3.5
Generate a DFT matrix ψ for n = 512.  We will use this as a basis for compressed sensing, and we will compute the incoherence of this basis and different measurement matrices. For p = 16, create different random measurement matrices **C** given by Gaussian random measurements, Bernoulli random measurements, and random single- pixel measurements. For each matrix, normalize the length of each row to 1. Now, for each measurement matrix type, compute the incoherence μ(**C, ψ**). Repeat this for many random instances of each C matrix type and compare the histogram of incoherence values for each matrix type. Further, compare the histogram of each inner product $\sqrt n (c_k, \psi_j)$ for each matrix type. Discuss any trends and the implications for compressed sensing with these measurement matrices. Are there other factors that are relevant for the sensing matrix?"

# ╔═╡ ab49b66b-9744-4f39-894f-e5b1ab657336


# ╔═╡ 63e19dbb-5845-412f-91ad-26a8f89af46c
md"## Exercise 3.6
This exercise will explore sparse representation from Section 3.6 to estimate a fluid flow field, following Callaham et al. [146]. Load the cylinder flow data set. Coarsen each flow field by a factor of 20 in each direction using imresize, and build a library of these coarsened measurements (i.e., a matrix, where each column contains these downsampled measurements). Plot a movie of the flow field in these new coordinates. Now, pick a column of the full flow field matrix and add Gaussian random noise to this field. Downsample the noisy field by a factor of 20 and use SRC to find the closest downsampled library element. Then use this column of the full flow field library as your reconstructed estimate. Try this approach with different levels of noise added to the original flow field. See how much noise is required before the method breaks. Try different approaches to creating a low-dimensional representation of the image (i.e., instead of downsampling, you can measure the flow field in a small 5 × 10 patch and use this as the low-dimensional feature for SRC).
"

# ╔═╡ 9ca403ad-95cd-4659-a7af-0e470fc3a034


# ╔═╡ 2a95f1b3-fe43-4cc3-942d-1ea443b6515d
md"## Exercise 3.7

This exercise will explore RPCA from Section 3.7 for robust flow field anal-
ysis, following Scherl et al.

1. Load the cylinder flow data set. Compute the SVD as in Exercise 1.7 and plot the movie of the flow. Also plot the singular values and first 10 singular vectors.

2. Now, contaminate a random 10% of the entries of the data matrix with salt-and- pepper noise. The contaminated points should not be the same for each column, and the salt-and-pepper noise should be ±5η, where η is the standard deviation of the entire data set. Compute the SVD of the contaminated matrix and plot the movie of the flow along with the singular values and first 10 singular vectors.

3. Clean the contaminated data set by applying RPCA and keeping the low-rank portion **L**. Again, compute the SVD of the decontaminated matrix L and plot the movie of  the flow along with the singular values and first 10 singular vectors. Compare these with the results from the original clean and contaminated data sets.

4. Try to clean the data by applying the Gavish–Donoho threshold to the data matrix
contaminated with salt-and-pepper noise. Does this work? Explain why or why not.
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


# ╔═╡ 6af0479e-ae0c-4c30-952f-a13702613fc0
md"## Exercise 3.9
In the exercise above, for p = 100, compare the reconstruction results using
the p = 100 QR sensors to reconstruct in the first r = 100 modes, versus using the same
sensors to reconstruct in the first r = 90 SVD modes. Is one more accurate than the
other? Compare the condition number of the 100 × 100 and 100 × 90 matrices obtained
by sampling the p rows of the r = 100 and r = 90 columns of $\hat U$ from the SVD.
"

# ╔═╡ 79a94dcb-2948-4ea0-adf2-8a8b108f7cf9


# ╔═╡ Cell order:
# ╠═03fc6b06-dd40-11ee-2447-c16ee3267487
# ╠═72a551d1-1282-4cce-9034-2dee7fabd6d7
# ╟─00f5bf70-27cd-47cb-af66-1963dea45fe5
# ╟─2098bbc5-0f7d-4c83-85d1-88a55cb878e2
# ╠═1dd2f677-ce86-4292-940f-581a20ee258b
# ╟─d4c32227-7f0b-4962-ba92-78df8978cd92
# ╠═902809dc-675d-4e0b-80ac-d99d88ecc27d
# ╟─38bbd6a1-c5b6-4e13-8ecd-f986b22bd4ec
# ╠═7265beb5-13af-4afa-ac7d-604029b9ff85
# ╟─01a25614-364d-473d-87dd-76f33174bb5e
# ╠═77d4e06b-5dfd-4404-b1a5-195c5ba34b35
# ╟─c2810938-988d-4558-9d61-271528ef024f
# ╠═ab49b66b-9744-4f39-894f-e5b1ab657336
# ╟─63e19dbb-5845-412f-91ad-26a8f89af46c
# ╠═9ca403ad-95cd-4659-a7af-0e470fc3a034
# ╟─2a95f1b3-fe43-4cc3-942d-1ea443b6515d
# ╠═1c4c379f-edeb-472c-8ffe-d78b58caa222
# ╟─705cd927-f900-45ff-9471-981495ad5a05
# ╠═7da43eff-c7c3-44b2-9276-1d48f2ad442a
# ╟─6af0479e-ae0c-4c30-952f-a13702613fc0
# ╠═79a94dcb-2948-4ea0-adf2-8a8b108f7cf9

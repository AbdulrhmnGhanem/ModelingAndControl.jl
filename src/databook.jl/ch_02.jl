### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ d2096044-dc99-11ee-1393-9f98f23c7bc0
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
    data_path = joinpath("..", "..", "books", "databook", "DATA")
end;

# ╔═╡ 86d324dc-ffc9-4e17-9c07-f77a98bad08a
using FileIO, ImageShow, ImageCore, FFTW, Plots, Compat

# ╔═╡ 76eb5bfb-1760-4c92-a75c-f67dd7b6f85e
md"# Chapter 2 - Fourier and Wavelet Transforms"

# ╔═╡ deb0d340-ef17-403c-9f7d-5467c6a0d5a7
md"## Exercise 2.1
Load the image dog.jpg and convert to grayscale. Use the FFT to compress
the image at different compression ratios. Plot the error between the compressed and actual image as a function of the compression ratio.
"

# ╔═╡ 0258268e-bd13-4c5a-80d6-c93e8b5cb3d2
begin
	function solve_one(ratios)
		A = joinpath(data_path, "dog.jpg") |> load .|> Gray |> channelview |> float
		B = fft(A)
		B_abs = abs.(B)
		B_abs_sorted = sort(B_abs[:])
		
		compressed = Array{Matrix{Float32}}(undef, length(ratios))
		errors = similar(ratios)
		for (i, ratio) in enumerate(ratios)
			threshold = B_abs_sorted[floor((1 - ratio) * length(B_abs_sorted)) |> Int]
			mask = B_abs .> threshold
			coeffecients = B .* mask
			compressed[i] = real.(ifft(coeffecients))
			errors[i] = norm(A - compressed[i])
		end
		log10.(B_abs), compressed, errors, ratios
	end
	
	function plot_one(sol)
		coefficients, compressed, errors, ratios = sol
		p1 = plot(Gray.(coefficients);
			axis=([], false),
			title="coefficients",
		)
		ps = [
			plot(Gray.(compressed[i]);
				axis=([], false),
				title=round(ratio;digits=4),
			) for (i, ratio) in enumerate(ratios) ]

		p3 = plot(ratios, errors;
			xlabel="compression ratio",
			ylabel="error",
			legend=false,
		)
		
		plot(p1, ps..., p3;
			size=(900, 900),
		)
	end
end

# ╔═╡ 8b23a4c4-fb49-4e73-b833-538c8a992d54
plot_one(solve_one(logrange(0.001, 0.9, 10)))

# ╔═╡ 02bfe122-c0fa-461c-9531-168b4f09ec1f
md"## Exercise 2.3
Use the FFT to solve the Korteweg–de Vries (KdV) equation,

$u_t + u_{xxx} - uu_x = 0$
on a large domain with an initial condition u(x, 0) = sech(x). Plot the evolution.
"

# ╔═╡ 5232e813-eb09-4af1-bb0f-931c6febf0d2


# ╔═╡ 686fcd81-b13b-4aac-a1b0-aba51ed6a1ee
md"## Exercise 2.4
Use the FFT to solve the Kuramoto–Sivashinsky (KS) equation,

$u_t + u_{xx} + u_{xxxx} + \frac{1}{2} u^2_{x}$

on a large domain with an initial condition u(x, 0) = sech(x). Plot the evolution.
"

# ╔═╡ 81bafae5-e88d-41d7-9835-f92592a71e0e


# ╔═╡ bc581a79-8649-4aa1-a604-5d5bee07ac13
md"## Exercise 2.5
Solve for the analytic equilibrium temperature distribution using the 2D
Laplace equation on an L × H sized rectangular domain with the following boundary
conditions.

(a) Left: ux(0, y) = 0 (insulating).

(b) Bottom: u(x, 0) = 0 (fixed temperature).

(c) Top: u(x, H) = f (x) (zero temperature).

(d) Right: ux(L, y) = 0 (insulating).


![image](https://raw.githubusercontent.com/AbdulrhmnGhanem/ModelingAndControl.jl/main/src/databook.jl/temperature.png)

Solve for a general boundary temperature f (x). Also solve for a particular temperature
distribution f (x); you may choose any non-constant distribution you like.
How would this change if the left and right boundaries were fixed at zero temperature?
(You do not have to solve this new problem, just explain in words what would change.)
"

# ╔═╡ 58901914-fe97-4675-a11c-2f8c9d603c0f


# ╔═╡ db6e0b60-67b0-469d-b03c-1b4f966701af
md"## Exercise 2.6
Now, compute the solution to the 2D heat equation on a circular disk through
simulation. Recall that the heat equation is given by

$u_t = α^2 ∇^2 u$

For this problem, we will solve the heat equation using a finite-difference scheme on a Cartesian grid. We will use a grid of 300 × 300 with the circular disk in the center. The radius of the circle is r = 1, α = 1, and the domain is [−1.5, 1.5] in x and [−1.5, 1.5] in y. You can impose the boundary conditions by enforcing the temperature at points that are outside of the disk at the beginning of each new time step. It should be easy to find points that are outside the disk, because they satisfy $x^2 + y^2 = 1$.

Simulate the unsteady heat equation for the following boundary conditions:

1. The left half of the boundary of the disk is fixed at a temperature of u = 1 and the right half of the boundary is fixed at u = 2. Try simulating this with zero initial conditions first. Next, try initial conditions inside the disk where the top half is u = −1 and the bottom half is u = 1.

2. The temperature at the boundary of the disk is fixed at u(θ) = cos(θ).


Include your code and show some plots of your solutions to the heat equation. Plot the temperature distribution for each case (1) early on in the diffusion process, (2) near steady state, and (3) somewhere in the middle
"

# ╔═╡ 1cef26cc-0e99-454f-ae1b-a8d566f7d256


# ╔═╡ 9da009b7-6c99-4514-b9f7-94a5f29ea0ae
md"## Exercise 2.7
Consider the PDE for a vibrating string of finite length L,

$u_{tt} = c^2 u_{xx}, \ 0 \leq x \leq l,$

with the initial conditions,

$u(x, 0) = 0,\ u_t(x, 0) = 0,$

and boundary conditions

$u(0, t) = 0, u_x(L, t) = f (t).$

Solve this PDE by using the Laplace transform. You may keep your solution in the fre-
quency domain, since the inverse transform is complicated. Please simplify as much as
possible using functions like sinh and cosh.
This PDE cannot be solved by separation of variables. Why not? (That is, try to solve
with separation of variables until you hit a contradiction.)
"

# ╔═╡ fb0caf3f-2ae1-4a74-ad81-400fa4fe7f4f


# ╔═╡ f01710e1-489e-4913-8476-2fca7cc9a0a6
md"## Exercise 2.8

Now, we will use the FFT to simultaneously compress and re-
master an audio file. Please download the file r2112.mat and load the audio data into
the matrix rush and the sample rate FS.

1. Listen to the audio signal (»sound(rush,FS);). Compute the FFT of this audio signal.

2. Compute the power spectral density vector. Plot this to see what the output looks like. Also plot the spectrogram.

3. Now, download r2112noisy.mat and load this file to initialize the variable rushnoisy. This signal is corrupted with high-frequency artifacts. Manually zero the last three-quarters of the Fourier components of this noisy signal (if n=length(rushnoisy), then zero-out all Fourier coefficients from n/4:n). Use this filtered frequency spectrum to reconstruct the clean audio signal. When reconstructing, be sure to take the real part of the inverse Fourier transform: cleansignal=real(ifft(filteredcoefs));.

Because we are only keeping the first quarter of the frequency data, you must
multiply the reconstructed signal by 2 so that it has the correct normalized power.
Be sure to use the sound command to listen to the pre- and post-filtered versions.
Plot the power spectral density and spectrograms of the pre- and post-filtered signals.
"

# ╔═╡ 131297f3-8cac-4dab-8f8c-4c4ad29efd56


# ╔═╡ 16374e5a-d05a-450c-bec7-4a65b83141a3
md"## Exercise 2.9

The convolution integral and the impulse response may be used to simulate how an audio signal would sound under various conditions, such as in a long hallway, in a concert hall, or in a swimming pool.

The basic idea is that you can record the audio response to an impulsive sound in a given location, like a concert hall. For example, imagine that you put a microphone in the most expensive seats in the hall and then record the sound from a shotgun blast up on the stage. (Do not try this!!) Then, if you have a “flat” studio recording of some other audio, you can simulate how it would have sounded in the concert hall by convolving the two signals.

Download and unzip sounds.zip to find various sounds and impulse-response filters. Convolve the various audio files (labeled sound1.wav, . . . ) with the various filters (labeled FilterXYZ.wav, . . . ). In MATLAB, use the wavread command to load and the conv command to convolve. It is best to add 10% of the filtered audio (also known as “wet” audio) to 90% of the original audio (also known as “dry” audio). Listen to the filtered audio, as well as the original audio and the impulse-response filters (note that each sound has a sampling rate of FS=11,025). However, you will need to be careful when adding the 10% filtered and 90% unfiltered signals, since the filtered audio will not necessarily have the same length as the filtered audio.

There is a great video explaining how to actually create these impulse responses:
[www.audioease.com/Pages/Altiverb/sampling.php](www.audioease.com/Pages/Altiverb/sampling.php)
"

# ╔═╡ 4c7de71a-7bb3-4c66-a49c-d8634570bfab


# ╔═╡ Cell order:
# ╠═d2096044-dc99-11ee-1393-9f98f23c7bc0
# ╠═86d324dc-ffc9-4e17-9c07-f77a98bad08a
# ╟─76eb5bfb-1760-4c92-a75c-f67dd7b6f85e
# ╟─deb0d340-ef17-403c-9f7d-5467c6a0d5a7
# ╠═0258268e-bd13-4c5a-80d6-c93e8b5cb3d2
# ╟─8b23a4c4-fb49-4e73-b833-538c8a992d54
# ╟─02bfe122-c0fa-461c-9531-168b4f09ec1f
# ╠═5232e813-eb09-4af1-bb0f-931c6febf0d2
# ╟─686fcd81-b13b-4aac-a1b0-aba51ed6a1ee
# ╠═81bafae5-e88d-41d7-9835-f92592a71e0e
# ╟─bc581a79-8649-4aa1-a604-5d5bee07ac13
# ╠═58901914-fe97-4675-a11c-2f8c9d603c0f
# ╟─db6e0b60-67b0-469d-b03c-1b4f966701af
# ╠═1cef26cc-0e99-454f-ae1b-a8d566f7d256
# ╟─9da009b7-6c99-4514-b9f7-94a5f29ea0ae
# ╠═fb0caf3f-2ae1-4a74-ad81-400fa4fe7f4f
# ╟─f01710e1-489e-4913-8476-2fca7cc9a0a6
# ╠═131297f3-8cac-4dab-8f8c-4c4ad29efd56
# ╟─16374e5a-d05a-450c-bec7-4a65b83141a3
# ╠═4c7de71a-7bb3-4c66-a49c-d8634570bfab

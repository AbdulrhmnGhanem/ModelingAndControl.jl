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
    data_path = joinpath("..", "..", "books", "DATA")
end;

# ╔═╡ 86d324dc-ffc9-4e17-9c07-f77a98bad08a
# ╠═╡ show_logs = false
using FileIO, ImageShow, ImageCore, FFTW, Plots, Compat, DifferentialEquations, PlutoUI, MAT, Sound, DSP

# ╔═╡ 76eb5bfb-1760-4c92-a75c-f67dd7b6f85e
md"# Chapter 2 - Fourier and Wavelet Transforms"

# ╔═╡ 514d111e-163b-44f6-ae98-69b92113be06
md"## Solving heat equation using FFT
$u_t = \alpha^2u_{xx}$

where $u(t,\ x)$ is the temperature distribution in time and space.

$\mathcal{F}(u(t,\ x)) = \hat u(t, \ \kappa)$

$u_x \xrightarrow{\mathcal{F}} i k \hat u$


$\hat u_t = - \alpha^2 k^2 \hat u$
"

# ╔═╡ de62eabe-e4de-4694-8652-2a5d64772084
function heat()
    α = 1.0
    L = 100.0
    N = 1000
    dx = L / N
    domain = range(-L / 2, stop=L / 2 - dx, length=N)

    u₀ = collect(0 * domain)
    u₀[400:600] .= 1
    u₀ = fft(u₀, 1)

    κ = (2π / L) * (-N/2:N/2-1)
    κ = fftshift(κ, 1)

    tspan = (0, 1000)
    t = tspan[1]:0.1:tspan[2]
    params = (α, κ)

    function rhs(dû, û, p, t)
        α, k = p
        dû .= -α^2 * (k .^ 2) .* û
    end
    prob = ODEProblem(rhs, u₀, tspan, params)
    û = solve(prob)
    u = zeros(Complex, size(û)...)

    for k in 1:length(eachcol(u))
        u[:, k] = ifft(û[:, k])
    end
    sol = real(u)

    @gif for i in 1:length(eachcol(sol))
        plot(domain, sol[:, i];
            ylims=(-0.01, 1),
            legend=false,
            title="t = $(round(û.t[i] *1000; digits=2)) ms",
        )
    end
end

# ╔═╡ bbc60ff6-5199-40cc-bf82-7cb15fffdae7
# ╠═╡ show_logs = false
heat()

# ╔═╡ 6e3f102a-91a3-4540-bdce-aa615a7f3c2d
md"## Solving Burgers’ equation using FFT

$u_t + u u_x = νu_{xx}$
"

# ╔═╡ 9b024312-3adc-481e-a924-18364235e74d
function burgers()
    ν = 0.001
    L = 20
    N = 1000
    dx = L / N
    domain = range(-L / 2, stop=L / 2 - dx, length=N)

    κ = (2π / L) * (-N/2:N/2-1)
    κ = fftshift(κ, 1)
    u₀ = sech.(collect(domain))

    tspan = (0, 2.5)
    params = (κ, ν)

    function rhs(dudt, u, p, t)
        κ, ν = p
        û = fft(u)
        dû = im * κ .* û
        ddû = -(κ .^ 2) .* û
        du = ifft(dû)
        ddu = ifft(ddû)
        dudt .= real(-u .* du + ν * ddu)
    end

    prob = ODEProblem(rhs, u₀, tspan, params)
    u = solve(prob, reltol=1e-10, abstol=1e-10)
    sol, t = u[1:end, 1:end], u.t

    @gif for i in 1:length(eachcol(sol))
        plot(domain, sol[:, i];
            ylims=(-0.01, 1.1),
            legend=false,
            title="t = $(round(t[i] *1000; digits=2)) ms",
        )
    end
end

# ╔═╡ 3e7d55c7-bd24-442f-bbd1-b4d90509ae99
# ╠═╡ show_logs = false
burgers()

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
            threshold = B_abs_sorted[floor((1 - ratio) * length(B_abs_sorted))|>Int]
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
                title=round(ratio; digits=4),
            ) for (i, ratio) in enumerate(ratios)]

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
on a large domain with an initial condition $u(x, 0) = sech(x)$. Plot the evolution.
"

# ╔═╡ 5232e813-eb09-4af1-bb0f-931c6febf0d2
function solve_three()
    L = 100.0
    N = 1000
    dx = L / N
    domain = range(-L / 2, stop=L / 2 - dx, length=N)

    κ = (2π / L) * (-N/2:N/2-1)
    κ = fftshift(κ, 1)
    u₀ = sech.(collect(domain))

    tspan = (0, 4)

    function rhs(dudt, u, p, t)
        û = fft(u)
        dû = im * κ .* û
        dddû = -im * (κ .^ 3) .* û
        du = ifft(dû)
        dddu = ifft(dddû)
        dudt .= real(u .* du + dddu)
    end
    t = tspan[1]:tspan[2]/500:tspan[2]

    prob = ODEProblem(rhs, u₀, tspan)
    u = solve(prob, saveat=t, reltol=1e-10, abstol=1e-10)
    sol, t = u[1:end, 1:end], u.t


    @gif for i in 1:length(eachcol(sol))
        plot(domain, sol[:, i];
            ylims=(-0.25, 1),
            legend=false,
            title="t = $(round(t[i] *1000; digits=2)) ms",
        )
    end
end

# ╔═╡ 24cecda6-49ba-4958-b83a-e0cd579eebf3
# ╠═╡ show_logs = false
solve_three()

# ╔═╡ 686fcd81-b13b-4aac-a1b0-aba51ed6a1ee
md"## Exercise 2.4
Use the FFT to solve the Kuramoto–Sivashinsky (KS) equation,

$u_t + u_{xx} + u_{xxxx} + \frac{1}{2} u^2_{x}$

on a large domain with an initial condition $u(x, 0) = sech(x)$. Plot the evolution.
"

# ╔═╡ 81bafae5-e88d-41d7-9835-f92592a71e0e
function solve_four()
    L = 100.0
    N = 1000
    dx = L / N
    domain = range(-L / 2, stop=L / 2 - dx, length=N)

    κ = (2π / L) * (-N/2:N/2-1)
    κ = fftshift(κ, 1)
    u₀ = sech.(collect(domain))

    tspan = (0, 5)

    function rhs(dudt, u, p, t)
        û = fft(u)
        dû = im * κ .* û
        ddû = -(κ .^ 2) .* û
        ddddû = (κ .^ 4) .* û
        du = ifft(dû)
        ddu = ifft(ddû)
        ddddu = ifft(ddddû)
        dudt .= real(0.5 .* du .^ 2 - ddu - ddddu)
    end
    t = tspan[1]:tspan[2]/200:tspan[2]

    prob = ODEProblem(rhs, u₀, tspan)
    u = solve(prob, saveat=t)
    sol, t = u[1:end, 1:end], u.t

    @gif for i in 1:length(eachcol(sol))
        plot(domain, sol[:, i];
            ylims=(-0.5, 1.75),
            legend=false,
            title="t = $(round(t[i] *1000; digits=2)) ms",
        )
    end
end

# ╔═╡ 3a1c3f2e-7b7c-498b-9fac-beb6c89176f8
# ╠═╡ show_logs = false
solve_four()

# ╔═╡ bc581a79-8649-4aa1-a604-5d5bee07ac13
md"## Exercise 2.5
Solve for the analytic equilibrium temperature distribution using the 2D
Laplace equation on an L × H sized rectangular domain with the following boundary
conditions.

(a) Left: $u_x(0, y) = 0$ (insulating).

(b) Bottom: $u(x, 0) = 0$ (fixed temperature).

(c) Top: $u(x, H) = f(x)$ (zero temperature).

(d) Right: $u_x(L, y) = 0$ (insulating).


![image](https://raw.githubusercontent.com/AbdulrhmnGhanem/ModelingAndControl.jl/main/src/databook.jl/temperature.png)

Solve for a general boundary temperature $f(x)$. Also solve for a particular temperature distribution $f(x)$; you may choose any non-constant distribution you like. How would this change if the left and right boundaries were fixed at zero temperature? (You do not have to solve this new problem, just explain in words what would change.)
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
master an audio file. Please download the file _r2112.mat_ and load the audio data into the matrix _rush_ and the sample rate FS.

1. Listen to the audio signal `(»sound(rush,FS);)`. Compute the FFT of this audio signal.

2. Compute the power spectral density vector. Plot this to see what the output looks like. Also plot the spectrogram.

3. Now, download `r2112noisy.mat` and load this file to initialize the variable `rushnoisy`. This signal is corrupted with high-frequency artifacts. Manually zero the last three-quarters of the Fourier components of this noisy signal (if `n=length(rushnoisy)`, then zero-out all Fourier coefficients from `n/4:n`). Use this filtered frequency spectrum to reconstruct the clean audio signal. When reconstructing, be sure to take the real part of the inverse Fourier transform: `cleansignal=real(ifft(filteredcoefs));`.

Because we are only keeping the first quarter of the frequency data, you must
multiply the reconstructed signal by 2 so that it has the correct normalized power.
Be sure to use the sound command to listen to the pre- and post-filtered versions.
Plot the power spectral density and spectrograms of the pre- and post-filtered signals.
"

# ╔═╡ 131297f3-8cac-4dab-8f8c-4c4ad29efd56
function solve_eight()
    r2112 = matread(joinpath(data_path, "r2112.mat"))
    r2112_noisy = matread(joinpath(data_path, "r2112noisy.mat"))
    rush = r2112["rush"]
    rush_noisy = r2112_noisy["rushnoisy"]
    N = length(rush)
    fₛ = r2112["FS"]

    RUSH = fft(rush)
    psd = abs2.(RUSH) / (N * fₛ)
    RUSH_NOISY = fft(rush_noisy)
    psd_noisy = abs2.(RUSH_NOISY) / (N * fₛ)
    RUSH_CLEAN = deepcopy(RUSH_NOISY)
    RUSH_CLEAN[N÷4:N] .= 0
    rush_clean = 2 * ifft(real(RUSH_CLEAN))
    psd_clean = abs2.(fft(rush_clean)) / (N * fₛ)


    p1 = plot(1:fₛ, real(RUSH[1:N÷2]);
        legend=false,
        xlabel="Freq (Hz)",
        ylabel="Amplitude",
        title="original spectrum",
    )
    p2 = plot(1:fₛ, amp2db.(psd[1:N÷2]);
        legend=false,
        xlabel="Freq (Hz)",
        ylabel="PSD (dB / Hz)",
        title="original PSD",
    )
    p3 = plot(1:fₛ, real(RUSH_NOISY[1:N÷2]);
        legend=false,
        ylabel="Amplitude",
        xlabel="Freq (Hz)",
        title="noisy spectrum",
    )
    p4 = plot(1:fₛ, amp2db.(psd_noisy[1:N÷2]);
        legend=false,
        xlabel="Freq (Hz)",
        ylabel="PSD (dB / Hz)",
        title="noisy PSD",
    )

    p5 = plot(1:fₛ, real(RUSH_CLEAN[1:N÷2]);
        legend=false,
        ylabel="Amplitude",
        xlabel="Freq (Hz)",
        title="clean spectrum",
    )
    p6 = plot(1:fₛ, amp2db.(psd_clean[1:N÷2]);
        legend=false,
        xlabel="Freq (Hz)",
        ylabel="PSD (dB / Hz)",
        title="clean PSD",
        ylim=(-400, 0),
    )
    plot(p1, p2, p3, p4, p5, p6;
        layout=(3, 2),
        size=(600, 800)
    )
end

# ╔═╡ 2d13b62d-6ebe-41af-8ba7-c02d1bebbfc0
solve_eight()

# ╔═╡ 10eb7f7f-5a06-4ee9-bfba-fb5b9cc4a99c
md"
!!! remarks
	$\text{power spectral density (PSD)} = \frac{\text{power spectrum}}{\text{Equivalent Noise bandwidth (ENB)}}$
	$\text{PSD (for rectangular window)} = \frac{|\text{fft}|^2}{Nf_s}$
"

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
# ╟─514d111e-163b-44f6-ae98-69b92113be06
# ╠═de62eabe-e4de-4694-8652-2a5d64772084
# ╟─bbc60ff6-5199-40cc-bf82-7cb15fffdae7
# ╟─6e3f102a-91a3-4540-bdce-aa615a7f3c2d
# ╠═9b024312-3adc-481e-a924-18364235e74d
# ╟─3e7d55c7-bd24-442f-bbd1-b4d90509ae99
# ╟─deb0d340-ef17-403c-9f7d-5467c6a0d5a7
# ╠═0258268e-bd13-4c5a-80d6-c93e8b5cb3d2
# ╟─8b23a4c4-fb49-4e73-b833-538c8a992d54
# ╟─02bfe122-c0fa-461c-9531-168b4f09ec1f
# ╠═5232e813-eb09-4af1-bb0f-931c6febf0d2
# ╠═24cecda6-49ba-4958-b83a-e0cd579eebf3
# ╟─686fcd81-b13b-4aac-a1b0-aba51ed6a1ee
# ╠═81bafae5-e88d-41d7-9835-f92592a71e0e
# ╠═3a1c3f2e-7b7c-498b-9fac-beb6c89176f8
# ╟─bc581a79-8649-4aa1-a604-5d5bee07ac13
# ╠═58901914-fe97-4675-a11c-2f8c9d603c0f
# ╟─db6e0b60-67b0-469d-b03c-1b4f966701af
# ╠═1cef26cc-0e99-454f-ae1b-a8d566f7d256
# ╟─9da009b7-6c99-4514-b9f7-94a5f29ea0ae
# ╠═fb0caf3f-2ae1-4a74-ad81-400fa4fe7f4f
# ╟─f01710e1-489e-4913-8476-2fca7cc9a0a6
# ╠═131297f3-8cac-4dab-8f8c-4c4ad29efd56
# ╟─2d13b62d-6ebe-41af-8ba7-c02d1bebbfc0
# ╟─10eb7f7f-5a06-4ee9-bfba-fb5b9cc4a99c
# ╟─16374e5a-d05a-450c-bec7-4a65b83141a3
# ╠═4c7de71a-7bb3-4c66-a49c-d8634570bfab

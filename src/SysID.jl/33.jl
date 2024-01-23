### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ c87ae1ea-a384-11ee-1623-f92e422d4ccb
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 5ff84bfc-d93c-4d81-b7ef-a9ff5b4d649f
using Plots, DSP, StatsBase, FFTW

# ╔═╡ 42809776-c132-4a13-b4e0-2e697adee2df
md"# Exercise 33: Exploting the periodic nature of signals (differentiation, integration, averaging, and filtering)

!!! purpose
	Using frequncy domain approach to find the derivative of a sawtooth function.

	- After adding white noise to the signal, the derivative is much noiser than the original siganl because the differentiation action amplified teh high frequency noise.

	- Averaging over M periods reduces the noise by $\sqrt M$. Eliminating all even frequencies reduces the noise one more with a factor of $\sqrt 2$.
"

# ╔═╡ fc431c55-1b9c-47ce-9739-2b3c244ac993
function butter(n, wn)
    digitalfilter(Lowpass(wn), Butterworth(n))
end

# ╔═╡ 9969c358-3964-457e-bdf5-77bcea4178b8
begin
    N = 2^10
    fₛ = N
    M = 10
    h = butter(2, 0.3)
    SNR = 0.001

    u₀ = [0:N/2-1; N/2:-1:1]
    u₀ .-= mean(u₀)
    u₀ /= std(u₀)
    u₀ = repeat(u₀, M + 1)
    u₀ = filt(h, u₀)
    deleteat!(u₀, 1:N)
    u = u₀ + SNR * randn(length(u₀))

    n = length(u)
    f = 0:n-1
    ω = 2π * f

    # derivative of original signal
    U₀ = fft(u₀)
    U₀[n÷2:end] .= 0
    dU₀ = U₀ .* ω .* im
    du₀ = 2ifft(dU₀) .|> real

    # derivative of the noisy signal
    U = fft(u)
    U[n÷2:end] .= 0
    dU = U .* ω .* im
    du = 2ifft(dU) .|> real

    du_averaged = mean(reshape(du, N, M); dims = 2)  # average over M periods

    dU_averaged = fft(du_averaged)
    dU_averaged[1:2:end] .= 0  # filter even harmonics
    du_averaged_filtered = ifft(dU_averaged) .|> real
end;

# ╔═╡ 49a9eac6-8657-4fb9-b907-3afa79097d4f
begin
    xlabel = "Time (ms)"
    plot(u₀[1:N]; title = "Noiseless", label = "Original", xlabel)
    plot!(du₀[1:N]; label = "Derivative")
end

# ╔═╡ 623f3382-29d0-4077-8631-910530059c77
begin
    plot(u[1:N]; title = "Noisy", label = "Origianl", xlabel)
    plot!(du[1:N]; label = "Derivative")
end

# ╔═╡ bb7c83a9-c80b-41e5-a0bf-59055007ea34
begin
    plot(u[1:N]; title = "Noisy (averaged)", label = "Original", xlabel)
    plot!(du_averaged[1:N], label = "Derivative")
end

# ╔═╡ 92ed3e59-08f0-44b9-aaef-7efbaab46d23
begin
    plot(
        u[1:N];
        title = "Noisy (averaged and only odd harmonics)",
        label = "Original",
        xlabel,
    )
    plot!(du_averaged_filtered[1:N], label = "Derivative")
end

# ╔═╡ Cell order:
# ╠═c87ae1ea-a384-11ee-1623-f92e422d4ccb
# ╠═5ff84bfc-d93c-4d81-b7ef-a9ff5b4d649f
# ╟─42809776-c132-4a13-b4e0-2e697adee2df
# ╠═fc431c55-1b9c-47ce-9739-2b3c244ac993
# ╠═9969c358-3964-457e-bdf5-77bcea4178b8
# ╟─49a9eac6-8657-4fb9-b907-3afa79097d4f
# ╟─623f3382-29d0-4077-8631-910530059c77
# ╟─bb7c83a9-c80b-41e5-a0bf-59055007ea34
# ╟─92ed3e59-08f0-44b9-aaef-7efbaab46d23

### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 62e757e9-5279-42c7-86ce-0dc0a96ff426
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ cd3fdf9d-7285-4b24-96db-22edcf015835
using Plots, DSP, FFTW, StatsBase, PlutoUI

# ╔═╡ cd6ea107-5923-4239-a2fe-6da1a6dffcd8
md"# Exercise 39: FRF measurement using a noise excitation and a Hanning window

!!! purpose
	- Using windowing makes the error significantly ($G_h - G_0$) samlller at most frequencies.
	- For low number of averages ($M = 1$), the stochastic leakage errors dominate.
	- By increasing the number of averages, the stochastic errors are averaged to zero (decreases the variance, and visually it becomes less scattered) and the remaining errors are dominated by the bias errors of the leakage effect.
	- The systematic errors are dominant around the resonance frequency.
	- The level of the systematic error has $O(N^{-2})$ and the variance has $O(M^{-1} N^{-2})$ where $N$ is window size, and $M$ is the number of averages.
"

# ╔═╡ d155c6e6-7a88-4b6b-9762-a148a1e6415c
function butter(n, wn)
    h = digitalfilter(Lowpass(wn), Butterworth(n))
    convert(PolynomialRatio, h)
end

# ╔═╡ 62c5feb4-250b-4341-9721-fd65f019f519
function cheby1(n, r, wp)
    h = digitalfilter(Lowpass(wp), Chebyshev1(n, r))
    convert(PolynomialRatio, h)
end

# ╔═╡ 94c099c2-a626-11ee-3cc4-4504b89da241
begin
    Ms = [1, 4, 16, 256, 1024, 4096]
    M = maximum(Ms)
    N = 1024
    fₛ = 128
    Nₜᵣₐₙₛ = 1024

    t = LinRange(0, (N - 1) / fₛ, N)
    Lines = 1:N÷2
    f = (Lines .- 1) / N * fₛ
    h = cheby1(2, 10, 0.2)
    h_gen = butter(2, 0.3 * 2)

    # random noise
    u = randn(N * M + Nₜᵣₐₙₛ)
    u = filt(h_gen, u)
    y = filt(h, u)
    # remove the transient response
    deleteat!(u, 1:Nₜᵣₐₙₛ)
    deleteat!(y, 1:Nₜᵣₐₙₛ)

    windows = repeat(hanning(N), 1, M)
    # repeated excitations
    y = reshape(y, N, M)
    Y = fft(y, 1) / sqrt(N / 2)
    Y = Y[Lines, :]

    Yₕ = fft(y .* windows, 1) / sqrt(N / 2)
    Yₕ = Yₕ[Lines, :]

    u = reshape(u, N, M)
    U = fft(u, 1) / sqrt(N / 2)
    U = U[Lines, :]

    Uₕ = fft(u .* windows, 1) / sqrt(N / 2)
    Uₕ = Uₕ[Lines, :]

    G₀, w = freqresp(h)

    G = zeros(ComplexF64, length(Lines), length(Ms))
    Gₕ = similar(G)
    UU = similar(G)
    UUₕ = similar(G)
    YU = similar(G)
    YUₕ = similar(G)

    for m = 1:length(Ms)
        UU[:, m] .= mean(abs2.(U[:, 1:Ms[m]]), dims = 2)
        YU[:, m] .= mean(Y[:, 1:Ms[m]] .* conj.(U[:, 1:Ms[m]]), dims = 2)
        G[:, m] .= YU[:, m] ./ UU[:, m]

        UUₕ[:, m] .= mean(abs2.(Uₕ[:, 1:Ms[m]]), dims = 2)
        YUₕ[:, m] .= mean(Yₕ[:, 1:Ms[m]] .* conj.(Uₕ[:, 1:Ms[m]]), dims = 2)
        Gₕ[:, m] .= YUₕ[:, m] ./ UUₕ[:, m]
    end
end

# ╔═╡ e2c01037-b881-458f-9fba-9129a7a90ddf
begin
    mapped_f = f[1:2:end-1]

    ps = []
    for m = 1:length(Ms)
        p = plot(
            mapped_f,
            G₀ .|> abs .|> amp2db;
            ylims = (-80, 0),
            legend = false,
            title = "M = $(Ms[m])",
            xlabel = m > 3 ? "Frequency (Hz)" : "",
            ylabel = m == 1 || m == 4 ? "Amplitude (dB)" : "",
            label = "M = $(Ms[m])",
        )
        scatter!(mapped_f, (Gₕ[1:2:end, m] - G₀[1:1:end-1]) .|> abs .|> amp2db;)
        p = scatter!(mapped_f, (G[1:2:end, m] - G₀[1:1:end-1]) .|> abs .|> amp2db)
        push!(ps, p)
    end
    legend = scatter(
        (1:3)';
        ylims = (0, 0.00001),
        framestyle = :none,
        label = reshape(["G₀", "Gₕ - G₀", "G - G₀"], 1, :),
    )
    push!(ps, legend)

    plot(ps...; layout = (3, 3), size = (700, 500))
end

# ╔═╡ 476b4eb8-b910-4498-81bc-3ffe3a342608
md" ### Here we study the effect of the window size when the number averages is constant."

# ╔═╡ 46714228-5063-422d-b568-c68f97c11338
function part2()
    # Just scoping the second part of the exercise to avoid overwrting variables
    M = 128
    Ns = [128, 256, 512, 1024, 2048, 4096]
    fₛ = 128
    Nₜᵣₐₙₛ = 1024

    for (i, N) in enumerate(Ns)
        t = LinRange(0, (N - 1) / fₛ, N)
        Lines = 1:N÷2
        f = (Lines .- 1) / N * fₛ
        h = cheby1(2, 10, 0.2)
        h_gen = butter(2, 0.3 * 2)

        # random noise
        u = randn(N * M + Nₜᵣₐₙₛ)
        u = filt(h_gen, u)
        y = filt(h, u)
        # remove the transient response
        deleteat!(u, 1:Nₜᵣₐₙₛ)
        deleteat!(y, 1:Nₜᵣₐₙₛ)

        windows = repeat(hanning(N), 1, M)
        # repeated excitations
        y = reshape(y, N, M)
        Y = fft(y, 1) / sqrt(N / 2)
        Y = Y[Lines, :]

        Yₕ = fft(y .* windows, 1) / sqrt(N / 2)
        Yₕ = Yₕ[Lines, :]

        u = reshape(u, N, M)
        U = fft(u, 1) / sqrt(N / 2)
        U = U[Lines, :]

        Uₕ = fft(u .* windows, 1) / sqrt(N / 2)
        Uₕ = Uₕ[Lines, :]

        G₀, w = freqresp(h)

        G = zeros(ComplexF64, length(Lines), length(Ns))
        Gₕ = similar(G)
        UU = similar(G)
        UUₕ = similar(G)
        YU = similar(G)
        YUₕ = similar(G)

        UU[:, i] .= mean(abs2.(U[:, 1:M]), dims = 2)
        UU[:, i] .= mean(abs2.(U[:, 1:M]), dims = 2)
        YU[:, i] .= mean(Y[:, 1:M] .* conj.(U[:, 1:M]), dims = 2)
        G[:, i] .= YU[:, i] ./ UU[:, i]

        UUₕ[:, i] .= mean(abs2.(Uₕ[:, 1:M]), dims = 2)
        YUₕ[:, i] .= mean(Yₕ[:, 1:M] .* conj.(Uₕ[:, 1:M]), dims = 2)
        Gₕ[:, i] .= YUₕ[:, i] ./ UUₕ[:, i]
    end
    return G₀, G, Gₕ, Ns
end

# ╔═╡ 58b89061-f26a-48c8-8137-b6953d7cea07
begin
    G₀2, G2, Gₕ2, Ns = part2()
    ps2 = []
    for n = 1:length(Ns)
        p = plot(
            mapped_f,
            G₀2 .|> abs .|> amp2db;
            ylims = (-80, 0),
            legend = false,
            title = "N = $(Ns[n])",
            xlabel = n > 3 ? "Frequency (Hz)" : "",
            ylabel = n == 1 || n == 4 ? "Amplitude (dB)" : "",
            label = "M = $(Ns[n])",
        )
        scatter!(mapped_f, (Gₕ2[1:2:end, n] - G₀2[1:1:end-1]) .|> abs .|> amp2db;)
        p = scatter!(mapped_f, (G2[1:2:end, n] - G₀2[1:1:end-1]) .|> abs .|> amp2db)
        push!(ps2, p)
    end
    push!(ps2, legend)

    plot(ps2...; layout = (3, 3), size = (700, 500))
end

# ╔═╡ Cell order:
# ╠═62e757e9-5279-42c7-86ce-0dc0a96ff426
# ╠═cd3fdf9d-7285-4b24-96db-22edcf015835
# ╟─cd6ea107-5923-4239-a2fe-6da1a6dffcd8
# ╠═d155c6e6-7a88-4b6b-9762-a148a1e6415c
# ╠═62c5feb4-250b-4341-9721-fd65f019f519
# ╠═94c099c2-a626-11ee-3cc4-4504b89da241
# ╠═e2c01037-b881-458f-9fba-9129a7a90ddf
# ╟─476b4eb8-b910-4498-81bc-3ffe3a342608
# ╠═46714228-5063-422d-b568-c68f97c11338
# ╟─58b89061-f26a-48c8-8137-b6953d7cea07

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
using Plots, DSP, FFTW, StatsBase

# ╔═╡ cd6ea107-5923-4239-a2fe-6da1a6dffcd8
md"# Exercise 39: FRF measurement using a noise excitation and a diff window
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

    # repeated excitations
    y = reshape(y, N, M)
    Y = fft(y, 1) / sqrt(N / 2)
    Y = Y[Lines, :]
    Y = diff(Y, dims = 1)

    u = reshape(u, N, M)
    U = fft(u, 1) / sqrt(N / 2)
    U = U[Lines, :]
    U = diff(U, dims = 1)

    G₀, w = freqresp(h)

    G = zeros(ComplexF64, length(Lines) - 1, length(Ms))
    UU = similar(G)
    YU = similar(G)

    for m = 1:length(Ms)
        # println(size(mean(abs2.(U[:, 1:Ms[m]]), dims=2)))

        UU[:, m] .= mean(abs2.(U[:, 1:Ms[m]]), dims = 2)
        YU[:, m] .= mean(Y[:, 1:Ms[m]] .* conj.(U[:, 1:Ms[m]]), dims = 2)
        G[:, m] .= YU[:, m] ./ UU[:, m]
    end
end

# ╔═╡ e2c01037-b881-458f-9fba-9129a7a90ddf
begin
    mapped_f = f[1:2:end-1]

    ps = []
    for m = 1:length(Ms)
        plot(
            mapped_f,
            G₀ .|> abs .|> amp2db;
            ylims = (-60, 0),
            legend = false,
            title = "M = $(Ms[m])",
            xlabel = m > 3 ? "Frequency (Hz)" : "",
            ylabel = m == 1 || m == 4 ? "Amplitude (dB)" : "",
            label = "M = $(Ms[m])",
        )
        p = scatter!(mapped_f, (G[1:2:end, m] - G₀[1:1:end-1]) .|> abs .|> amp2db)
        push!(ps, p)
    end
    legend = scatter(
        (1:2)';
        ylims = (0, 0.00001),
        framestyle = :none,
        label = reshape(["G₀", "G - G₀"], 1, :),
    )
    push!(ps, legend)

    plot(ps...; layout = (3, 3), size = (700, 500))
end

# ╔═╡ Cell order:
# ╠═62e757e9-5279-42c7-86ce-0dc0a96ff426
# ╠═cd3fdf9d-7285-4b24-96db-22edcf015835
# ╟─cd6ea107-5923-4239-a2fe-6da1a6dffcd8
# ╠═d155c6e6-7a88-4b6b-9762-a148a1e6415c
# ╠═62c5feb4-250b-4341-9721-fd65f019f519
# ╠═94c099c2-a626-11ee-3cc4-4504b89da241
# ╠═e2c01037-b881-458f-9fba-9129a7a90ddf

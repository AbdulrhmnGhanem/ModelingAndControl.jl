### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ dc2d1942-accd-11ee-34f8-ef2f689c65b2
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 6ce6c76b-c24b-4a6b-8297-f62cdce1191e
using Plots, FFTW, DSP, Distributions

# ╔═╡ 40338bfb-3c0d-4990-b84b-0ca6d42373c0
md"## Exercise 43: Measurement of the FRF using a random noise sequence and a random phase multisine in the presence of output noise.

!!! purpose
	Averaging over multiple realizations reduces the leakage error.
"

# ╔═╡ 0345624c-b9f8-4a13-9371-f08f428523ce
function cheby1(n, r, wp)
    digitalfilter(Lowpass(wp), Chebyshev1(n, r))
end

# ╔═╡ 7c4713b0-e397-4f06-9abf-77e73480a671
function random_phase_mutli_sine(excited_harm, N, total_size)
    u = Vector{Float64}(undef, total_size)  # Pre-allocate space for u
    S = zeros(ComplexF64, N)

    for i = 1:total_size÷N
        S[excited_harm.+1] .= exp.(im .* 2 * π * rand(size(excited_harm)))
        r = 2 * real(ifft(S))
        r /= std(r)  # rms value = 1
        u[(i-1)*N+1:i*N] .= r
    end
    u
end

# ╔═╡ 3d1a9815-e856-4c2d-a7ab-2ff2d5daadd2
begin
    period = 1024
    fₛ = 128
    t = LinRange(0, (period - 1) / fₛ, period)
    Lines = 1:period÷2
    f = (Lines .- 1) / period * fₛ
    Ms = [1, 4, 16, 64]
    A = 0.33
    σ = 0.02
    Nₜᵣₐₙₛ = period

    filter_order = 2
    filter_resonance = 10
    filter_cutoff = 0.1fₛ

    h = cheby1(filter_order, filter_resonance, 2filter_cutoff / fₛ)
    H, w = freqresp(h)


    u = []
    y = []
    yₙ = []
    GDiff = zeros(ComplexF64, length(Lines) - 1, length(Ms))
    G₀Diff = similar(GDiff)
    GSin = zeros(ComplexF64, length(Lines), length(Ms))

    for (i, m) in enumerate(Ms)
        N = m * period
        # Excite the system with random input
        u = A * randn(N + Nₜᵣₐₙₛ)
        y₀ = filt(h, u)
        # Add noise to the output
        yₙ = σ * randn(size(y₀))
        y = y₀ + yₙ
        # Remove transient
        deleteat!(u, 1:Nₜᵣₐₙₛ)
        deleteat!(y₀, 1:Nₜᵣₐₙₛ)
        deleteat!(y, 1:Nₜᵣₐₙₛ)
        # Average over `m` realization
        u = reshape(u, period, m)
        y = reshape(y, period, m)
        y₀ = reshape(y₀, period, m)

        # FRF
        U = fft(u, 1) / (period / 2)
        U = U[Lines, :]
        U = diff(U, dims = 1)

        Y = fft(y, 1) / (period / 2)
        Y = Y[Lines, :]
        Y = diff(Y, dims = 1)

        Y₀ = fft(y₀, 1) / (period / 2)
        Y₀ = Y₀[Lines, :]
        Y₀ = diff(Y₀, dims = 1)

        UU = mean(abs2.(U[:, 1:m]), dims = 2)
        YU = mean(Y[:, 1:m] .* conj.(U[:, 1:m]), dims = 2)
        Y₀U = mean(Y₀[:, 1:m] .* conj.(U[:, 1:m]), dims = 2)
        GDiff[:, i] .= YU ./ UU
        G₀Diff[:, i] .= Y₀U ./ UU

        # Excite the system with multisine
        u_sin = A * random_phase_mutli_sine(2:period÷2-1, period, N + Nₜᵣₐₙₛ)
        y₀_sin = filt(h, u_sin)
        # Add noise to the output
        yₙ = σ * randn(size(y₀_sin))
        y_sin = y₀_sin + yₙ
        # Remove transient
        deleteat!(u_sin, 1:Nₜᵣₐₙₛ)
        deleteat!(y₀_sin, 1:Nₜᵣₐₙₛ)
        deleteat!(y_sin, 1:Nₜᵣₐₙₛ)
        # Average over `m` realization
        u_sin = reshape(u_sin, period, m)
        y_sin = reshape(y_sin, period, m)
        y₀_sin = reshape(y₀_sin, period, m)

        Y_sin = fft(y_sin) / (period / 2)
        Y_sin = Y_sin[Lines, :]
        U_sin = fft(u_sin) / (period / 2)
        U_sin = U_sin[Lines, :]
        GSin[:, i] = mean(Y_sin, dims = 2) ./ mean(U_sin, dims = 2)
    end
end;

# ╔═╡ 465cf76b-5322-4cbc-8da4-756a8c167edf
md"### Time plot"

# ╔═╡ e98ec1e2-f689-4dce-aaa5-6d74fe5c2540
begin
    p = plot(
        t,
        u[:, 1];
        xlabel = "Time (s)",
        ylabel = "Amplitude",
        ylims = (-1, 1),
        labels = reshape(["u(t)"], 1, :),
    )
    p2 = plot(
        t,
        [y[:, 1], yₙ[1:period]];
        xlabel = "Time (s)",
        ylims = (-1, 1),
        labels = reshape(["y(t)", "r(t)"], 1, :),
    )
    plot(p, p2; size = (600, 250))
end

# ╔═╡ 0c72b0c6-4e86-4222-a125-b4b7cb0471ec
md"### FRF"

# ╔═╡ 016c4ff5-5313-4bcd-94ba-a513aaa1b465
begin
    mapped_f = f[1:end-1]
    mapped_f_sin = mapped_f[1:2:end]
    G₀ = H[1:end-1]

    ps = []
    for m = 1:length(Ms)
        local p = scatter(
            mapped_f,
            (GDiff[:, m] - G₀Diff[:, m]) .|> abs .|> amp2db;
            ylims = (-60, 0),
            legend = false,
            title = "M = $(Ms[m])",
            xlabel = m > 2 ? "Frequency (Hz)" : "",
            ylabel = m == 1 || m == 3 ? "Amplitude (dB)" : "",
            label = "M = $(Ms[m])",
        )
        scatter!(mapped_f_sin, (GSin[:, m][1:end-1] - G₀Diff[:, 4]) .|> abs .|> amp2db;)
        plot!(mapped_f_sin, G₀ .|> abs .|> amp2db)
        push!(ps, p)
    end
    legend = scatter(
        (1:3)';
        ylims = (0, 0.00001),
        framestyle = :none,
        label = reshape(["G rand with diff win - G₀", "G multisinse - G₀", "G₀"], 1, :),
    )
    push!(ps, legend)

    plot(ps...; size = (700, 500), layout = (3, 2))
end

# ╔═╡ Cell order:
# ╠═dc2d1942-accd-11ee-34f8-ef2f689c65b2
# ╠═6ce6c76b-c24b-4a6b-8297-f62cdce1191e
# ╟─40338bfb-3c0d-4990-b84b-0ca6d42373c0
# ╠═0345624c-b9f8-4a13-9371-f08f428523ce
# ╠═7c4713b0-e397-4f06-9abf-77e73480a671
# ╠═3d1a9815-e856-4c2d-a7ab-2ff2d5daadd2
# ╟─465cf76b-5322-4cbc-8da4-756a8c167edf
# ╟─e98ec1e2-f689-4dce-aaa5-6d74fe5c2540
# ╟─0c72b0c6-4e86-4222-a125-b4b7cb0471ec
# ╠═016c4ff5-5313-4bcd-94ba-a513aaa1b465

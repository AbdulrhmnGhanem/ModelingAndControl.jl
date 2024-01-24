### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

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
md"# Exercise 37: FRF measurement using a noise excitation and rectangular window

!!! purpose
	- The errors in the FRF are completely due to leakage.
	- The initial error drops fast with growing `M` (number of average realizations).
	- For larger `M` values the error becomse proportional to $\frac{1}{\sqrt M}$.
	- Normalizing the FFT output slightly improves the FRF when `M` is small.
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

# ╔═╡ 665c8045-7d0b-4dac-8485-fdc458882510
md"### Normalize FFT"

# ╔═╡ b2165dd4-d9be-4add-9039-c796b58ec95d
@bind normalize CheckBox(default = true)

# ╔═╡ 94c099c2-a626-11ee-3cc4-4504b89da241
begin
    Ms = [1, 4, 16, 256, 1024, 4096]
    M = maximum(Ms)
    N = 128
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
    Y = fft(y, 1) / (normalize ? sqrt(N / 2) : 1)
    Y = Y[Lines, :]
    u = reshape(u, N, M)
    U = fft(u, 1) / (normalize ? sqrt(N / 2) : 1)
    U = U[Lines, :]

    G₀, w = freqresp(h)

    G = zeros(ComplexF64, length(Lines), length(Ms))
    UU = similar(G)
    YU = similar(G)
    for m = 1:length(Ms)
        UU[:, m] .= mean(abs2.(U[:, 1:Ms[m]]), dims = 2)
        YU[:, m] .= mean(Y[:, 1:Ms[m]] .* conj.(U[:, 1:Ms[m]]), dims = 2)
        G[:, m] .= YU[:, m] ./ UU[:, m]
    end
end

# ╔═╡ c6cb6957-2c2b-4939-bdfe-faa3fa049850
md"## FRF"

# ╔═╡ e2c01037-b881-458f-9fba-9129a7a90ddf
begin
    G₀_mapped_to_f = G₀[1:4:end-1]

    ps = []
    for m = 1:length(Ms)
        plot(f, G₀_mapped_to_f .|> abs .|> amp2db; ylims = (-60, 0), legend = false)
        scatter!(f, G[:, m] .|> abs .|> amp2db, label = "M = $(Ms[m])")
        p = scatter!(
            f,
            (G[:, m] - G₀_mapped_to_f) .|> abs .|> amp2db;
            title = "M = $(Ms[m])",
            xlabel = m > 3 ? "Frequency (Hz)" : "",
            ylabel = m == 1 || m == 4 ? "Amplitude (dB)" : "",
        )
        push!(ps, p)
    end
    legend = scatter(
        (1:3)';
        ylims = (0, 0.00001),
        framestyle = :none,
        label = reshape(["G₀", "G", "G - G₀"], 1, :),
    )
    push!(ps, legend)

    plot(ps...; layout = (3, 3), size = (1000, 500))
end

# ╔═╡ da756326-8389-4747-b4de-4747e0f6b4d1
md"## Input power spectrum averaged over `M` realizations"

# ╔═╡ cfcea01e-c5f3-4f96-9beb-b3f88156c938
begin
    ps2 = []
    for m = 1:length(Ms[1:4])
        # why are we dividing by 2?
        p = plot(
            f,
            UU[:, m] / 2 .|> abs .|> amp2db;
            title = "M = $(Ms[m])",
            xlabel = m > 2 ? "Frequency (Hz)" : "",
            ylabel = m == 1 || m == 3 ? "Amplitude (dB)" : "",
        )
        push!(ps2, p)
    end
    plot(ps2...; legend = false)
end

# ╔═╡ Cell order:
# ╠═62e757e9-5279-42c7-86ce-0dc0a96ff426
# ╠═cd3fdf9d-7285-4b24-96db-22edcf015835
# ╟─cd6ea107-5923-4239-a2fe-6da1a6dffcd8
# ╠═d155c6e6-7a88-4b6b-9762-a148a1e6415c
# ╠═62c5feb4-250b-4341-9721-fd65f019f519
# ╠═94c099c2-a626-11ee-3cc4-4504b89da241
# ╟─665c8045-7d0b-4dac-8485-fdc458882510
# ╟─b2165dd4-d9be-4add-9039-c796b58ec95d
# ╟─c6cb6957-2c2b-4939-bdfe-faa3fa049850
# ╟─e2c01037-b881-458f-9fba-9129a7a90ddf
# ╟─da756326-8389-4747-b4de-4747e0f6b4d1
# ╟─cfcea01e-c5f3-4f96-9beb-b3f88156c938

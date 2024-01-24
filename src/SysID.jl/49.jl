### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 42c6fcac-b935-11ee-1e01-cde2ae4b39e8
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ b12c3689-fbed-4803-b112-eb19c7113696
using Plots, DSP, FFTW, StatsBase

# ╔═╡ 6597d608-4a5d-4ef8-9fff-57b690333c0f
begin
    function cheby1(n, r, wp)
        h = digitalfilter(Lowpass(wp), Chebyshev1(n, r))
        tf = convert(PolynomialRatio, h)
        coefb(tf), coefa(tf)
    end

    function freqz(
        filt::Union{FilterCoefficients,Vector},
        fₛ::Int,
        lower::Int,
        upper::Int,
        NPeriod::Int,
    )
        frequencies = range(lower, upper, length = NPeriod)
        if filt isa Vector
            @assert length(filt) == 2
            filt = PolynomialRatio(filt[1], filt[2])
        end
        H = freqresp(filt, frequencies * 2π / fₛ)
        return frequencies, H
    end

    function random_phase_mutli_sine(excited_harm, N, total_size)
        u = Vector{Float64}(undef, total_size)
        S = zeros(ComplexF64, N)

        for i = 1:total_size÷N
            S[excited_harm.+1] .= exp.(im .* 2 * π * rand(size(excited_harm)))
            r = 2 * real(ifft(S))
            r /= std(r)  # rms value = 1
            u[(i-1)*N+1:i*N] .= r
        end
        u
    end
    md"### Matlab compatability utils"
end

# ╔═╡ 4c3e735b-2cd4-47dd-a3fc-a823e5ab6176
md"## Exercise 49: Direct measurement of FRF under feedback conditions
!!! purpose
	- The problem with identifying a with system with feedback configuration is that the process noise is feedback to the input of the system. This creates nonzero cross correlation between the input and the output.

	- In many cases the exact feedback configuration is not known or the user might be unaware of that the system is operating under feedback conditions.

	- If feedback isn't considered it results in systematic error.

	- In this exercise we use direct method which results in an approxiamated system with an expected value of
	
	$\mathbb{E}\{\hat G_{direct}(k) \} = \frac{G_{FF}(k)S_{rr}(k) - G^*_{FB}(k)S_{vv}(k)}{S_{rr}(k) + |G_{FB}|^2 S_{vv}(k)}$ 

	`Sᵥᵥ`, `Sᵣᵣ` are the power spectrum of the process noise and reference signal respectively.

	- The feedforwad $G_{FF}$ is identified when the `r` dominates over `v`
	- The inverse of feedback $G^{-1}_{FB}$ is identified when `v` dominates.
	- Usually a mixture of both is obtained resulting in a baised view of $G_{FF}$.
"

# ╔═╡ 05820097-a0b2-4196-9d85-16501b9d9777
begin
    bFF, aFF = cheby1(2, 20, 0.5)
    bFF[1] = 0
    bFF *= 2
    br = copy(bFF)
    ar = aFF + bFF
    bv = copy(aFF)
    av = aFF + bFF


    fₛ = 128
    NPeriod = 1024
    M = 256
    N = NPeriod * M
    wo, GD0 = freqz([bFF, aFF], fₛ, 0, 55, NPeriod)  # exact FRF

    As = [0.0, 0.2, 0.4, 0.8]
    σ = 0.2
    NTrans = 1024
    Lines = 1:NPeriod÷2
    f = (Lines .- 1) / NPeriod * fₛ
    fD = f .+ 0.5 / NPeriod * fₛ

    GD = zeros(ComplexF64, NPeriod ÷ 2, length(As))
    for (s, A) in enumerate(As)
        r = randn(N + NTrans) * A
        v = randn(N + NTrans) * σ

        y = filt(br, ar, r) + filt(bv, av, v)
        u = r - y

        deleteat!(y, 1:NTrans)
        deleteat!(u, 1:NTrans)

        u = reshape(u, NPeriod, M)
        y = reshape(y, NPeriod, M)

        Y = fft(y, 1) / (NPeriod / 2)
        U = fft(u, 1) / (NPeriod / 2)
        Y = diff(Y; dims = 1)
        U = diff(U; dims = 1)
        YD = Y[Lines, :]
        UD = U[Lines, :]

        UU = mean(abs2.(UD[:, 1:M]), dims = 2)
        YU = mean(YD[:, 1:M] .* conj.(UD[:, 1:M]), dims = 2)
        GD[:, s] = YU ./ UU
    end

    GMulti = zeros(ComplexF64, NPeriod ÷ 2 - 2, length(As))

    Lines = 2:NPeriod÷2-1
    f = (Lines .- 1) / NPeriod * fₛ

    for (s, A) in enumerate(As)
        r = random_phase_mutli_sine(Lines, NPeriod, N + NTrans) * A
        v = randn(N + NTrans) * σ

        y = filt(br, ar, r) + filt(bv, av, v)
        u = r - y

        deleteat!(y, 1:NTrans)
        deleteat!(u, 1:NTrans)

        u = reshape(u, NPeriod, M)
        y = reshape(y, NPeriod, M)

        Y = fft(y, 1) / (NPeriod / 2)
        U = fft(u, 1) / (NPeriod / 2)
        Y = Y[Lines, :]
        U = U[Lines, :]

        Um = mean(U, dims = 2)
        Ym = mean(Y, dims = 2)
        GMulti[:, s] = Ym ./ Um
    end
end

# ╔═╡ 8f2ed482-5fb3-40d8-8f41-d2dfa8f9c23a
begin
    ps = []

    db(s) = amp2db.(abs.(s))
    for i = 1:length(As)
        p = plot(
            wo[1:2:end],
            db.([GD0[1:2:end]]);
            ylims = (-40, 20),
            xlims = (1, 55),
            legend = false,
            title = "Amplitude $(As[i])",
        )
        plot!(
            wo[1:2:end][30:end-1], # shifted because how `diff` works in julia
            db.(GD[:, i][1:end-30]);
        )
        plot!(fD[1:end-2], GMulti[:, i] .|> abs .|> amp2db)
        push!(ps, p)
    end
    push!(
        ps,
        scatter(
            (1:3)';
            ylims = (0, 0.00001),
            framestyle = :none,
            label = reshape(["G₀", "GD", "GMul"], 1, :),
        ),
    )
    plot(ps...; layout = (3, 2))
end

# ╔═╡ 4601db2f-610b-495b-8c60-4795d05b676b
begin
    ps2 = []
    phase(s) = rad2deg.(angle.(s))
    for i = 1:length(As)
        p = plot(
            wo[1:2:end],
            phase.([GD0[1:2:end]]);
            ylims = (-180, 180),
            xlims = (1, 55),
            legend = false,
            title = "Amplitude $(As[i])",
        )
        scatter!(
            wo[1:2:end][30:end-1], # shifted because how `diff` works in julia
            phase.(GD[:, i][1:end-30]);
            ms = 2,
        )
        scatter!(fD[1:end-2], phase(GMulti[:, i]); ms = 2)
        push!(ps2, p)
    end
    push!(
        ps2,
        scatter(
            (1:3)';
            ylims = (0, 0.00001),
            framestyle = :none,
            label = reshape(["G₀", "GD", "GMul"], 1, :),
        ),
    )
    plot(ps2...; layout = (3, 2))
end

# ╔═╡ Cell order:
# ╠═42c6fcac-b935-11ee-1e01-cde2ae4b39e8
# ╠═b12c3689-fbed-4803-b112-eb19c7113696
# ╠═6597d608-4a5d-4ef8-9fff-57b690333c0f
# ╟─4c3e735b-2cd4-47dd-a3fc-a823e5ab6176
# ╠═05820097-a0b2-4196-9d85-16501b9d9777
# ╟─8f2ed482-5fb3-40d8-8f41-d2dfa8f9c23a
# ╟─4601db2f-610b-495b-8c60-4795d05b676b

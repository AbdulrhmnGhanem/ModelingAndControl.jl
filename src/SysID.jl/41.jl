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
md"# Exercise 41: FRF measurement using burst excitation

!!! purpose
	- On exciting a system with random noise, the leakage effects are due to the non periodic nature of the random noise.
	- An alternative is to use a burst excitation, that starts at the beginning of the measurement window, and falls back to zero well before the end of the window. This makes it possible to include the transient effects completely in the window (the begin and end transients), thus eliminating all the leakage effects.
	- A white burst signal `burst(N, n)` has size `N` and the first `n` elements
	  are white noise.
"

# ╔═╡ 62c5feb4-250b-4341-9721-fd65f019f519
function cheby1(n, r, wp)
    h = digitalfilter(Lowpass(wp), Chebyshev1(n, r))
    tf = convert(PolynomialRatio, h)
    coefb(tf), coefa(tf)
end

# ╔═╡ 94c099c2-a626-11ee-3cc4-4504b89da241
begin
    N = 256
    N_burst = [64, 128, 192]
    fₛ = N

    b, a = cheby1(2, 20, 0.2)
    b *= 10  # increasing the filter gain

    t = (0:N-1) / fₛ
    Lines = 1:N÷2
    f = (Lines .- 1) / N * fₛ

    G₀, w = freqresp(PolynomialRatio(b, a))

    ps = []
    for n in N_burst
        u = zeros(N, 1)
        # a white burst signal burst(n, N) has size `N` and the first `n` elements
        # are white noise.
        u[1:n] = randn(n, 1)
        U = fft(u) / √N
        U = U[Lines]

        y = filt(b, a, u)
        Y = fft(y) / √N
        Y = Y[Lines]

        G = Y ./ U

        time_plot = plot(
            t,
            [u, y];
            ylims = (-5, 5),
            xlabel = "Time (s)",
            ylabel = n == N_burst[1] ? "Amplitude" : "",
            legend = false,
        )
        push!(ps, time_plot)
        spec_plot = plot(
            f,
            G₀[1:2:end] .|> abs .|> amp2db;
            xlims = (0, 60),
            ylims = (-40, 20),
            xlabel = "Frequency (Hz)",
            ylabel = n == N_burst[1] ? "Amplitude (dB)" : "",
            legend = false,
        )
        spec_plot = scatter!((G₀[1:2:end-1] - G) .|> abs .|> amp2db)
        push!(ps, spec_plot)
    end

end

# ╔═╡ 8ea13547-6b1a-4223-865c-8cf0e0d8764c
begin
    legend = scatter(
        (1:2)';
        ylims = (0, 0.00001),
        framestyle = :none,
        label = reshape(["u(t)", "y(t)"], 1, :),
    )
    plot(ps[1:2:end]..., legend; layout = (2, 3))
end

# ╔═╡ a20d9367-dbee-4586-854b-edfb1dbcc1fd
begin
    legend2 = scatter(
        (1:2)';
        ylims = (0, 0.00001),
        framestyle = :none,
        label = reshape(["G₀", "G₀ - G"], 1, :),
    )
    plot(ps[2:2:end]..., legend2; layout = (2, 3))
end

# ╔═╡ Cell order:
# ╠═62e757e9-5279-42c7-86ce-0dc0a96ff426
# ╠═cd3fdf9d-7285-4b24-96db-22edcf015835
# ╟─cd6ea107-5923-4239-a2fe-6da1a6dffcd8
# ╠═62c5feb4-250b-4341-9721-fd65f019f519
# ╠═94c099c2-a626-11ee-3cc4-4504b89da241
# ╟─8ea13547-6b1a-4223-865c-8cf0e0d8764c
# ╟─a20d9367-dbee-4586-854b-edfb1dbcc1fd

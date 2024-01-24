### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 3b4420fe-a053-11ee-1fdf-a9dd111cdac1
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 146a58f9-2848-414b-a87c-9eb15d1e3a08
using Plots, FFTW, DSP

# ╔═╡ fe2809da-37d0-4cd2-9143-28ec17b4b786
md"# Exercise 16: Generate a sine wave, noninteger number of periods

!!! purpose
	Leakage errors appears if the number of measured periods (depends on ω) is non interger.
"

# ╔═╡ 1afece58-a648-4285-95e8-1d7a68bec830
begin
    f = 80
    ω = 2π * f
    A = 1
    ϕ = 0.5π
    fₛ = 1000
    Tₛ = 1 / fₛ
    N = 16
	# from zero to N-1 because the Nth point belongs to the next period
	# it's an open interval [0, NTₛ[
    interval = 0:N-1

    u = map(t -> A * sin(ω * t * Tₛ + ϕ), interval)
    # we choose `1/N` as the scaling factor because the number of the components of the sine wave in the frequencey domain doesn't depend on `N`.
    U = fft(u) / N

    freqs = fₛ * (interval) / N  # how to transform DFT line numbers to frequencies.
    t = Tₛ * interval        # how to convert sampling points to time scale.
end;

# ╔═╡ f55f44f3-006b-4691-a50e-21c02731c49f
begin
    p1 = scatter(t, u; ylabel = "y(t)", xlabel = "Time (s)")

    p2 = scatter(
        abs.(U);
        line = :stem,
        marker = :circle,
        ylims = (0, 1),
        xlims = (0, 17),
        ylabel = "Amplitude (Linear)",
        xlabel = "DFT line number",
    )

    p3 = scatter(
        freqs,
        abs.(U);
        line = :stem,
        marker = :circle,
        ylims = (0, 1),
        xlims = (-70, 1000),
        ylabel = "Amplitude",
        xlabel = "Hz",
    )

    p4 = scatter(
        freqs,
        amp2db.(abs.(U));
        marker = :x,
        ylims = (-40, 0),
        xlims = (-70, 1000),
        ylabel = "Amplitude (dB)",
        xlabel = "Hz",
    )
    plot(p1, p4, p2, p3; legend = false)
end

# ╔═╡ Cell order:
# ╠═3b4420fe-a053-11ee-1fdf-a9dd111cdac1
# ╠═146a58f9-2848-414b-a87c-9eb15d1e3a08
# ╟─fe2809da-37d0-4cd2-9143-28ec17b4b786
# ╠═1afece58-a648-4285-95e8-1d7a68bec830
# ╟─f55f44f3-006b-4691-a50e-21c02731c49f

### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 47dcc054-a356-11ee-2f7c-91a157fb9189
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ c7860918-1e9f-4752-89f8-2a67d7b40ca2
using Plots, StatsBase, FFTW, DSP

# ╔═╡ 26680e4b-b217-44ea-b99d-099a51873feb
md"# Exercise 30: Smoothing the amplitude spectrum of a random excitation

!!! purpose
	- Averaging over multiple realizations allows to average the power spectrum.
	- The measured power spectrum ix Χ²-distributed with mean of Sᵤᵤ(ωₖ) and variance 2Sᵤᵤ(ωₖ) / M.
"

# ╔═╡ 9d5eb21e-f8fb-46bb-b572-c89204b474ae
begin
    N = 2^12
    Ms = [1, 16, 64]
    fₛ = 1000
    freqs = 1:N÷2
    f = (freqs .- 1) / N * fₛ

    Us = Vector(undef, 6)

    for i in eachindex(Ms)
        u = randn(N, Ms[i])  # one column is one realization
        U = fft(u, 1) / √N   # the `fft(u, 1) means apply fft along each columns
        U = mean(abs.(U[freqs, :]) .^ 2, dims=2) .|> sqrt

        Us[i] = scatter(f, U .|> amp2db;
            xlabel="Frequency (Hz)",
            ylabel="Amplitude (dB)",
            ylims=(-30, 20),
            legend=false,
            title="Number of averages $(Ms[i])"
        )

        Us[i+3] = stephist(U;
            legend=false,
            bins=:sqrt,
            xlims=(0, 4),
            ylims=(0, 200),
            xlabel="Amplitude (linear)",
        )
    end
end

# ╔═╡ 7dcaf04b-ddea-4ea7-80ba-7c72b1ab9c7e
plot(Us[1], Us[1+3], layout=(2, 1))

# ╔═╡ 3552b590-1be3-40b0-b421-71651aabb0bb
plot(Us[2], Us[1+3], layout=(2, 1))

# ╔═╡ 658ed415-9159-413c-b84f-189c51b32d7c
plot(Us[3], Us[3+3], layout=(2, 1))

# ╔═╡ Cell order:
# ╠═47dcc054-a356-11ee-2f7c-91a157fb9189
# ╟─26680e4b-b217-44ea-b99d-099a51873feb
# ╠═c7860918-1e9f-4752-89f8-2a67d7b40ca2
# ╠═9d5eb21e-f8fb-46bb-b572-c89204b474ae
# ╟─7dcaf04b-ddea-4ea7-80ba-7c72b1ab9c7e
# ╟─3552b590-1be3-40b0-b421-71651aabb0bb
# ╟─658ed415-9159-413c-b84f-189c51b32d7c

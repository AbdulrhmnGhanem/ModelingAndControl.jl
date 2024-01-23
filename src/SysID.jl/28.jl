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
using Plots, Statistics, FFTW, DSP

# ╔═╡ 26680e4b-b217-44ea-b99d-099a51873feb
md"# Exercise 28: Repeated realizations of a white random noise excitation with fixed lenght

!!! purpose
	White noise, is random noise with zero mean and finite variance. Its power specteral density is constant $G_x(f) = A$.

	Although we can see in this exercise that individual realizations is not flat. There are large spikes and dibs. This leads to poor SNR.
"

# ╔═╡ 9d5eb21e-f8fb-46bb-b572-c89204b474ae
begin
    N = 128
    fₛ = 1000
    interval = 1:N/2
    freqs = interval / N * fₛ

    Us = Vector(undef, 4)

    for i in eachindex(Us)
        u = randn(N) |> u -> u / std(u)
        U = fft(u) / √N .|> abs .|> amp2db
        Us[i] = scatter(
            freqs,
            U;
            xlabel = "Frequency (Hz)",
            ylabel = "Amplitude (dB)",
            ylims = (-30, 10),
            legend = false,
        )
    end
end

# ╔═╡ 7dcaf04b-ddea-4ea7-80ba-7c72b1ab9c7e
plot(Us...)

# ╔═╡ Cell order:
# ╠═47dcc054-a356-11ee-2f7c-91a157fb9189
# ╟─26680e4b-b217-44ea-b99d-099a51873feb
# ╠═c7860918-1e9f-4752-89f8-2a67d7b40ca2
# ╠═9d5eb21e-f8fb-46bb-b572-c89204b474ae
# ╟─7dcaf04b-ddea-4ea7-80ba-7c72b1ab9c7e

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
md"# Exercise 31: Generation of random noise excitations with a user imposed power spectrum

!!! purpose
	Shaping white noise power spectrum using user defined digital filter.
"

# ╔═╡ 337395e9-200d-434d-ab3d-dc2d7192551e
function butter(n, wn)
    digitalfilter(Lowpass(wn), Butterworth(n))
end

# ╔═╡ 9d5eb21e-f8fb-46bb-b572-c89204b474ae
begin
    N = 2^12
    Ms = [1, 16, 64]
    fₛ = 1000
    freqss = 1:N÷2
    f = (freqss .- 1) / N * fₛ
    h = butter(6, 0.1 * 2)

    Us = Vector(undef, 3)

    for i in eachindex(Ms)
        u = randn(N, Ms[i])  # one column is one realization
        u = filt(h, u)
        # the `fft(u, 1) means apply fft along each columns
        U_react = fft(u, 1) / √N
        U_react = mean(abs.(U_react[freqss, :]) .^ 2, dims = 2) .|> sqrt

        scale = mean(hanning(N) .^ 2) .|> sqrt
        U_hann = fft(hanning(N) .* u, 1) / √N
        U_hann = (mean(abs.(U_hann[freqss, :]) .^ 2, dims = 2) .|> sqrt) / scale

        scatter(
            f,
            U_react .|> amp2db;
            xlabel = "Frequency (Hz)",
            ylabel = "Amplitude (dB)",
            ylims = (-80, 20),
            # legend=false,
            label = "Rectangular",
            title = "Number of averages $(Ms[i])",
            legendtitle = "window",
        )

        Us[i] = scatter!(f, U_hann .|> amp2db; label = "Hanning")
    end
end

# ╔═╡ 486b78b4-b0c5-49b0-9b07-6f996c11b228
plot(Us[1])

# ╔═╡ a2cb3f1f-3260-4259-ba0e-3b3b0af9346a
plot(Us[2])

# ╔═╡ 4ec4c6fe-5a04-4c92-8ca3-e4b00944f635
plot(Us[3])

# ╔═╡ Cell order:
# ╠═47dcc054-a356-11ee-2f7c-91a157fb9189
# ╟─26680e4b-b217-44ea-b99d-099a51873feb
# ╠═c7860918-1e9f-4752-89f8-2a67d7b40ca2
# ╠═337395e9-200d-434d-ab3d-dc2d7192551e
# ╠═9d5eb21e-f8fb-46bb-b572-c89204b474ae
# ╟─486b78b4-b0c5-49b0-9b07-6f996c11b228
# ╟─a2cb3f1f-3260-4259-ba0e-3b3b0af9346a
# ╟─4ec4c6fe-5a04-4c92-8ca3-e4b00944f635

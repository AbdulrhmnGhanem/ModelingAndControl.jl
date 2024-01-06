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
md"# Exercise 42: Impulse response function measurements in the presence of output noise

!!! purpose
	Increasing the record length does not result in a proportional increase of information about the system because the noise dominates completely at the end of
	the window. After a while, only noise is added to the record which will deteriorate the SNR.
"

# ╔═╡ 0345624c-b9f8-4a13-9371-f08f428523ce
function cheby1(n, r, wp)
    digitalfilter(Lowpass(wp), Chebyshev1(n, r))
end

# ╔═╡ 3d1a9815-e856-4c2d-a7ab-2ff2d5daadd2
begin
	N = 128
	Nᵣ = 100
    fₛ = 256
    t = (0:N-1) / fₛ
    freqs = (0:N-1) / N * fₛ
    freqs_lines = 1:N÷2
    f_cuttoff = 0.1fₛ
    order = 2
    reseonance = 10

    h = cheby1(order, reseonance, 2f_cuttoff / fₛ)
    # Calculate FRF directly.
    H, w = freqresp(h)

    # Calculate FRF using impulse reseponse
    u = [1; zeros(N - 1)]
    y₀ = filt(h, u)
	ỹₙ = Normal(0, 0.02)
	y = zeros(N, Nᵣ)
	
	for n in 1:Nᵣ
		y[:, n] = y₀ + rand(ỹₙ, size(y₀))
	end

    U = fft(u) / N
    Y₀ = fft(y₀)
	Y = fft(y, 2)
end;

# ╔═╡ 13ff628a-5158-4a6e-be0d-9e07defe5757
begin
    p1 = scatter(t, u;
        line=:stem,
        marker=:circle,
        title="Time domain",
        ylabel="Amplitude",
        legend=false,
        ylims=(0, 1.5),
    )
    p2 = plot(t, [y₀, rand(ỹₙ, size(y₀))];
        xlabel="Time (s)",
        ylabel="Amplitude",
        legend=false,
        ylims=(-0.2, 0.2),
    )
    plot(p1, p2; layout=(2, 1))
end

# ╔═╡ de29e0ad-7a06-40dd-ab69-531496071a81
begin
	abs_db(s) = abs.(s) .|> amp2db
	 
	plot(freqs[freqs_lines], abs_db(Y₀[freqs_lines]);
	        xlabel="Frequency (Hz)",
	        ylabel="Amplitude",
			ylims=(-60, 0),
			label="Y₀",
			title="Single experiment",
	    )
	scatter!(freqs[freqs_lines], 
			abs_db.([Y₀[freqs_lines] - Y[end, :][freqs_lines], Y[end, :][freqs_lines]]);
			labels=reshape(["Y₀ - Y", "Y"], 1,:)
	)
end

# ╔═╡ a987c10e-03cc-49f7-ba6b-b5542454e6ea
begin
	Ȳ = mean(Y, dims=1)
	σ = std(Y, dims=1)
	plot(freqs[freqs_lines], abs_db(Y₀[freqs_lines]);
	        xlabel="Frequency (Hz)",
	        ylabel="Amplitude",
			ylims=(-60, 0),
			label="Y₀",
			title="Averaged output spectrum",
	    )
	scatter!(freqs[freqs_lines], 
			abs_db.([Y₀[freqs_lines] - Ȳ[freqs_lines], Ȳ[freqs_lines], σ[freqs_lines] / √Nᵣ]);
			labels=reshape(["Y₀ - Y", "Y", "σ"], 1,:)
	)
end

# ╔═╡ Cell order:
# ╠═dc2d1942-accd-11ee-34f8-ef2f689c65b2
# ╠═6ce6c76b-c24b-4a6b-8297-f62cdce1191e
# ╟─40338bfb-3c0d-4990-b84b-0ca6d42373c0
# ╠═0345624c-b9f8-4a13-9371-f08f428523ce
# ╠═3d1a9815-e856-4c2d-a7ab-2ff2d5daadd2
# ╟─13ff628a-5158-4a6e-be0d-9e07defe5757
# ╟─de29e0ad-7a06-40dd-ab69-531496071a81
# ╟─a987c10e-03cc-49f7-ba6b-b5542454e6ea

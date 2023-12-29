### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 2707e0f5-bccf-4089-8856-b7a6915403eb
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ fcc0a096-006a-411f-b0ff-af7362fdedfe
using Plots, DSP, FFTW

# ╔═╡ 193f46d0-7eb2-451b-be9e-fff2d25a4c26
md"# Exercise 34: Impulse response function measurements"

# ╔═╡ 1f06e3d0-8ed7-4016-ad08-fe0411e394ea
function cheby1(n, r, wp)
    h = digitalfilter(Lowpass(wp), Chebyshev1(n, r))
    convert(PolynomialRatio, h)
end

# ╔═╡ 8a26f77e-4815-4760-8deb-62b7ef10e299
begin
    N = 128
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
    y = filt(h, u)

    U = fft(u) / √N
    Y = fft(y) / √N
end;

# ╔═╡ 6b38176b-15d3-4a6d-8162-952e74012904
plot(w, abs.(H);
    title="FRF",
    legend=false,
    xlabel="Normalized Frequency",
    ylabel="Amplitude",
    ylims=(0, 1.5),
)

# ╔═╡ 2df0cc31-f8a2-4627-a4d9-82cddf02ec21
begin
    p1 = scatter(t, u;
        line=:stem,
        marker=:circle,
        title="Time domain",
        ylabel="Amplitude",
        legend=false,
        ylims=(0, 1.5),
    )
    p2 = plot(t, y;
        xlabel="Time (s)",
        ylabel="Amplitude",
        legend=false,
        ylims=(-0.2, 0.2),
    )
    plot(p1, p2; layout=(2, 1))
end

# ╔═╡ c3acd5df-7e84-4b35-a3f3-48f8f84c937b
begin
    p3 = plot(freqs[freqs_lines], abs.(U[freqs_lines]);
        title="Frequency domain",
        ylabel="Amplitude",
        legend=false,
        ylims=(0, 0.12),
    )
    p4 = plot(freqs[freqs_lines], abs.(Y[freqs_lines] / U[1]); # why are we dividing by U[1]?
        xlabel="Frequency (Hz)",
        ylabel="Amplitude",
        legend=false,
        ylims=(0, 1.5),
    )
    plot(p3, p4; layout=(2, 1))
end

# ╔═╡ Cell order:
# ╠═2707e0f5-bccf-4089-8856-b7a6915403eb
# ╠═fcc0a096-006a-411f-b0ff-af7362fdedfe
# ╟─193f46d0-7eb2-451b-be9e-fff2d25a4c26
# ╠═1f06e3d0-8ed7-4016-ad08-fe0411e394ea
# ╠═8a26f77e-4815-4760-8deb-62b7ef10e299
# ╟─6b38176b-15d3-4a6d-8162-952e74012904
# ╟─2df0cc31-f8a2-4627-a4d9-82cddf02ec21
# ╟─c3acd5df-7e84-4b35-a3f3-48f8f84c937b

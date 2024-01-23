### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 475cda4c-fcf1-41e9-a810-966f87dcfaa6
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 08b06706-8ca8-4ef3-852d-9146228351a7
using Plots, FFTW, DSP

# ╔═╡ c66e7fca-0a21-4c24-81e0-466c73c40f50
md"""# Exercise 21: the swept sine signal
!!! purpose
	Generate a sweeping siganl using sweeping sine. Compared to the Schroeder multisine, not all power is in the frequency band of interest.
"""

# ╔═╡ 13fea8c0-e2a1-4525-9ab9-71573ccb98ce
begin
    f₀ = 1  # a period of one second
    fₛ = 1000
    N = 1000
    t = (0:N-1) / fₛ
    k₁ = 50   # the lower bound of the desired power spectera
    k₂ = 200  # the upper bound of the desired power spectera
    a = π * (k₂ - k₁) * f₀
    b = 2π * k₁ * f₀
    u = sin.((a .* t .+ b) .* t)
    U = amp2db.(abs.(fft(u) / sqrt(N)))
end;

# ╔═╡ d930593f-b3a7-4825-ae4e-f20f53ec9d52
plot(
    t,
    [u, U];
    xlabel = reshape(["Time (s)", "Frequency (Hz)"], 1, :),
    ylabel = reshape(["Amplitude", "Amplitude (dB)"], 1, :),
    layout = (2, 1),
    legend = false,
)

# ╔═╡ 3c823b94-12aa-4d2b-8a0a-7ddc3a9ff8d3
md"""
# Crest factor

cr(u) = $(maximum(u) / rms(u))
"""

# ╔═╡ Cell order:
# ╠═475cda4c-fcf1-41e9-a810-966f87dcfaa6
# ╠═08b06706-8ca8-4ef3-852d-9146228351a7
# ╟─c66e7fca-0a21-4c24-81e0-466c73c40f50
# ╠═13fea8c0-e2a1-4525-9ab9-71573ccb98ce
# ╟─d930593f-b3a7-4825-ae4e-f20f53ec9d52
# ╟─3c823b94-12aa-4d2b-8a0a-7ddc3a9ff8d3

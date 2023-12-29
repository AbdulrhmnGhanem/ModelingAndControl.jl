### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ 4658abb8-a51c-11ee-094d-830c5934815a
# ╠═╡ show_logs = false
begin 
 # If you are running this notebook as a stannalone notebook disable this cell.
 import Pkg 
 Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ f4535fab-50f6-4a5b-babe-43cfdc14158a
using Plots, DSP, FFTW, PlutoUI

# ╔═╡ e132bd52-784f-41b7-938a-a5f51f829328
md"# Exercise 35: Study of the sine response of a linear system (transient and steady state)

!!! purpose
	- It can be concluded that the steady-state is reached, only when teh transient is sufficiently samll at the end of the measurement interval.

	- Increasing the record length, decreases amplitude of the leakage skirt (you can verify that by adjusting the slider on top of second figure).
"

# ╔═╡ c3b41dde-1544-49ca-93c4-c43068b7a6da
function cheby1(n, r, wp)
    h = digitalfilter(Lowpass(wp), Chebyshev1(n, r))
    convert(PolynomialRatio, h)
end

# ╔═╡ 2a1ad35f-5c55-47a5-bb6f-359dd8d5e6a9
@bind record_length_multiplier Slider(1:2:32)

# ╔═╡ a570dbcf-0d97-426d-a45d-b66c0918ff72
begin
	N = 128 * record_length_multiplier
    fₛ = 256
    t = (0:N-1) / fₛ
    freqs = (0:N-1) / N * fₛ
    freqs_lines = 1:N÷2
    f_cuttoff = 0.1fₛ
    order = 2
    reseonance = 10
	f = 16
	ω = 2π * f
	periods = N ÷ f
    
	h = cheby1(order, reseonance, 2f_cuttoff / fₛ)
	# everything so far is exactly like the previous exercise.
	# Now, we excite the system with a cosine instead of impulse signal.
	u = cos.(ω .* t)
	y = filt(h, u)

	# we consider that all transient effects have vanished in the final period.
	yₛₛ = repeat(y[end-f+1:end], periods)

	# the spectrum
	U = fft(u) / N .|> abs .|> amp2db
	Y = fft(y) / N .|> abs .|> amp2db
	Yₛₛ = fft(yₛₛ) / N .|> abs .|> amp2db
end;

# ╔═╡ 06db396f-bc53-4211-9a80-750f60d7d71b
plot(t, [y yₛₛ];
	xlabel="Time (s)",
	ylabel="Amplitude",
	labels=reshape(["y", "yₛₛ"], 1, :),
	xlims=(0, 0.25),
)

# ╔═╡ a0f57a25-381e-42f9-b0cc-8863e54029fe
scatter(freqs[freqs_lines], [U[freqs_lines], Y[freqs_lines], Yₛₛ[freqs_lines]];
	xlabel="Frequency (Hz)",
	ylabel="Amplitude (dB)",
	ylims=(-60, 0),
	legend=:topright,
	labels=reshape(["u", "y", "yₛₛ",], 1, :),
)

# ╔═╡ Cell order:
# ╠═4658abb8-a51c-11ee-094d-830c5934815a
# ╠═f4535fab-50f6-4a5b-babe-43cfdc14158a
# ╟─e132bd52-784f-41b7-938a-a5f51f829328
# ╠═c3b41dde-1544-49ca-93c4-c43068b7a6da
# ╠═a570dbcf-0d97-426d-a45d-b66c0918ff72
# ╟─06db396f-bc53-4211-9a80-750f60d7d71b
# ╠═2a1ad35f-5c55-47a5-bb6f-359dd8d5e6a9
# ╟─a0f57a25-381e-42f9-b0cc-8863e54029fe

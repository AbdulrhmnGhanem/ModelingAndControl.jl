### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 89bb7956-a37f-11ee-1ca1-3305b8e087f1
# ╠═╡ show_logs = false
begin 
 # If you are running this notebook as a stannalone notebook disable this cell.
 import Pkg 
 Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 1a2a6937-e9fc-4c3d-a771-dfe6d6d81c38
using Plots, DSP

# ╔═╡ d453b6b9-1edd-48c0-b6c4-16ceb08f9763
md"# Exercise 32: Amplitude distribution of filtered noise

!!! purpose
	Show that filtering the noise results in a signal with noramlly distributed amplitude. This is due to ceneral limit theorem.
"

# ╔═╡ 00c56ca8-6287-4a12-8fb9-02c3a8e86ffa
function butter(n, wn)
    digitalfilter(Lowpass(wn), Butterworth(n))
end

# ╔═╡ c427462b-621b-4727-a20d-76cf2f8ca5bd
begin
	N = 100_000
	u = sign.(randn(N))
	h = butter(6, 0.1 * 2)
	y = filt(h, u)
end;

# ╔═╡ 02931245-5645-470b-97e6-12d13980e9a6
histogram(u;
	legend=false,
	xlabel="Amplitude",
	ylabel="Count",
	title="Binary input",
)

# ╔═╡ 65609025-db31-4ec4-b8da-feb3dd50fe4e
histogram(y;
	legend=false,
	xlabel="Amplitude",
	ylabel="Count",
	title="Filtered binary input",
)

# ╔═╡ Cell order:
# ╠═89bb7956-a37f-11ee-1ca1-3305b8e087f1
# ╠═1a2a6937-e9fc-4c3d-a771-dfe6d6d81c38
# ╟─d453b6b9-1edd-48c0-b6c4-16ceb08f9763
# ╠═00c56ca8-6287-4a12-8fb9-02c3a8e86ffa
# ╠═c427462b-621b-4727-a20d-76cf2f8ca5bd
# ╟─02931245-5645-470b-97e6-12d13980e9a6
# ╟─65609025-db31-4ec4-b8da-feb3dd50fe4e

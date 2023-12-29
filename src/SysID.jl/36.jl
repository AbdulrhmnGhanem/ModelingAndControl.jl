### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 6b94407e-3945-4fe5-b654-f178fff9c300
# ╠═╡ show_logs = false
begin 
 # If you are running this notebook as a stannalone notebook disable this cell.
 import Pkg 
 Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ e1197ba0-13a3-4d39-8da7-a9b9a3fb18ff
using Plots, DSP, FFTW, Distributions

# ╔═╡ 9fd5daa6-4493-45da-9487-544f4c4d0d28
md"# Exercise 36: Study of multisine response of a linear system (transient and steady state)

!!! purpose
	Exciting the system persistenly improves the SNR.
"

# ╔═╡ 334d4bdf-a5e2-4c85-9ce2-87ed45671105
function cheby1(n, r, wp)
    h = digitalfilter(Lowpass(wp), Chebyshev1(n, r))
    convert(PolynomialRatio, h)
end

# ╔═╡ 4714a432-9377-4a6d-b99f-65b7d9c8bc90
begin
	N = 128
	M = 3 # periods of excitation
    fₛ = 256
    t = (0:M*N-1) / fₛ
    freqs = (0:N-1) / N * fₛ
    freqs_lines = 1:N÷2
    f_cuttoff = 0.1fₛ
    order = 2
    reseonance = 10
	
	h = cheby1(order, reseonance, 2f_cuttoff / fₛ)
	# everything so far is exactly like the previous exercise.
	# Now, we excite the system with a multisine instead of cosine signal.
	F = 100  # Exciting the system up to 100 Hz
	ϕ_rand = rand(Uniform(0, 2π), N)
	U = zeros(Complex{Float64}, N, 1)
	U[2:F÷2-1] = exp.(im * ϕ_rand[2:F÷2-1])
	u = real(ifft(U))
	u = u / std(u)
	u = repeat(u, M)
	y = filt(h, u)

	# we consider that all transient effects have vanished in the final period.
	yₛₛ = repeat(y[end-N+1:end], M)
	yₜᵣ = y - yₛₛ
end;

# ╔═╡ e3685bb2-ea31-42e7-8b99-f306da900cad
plot(t, y;
	title="Output",
	xlabel="Time (s)",
	ylabel="y(t)",
	legend=false,
)

# ╔═╡ 5e4b1315-4fed-4856-8ebd-f03d43e211a9
plot(t, yₜᵣ;
	xlabel="Time (s)",
	ylabel="yₜᵣ(t)",
	legend=false,
)

# ╔═╡ c200d759-66c7-47e3-b56a-1060854dbbbe
plot(t[begin:255], abs.(yₜᵣ)[begin:255];
	xlabel="Time (s)",
	ylabel="yₜᵣ(t)",
	legend=false,
	yaxis=:log10,
	xlims=(0, 1.5),
)

# ╔═╡ Cell order:
# ╠═6b94407e-3945-4fe5-b654-f178fff9c300
# ╠═e1197ba0-13a3-4d39-8da7-a9b9a3fb18ff
# ╟─9fd5daa6-4493-45da-9487-544f4c4d0d28
# ╠═334d4bdf-a5e2-4c85-9ce2-87ed45671105
# ╠═4714a432-9377-4a6d-b99f-65b7d9c8bc90
# ╟─e3685bb2-ea31-42e7-8b99-f306da900cad
# ╠═5e4b1315-4fed-4856-8ebd-f03d43e211a9
# ╟─c200d759-66c7-47e3-b56a-1060854dbbbe

### A Pluto.jl notebook ###
# v0.19.27

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
md"# Exercise 19: Generate a sine wave using `ifft`

!!! purpose
	For a sum of harmoically related sines with many components, using `ifft` reduces the calculation time.
"

# ╔═╡ 1afece58-a648-4285-95e8-1d7a68bec830
begin
	A = 1
	ϕ = 0
	fₛ = 1000
	Tₛ = 1/fₛ
	N = 16
	f = fₛ / N
	ω = 2π * f
	interval = 0:N-1
	
	U₁ = zeros(Complex{Float64}, N, 1)
	U₁[2] = exp(-0.5π * im)
	U₁[N] = exp(0.5π * im)
	u₁ = ifft(U₁) * 0.5N

	U₂ = zeros(Complex{Float64}, N, 1)
	U₂[2] = exp(-0.5π * im)
	# here we are mutliplyin by `N` instead of `N/2` because we're using a single frequencey component
	u₂ = ifft(U₂) * N
	
	t = Tₛ * interval        # how to convert sampling points to time scale.
end;

# ╔═╡ f55f44f3-006b-4691-a50e-21c02731c49f
begin
	scatter(t, real.(u₁);
		ylabel="y(t)", 
		xlabel="Time (s)",
		label="Two components"
	)
	scatter!(t, real.(u₂);
		marker=:x,
		label="One component",
	)
end

# ╔═╡ Cell order:
# ╠═3b4420fe-a053-11ee-1fdf-a9dd111cdac1
# ╠═146a58f9-2848-414b-a87c-9eb15d1e3a08
# ╟─fe2809da-37d0-4cd2-9143-28ec17b4b786
# ╠═1afece58-a648-4285-95e8-1d7a68bec830
# ╟─f55f44f3-006b-4691-a50e-21c02731c49f

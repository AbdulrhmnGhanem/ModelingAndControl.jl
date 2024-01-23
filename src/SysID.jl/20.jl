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
using Plots, FFTW, DSP, StatsBase, Distributions

# ╔═╡ fe2809da-37d0-4cd2-9143-28ec17b4b786
md"# Exercise 20: Generation of a multisine with flat amplitude spectrum

!!! purpose
	Generate multisines signals with three different pahse shifts: `linear pahse (impulse)`, `random phase`, `Schroeder`.

$u(t) = \sum{}_{k = 1}^{F} A cos(2\pi f_0 k t T_s + \phi_r)$
`F`: number of sines, we select it depending on the frequencey band, in this exercise the frequencey band is `[0, 0.1fₛ]` so `F = 0.1fₛ`.


"

# ╔═╡ 1afece58-a648-4285-95e8-1d7a68bec830
begin
    A = 1
    ϕ = 0.5π
    fₛ = 1000
    Tₛ = 1 / fₛ
    N = 1000
    f = fₛ / N
    ω = 2π * f
    interval = 0:N-1

    F = Int(0.1N)
    freqs = fₛ * interval / N

    # impluse multisine
    τ = 0.3
    ϕᵣ_i = -τ * 2π * freqs
    Uᵢ = zeros(Complex{Float64}, N, 1)
    Uᵢ[2:F+1] = exp.(im * ϕᵣ_i[2:F+1])
    uᵢ = 2real(ifft(Uᵢ))
    uᵢ = uᵢ / std(uᵢ)

    # random phase multisine
    ϕ_rand = rand(Uniform(0, 2π), N)
    Uᵣ = zeros(Complex{Float64}, N, 1)
    Uᵣ[2:F+1] = exp.(im * ϕ_rand[2:F+1])
    uᵣ = 2real(ifft(Uᵣ))
    uᵣ = uᵣ / std(uᵣ)

    # shroeder multisine
    k = 1:F
    ϕ_s = -π / F * k .* (k .- 1)
    Uₛ = zeros(Complex{Float64}, N, 1)
    Uₛ[2:F+1] = exp.(im * ϕ_s)  # this is a different way from others, why?
    uₛ = 2real(ifft(Uₛ))
    uₛ = uₛ / std(uₛ)
end;

# ╔═╡ cf2c147f-8971-4990-b7ea-1f3b6a2dc71b
begin
    ylims = (-5, 15)
    p1 = plot(
        interval,
        [uᵣ, uᵢ];
        ylims,
        labels = reshape(["random ϕ", "impluse ϕ"], 1, :),
        xlabel = "Time (ms)",
        ylabel = "Amplitude",
    )
    p2 = plot(
        interval,
        [uᵣ, uₛ];
        ylims,
        labels = reshape(["random ϕ", "schroeder ϕ"], 1, :),
        xlabel = "Time (ms)",
    )
    plot(p1, p2)
end

# ╔═╡ ee42a2a3-e8a8-4f42-afc7-d334eecbacf1
md"## Crest factors"

# ╔═╡ b43a0963-4418-4d8a-a747-a8401e11bbce
begin
    crest(u) = maximum(u) / rms(u)

    crᵢ = crest(uᵢ)
    crᵣ = crest(uᵣ)
    crₛ = crest(uₛ)

    md"""
    Cr(uᵢ) = $crᵢ, 

    Cr(uᵣ) = $crᵣ

    Cr(uₛ) = $crₛ
    """
end

# ╔═╡ Cell order:
# ╠═3b4420fe-a053-11ee-1fdf-a9dd111cdac1
# ╠═146a58f9-2848-414b-a87c-9eb15d1e3a08
# ╟─fe2809da-37d0-4cd2-9143-28ec17b4b786
# ╠═1afece58-a648-4285-95e8-1d7a68bec830
# ╟─cf2c147f-8971-4990-b7ea-1f3b6a2dc71b
# ╟─ee42a2a3-e8a8-4f42-afc7-d334eecbacf1
# ╟─b43a0963-4418-4d8a-a747-a8401e11bbce

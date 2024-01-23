### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 0128b83c-986d-11ee-38be-2fd8b3ac03a0
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ 4b1d5b72-2d3b-4c7e-a456-1beace344605
using Random, Plots, DSP, StatsBase

# ╔═╡ ea5e2791-aef5-42ed-a3a1-57c5a9afd616
md"# Exercise 12.a: The effect of filtering input noise with (varying cutoff frequency)

!!! purpose
	Show the effects of white noise on the measurements, and how filtering can introduce a bias in the identified system.

	* On removing current noise `nᵢ` the bais disappears for all configurations.
	* Removing the voltage noise doesn't affect the bias in these configurations.

	LS → is baised due to the noise on the input `nᵢ` (filtering the input doesn't have an effect when using LS).

	IV → is more complicated. For the white noise situation, no bias is visible. However, once the output noise is filtered, a bias becomes visible. The relative bias is proportional to the ratio of the autocorrelation functions of the noise and the current `Rₙᵢₙᵢ(s)/Rᵢ₀ᵢ₀(s)` 👇️

	![Textbook](https://i.imgur.com/Dcovnk3.png)
"

# ╔═╡ 2432fede-195f-4f21-9c56-2bedc88eecff
function butter(n, wn)
    h = digitalfilter(Lowpass(wn), Butterworth(n))
    tf = convert(PolynomialRatio, h)
    coefb(tf), coefa(tf)
end

# ╔═╡ 12f567b2-1eb7-44ff-be4d-0f47c54e7a51
begin
    N = 5000
    Nᵣ = 10_000
    R₀ = 1000
    iₘₐₓ = 0.01
    f_gen = 0.05
    fₙ = [0.4995, 0.475, 0.3]
    Nₜᵣₐₙₛ = 1000

    b_gen, a_gen = butter(1, 2f_gen)
    i₀ =
        randn(N + Nₜᵣₐₙₛ, 1) |>
        prev ->
            filt(b_gen, a_gen, prev) |>
            prev -> prev[Nₜᵣₐₙₛ+1:end] |> prev -> prev * iₘₐₓ / std(prev)

    u₀ = R₀ * i₀

    lag = 1

    LS = zeros(Nᵣ, length(fₙ))
    IV = zeros(Nᵣ, length(fₙ))

    for r = 1:length(fₙ)
        b, a = butter(2, 2fₙ[r])

        for s = 1:Nᵣ
            nᵤ = randn(N)
            nᵢ =
                randn(N + Nₜᵣₐₙₛ) |>
                prev ->
                    filt(b, a, prev) |>
                    prev -> prev[Nₜᵣₐₙₛ+1:end] |> prev -> prev / std(prev) * iₘₐₓ

            i = i₀ + nᵢ
            u = u₀ + nᵤ

            LS[s, r] = i \ u

            iShift = copy(i)
            deleteat!(iShift, 1:lag)
            deleteat!(i, N-lag+1:N)
            deleteat!(u, N-lag+1:N)

            IV[s, r] = (u' * iShift) / (iShift' * i)
        end

    end
end

# ╔═╡ 8256a5bc-9706-4a1c-b723-6937d81eb803
stephist(
    [LS IV];
    normalize = :pdf,
    xlims = (0, 1550),
    ylims = (0, 0.101),
    xlabel = "R (Ω)",
    ylabel = "PDF",
    labels = reshape(
        [
            "LS (fₙ = $(fₙ[1]))",
            "LS (fₙ = $(fₙ[2]))",
            "LS (fₙ = $(fₙ[3]))",
            "IV (fₙ = $(fₙ[1]))",
            "IV (fₙ = $(fₙ[2]))",
            "IV (fₙ = $(fₙ[3]))",
        ],
        1,
        :,
    ),
)

# ╔═╡ Cell order:
# ╠═0128b83c-986d-11ee-38be-2fd8b3ac03a0
# ╠═4b1d5b72-2d3b-4c7e-a456-1beace344605
# ╟─ea5e2791-aef5-42ed-a3a1-57c5a9afd616
# ╟─2432fede-195f-4f21-9c56-2bedc88eecff
# ╠═12f567b2-1eb7-44ff-be4d-0f47c54e7a51
# ╟─8256a5bc-9706-4a1c-b723-6937d81eb803

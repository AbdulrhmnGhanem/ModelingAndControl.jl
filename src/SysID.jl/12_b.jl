### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 2f54993a-9bac-11ee-392a-735f95358d60
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ d92b47cc-60e8-4cb7-b5da-948c52dff469
using Random, Plots, DSP, StatsBase

# ╔═╡ 5d8eb239-e25d-49d8-9c78-0062e11da2bf
md"# Exercise 12.b: The effect of filtering input noise with (varying IV lag)

!!! purpose
	See Exercise `12.a`.

	LS → is baised due to the noise on the input `nᵢ` (filtering the input doesn't have an effect when using LS).

	IV → The bias becomes smaller with increasing the lag, but the std increases with it. **The IV works well if the bandwidth of the generator signal is much smaller than the noise distribution**.
"

# ╔═╡ 20b9d11e-58a5-428c-af6b-820d5fec5705
function butter(n, wn)
    h = digitalfilter(Lowpass(wn), Butterworth(n))
    tf = convert(PolynomialRatio, h)
    coefb(tf), coefa(tf)
end

# ╔═╡ e0ee4985-a0b9-4831-aee4-a83c330ef406
begin
    N = 5000
    Nᵣ = 10_000
    R₀ = 1000
    iₘₐₓ = 0.01
    f_gen = 0.05
    fₙ = 0.3
    Nₜᵣₐₙₛ = 1000

    b_gen, a_gen = butter(1, 2f_gen)
    i₀ =
        randn(N + Nₜᵣₐₙₛ, 1) |>
        prev ->
            filt(b_gen, a_gen, prev) |>
            prev -> prev[Nₜᵣₐₙₛ+1:end] |> prev -> prev * iₘₐₓ / std(prev)

    u₀ = R₀ * i₀

    lags = [1, 2, 5]

    LS = zeros(Nᵣ, length(lags))
    IV = zeros(Nᵣ, length(lags))

    for r = 1:length(lags)
        b, a = butter(2, 2fₙ)

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
            deleteat!(iShift, 1:lags[r])
            deleteat!(i, N-lags[r]+1:N)
            deleteat!(u, N-lags[r]+1:N)

            IV[s, r] = (u' * iShift) / (iShift' * i)
        end
    end
end

# ╔═╡ 61939ecb-79a9-4471-96be-7354c31a30d5
stephist(
    [LS IV];
    normalize = :pdf,
    xlims = (0, 1550),
    ylims = (0, 0.101),
    xlabel = "R (Ω)",
    ylabel = "PDF",
    labels = reshape(
        [
            "LS (s = $(lags[1]))",
            "LS (s = $(lags[2]))",
            "LS (s = $(lags[3]))",
            "IV (s = $(lags[1]))",
            "IV (s = $(lags[2]))",
            "IV (s = $(lags[3]))",
        ],
        1,
        :,
    ),
)

# ╔═╡ Cell order:
# ╟─2f54993a-9bac-11ee-392a-735f95358d60
# ╠═d92b47cc-60e8-4cb7-b5da-948c52dff469
# ╟─5d8eb239-e25d-49d8-9c78-0062e11da2bf
# ╟─20b9d11e-58a5-428c-af6b-820d5fec5705
# ╠═e0ee4985-a0b9-4831-aee4-a83c330ef406
# ╟─61939ecb-79a9-4471-96be-7354c31a30d5

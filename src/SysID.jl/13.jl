### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 233379ee-9b92-11ee-09b4-550db5f44a98
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ c94ced4d-360b-4978-af40-d993fd1dfc41
using Distributions, Plots

# ╔═╡ 79646165-904c-498c-9e3e-0acec016a2d7
md"# Exercise 13: errors-in-variables method

!!! purpose
	Use EiV methods to estimate, the drawback of this method is the covariance matrix of the input and output must be known.
"

# ╔═╡ 5e7b9196-0a3c-44cb-adf9-58cedc104d68
begin
    Nᵣ = 100000
    N = 5000
    R₀ = 1000
    ĩ₀ = Normal(0, 0.01)
    ñᵢ = Normal(0, 0.001)
    ñᵤ = Normal(0, 1)

    i₀ = rand(ĩ₀, N)
    u₀ = i₀ * R₀

    LS = zeros(Nᵣ)
    EiV = zeros(Nᵣ)
    IV = zeros(Nᵣ)
    lag = 1

    Threads.@threads for r = 1:Nᵣ
        i = i₀ + rand(ñᵢ, N)
        u = u₀ + rand(ñᵤ, N)

        LS[r] = i \ u

        z1 = (u' * u / ñᵤ.:σ^2 - i' * i / ñᵢ.:σ^2)
        z2 = u' * i / ñᵤ.:σ^2
        EiV[r] = (z1 + sqrt(z1^2 + 4z2^2 / ñᵢ.:σ^2)) / 2z2

        iShift = copy(i)
        deleteat!(iShift, 1:lag)
        deleteat!(i, N-lag+1:N)
        deleteat!(u, N-lag+1:N)

        IV[r] = (u' * iShift) / (iShift' * i)

    end
end

# ╔═╡ a395b0c6-c350-472a-a47a-755fb52a77e8
stephist(
    [LS EiV];
    normalize = :pdf,
    xlabel = "R (Ω)",
    ylabel = "PDF",
    labels = reshape(["LS", "EiV"], 1, :),
)

# ╔═╡ cd8f3fc7-a493-49d7-b9d2-eabfcf15b982
md"## The IV method failure

The IV methid fails for this case because ACF for the current is too small so that results in a huge std

σ(IV) = $(round(std(IV); digits=3))

σ(EiV) = $(round(std(EiV); digits=3))

σ(LS) = $(round(std(LS); digits=3))
"

# ╔═╡ Cell order:
# ╟─233379ee-9b92-11ee-09b4-550db5f44a98
# ╠═c94ced4d-360b-4978-af40-d993fd1dfc41
# ╟─79646165-904c-498c-9e3e-0acec016a2d7
# ╠═5e7b9196-0a3c-44cb-adf9-58cedc104d68
# ╟─a395b0c6-c350-472a-a47a-755fb52a77e8
# ╟─cd8f3fc7-a493-49d7-b9d2-eabfcf15b982

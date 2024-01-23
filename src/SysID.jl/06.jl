### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 6c9c5140-969d-11ee-30b0-5d1f501a232b
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ aa954823-05a5-4343-a78d-10f034dddf56
using LinearAlgebra, Plots, Polynomials

# ╔═╡ 8a085c11-88d5-43e1-807a-865a92315440
md"# Exercise 6: least squares estimation of models that are linear in the parameters"

# ╔═╡ 92c8f99e-b3ea-430f-9f60-1ca9f0d4a6a6
begin
    y(u) = tan(u * 0.9 * pi / 2)
    u(n) = LinRange(0, 1, n)
    n = range(1, 20)
    N = 100
    u0 = u(N)
    y0 = y.(u0)
    KAll = zeros((N, length(n)))  # what is KAll???
    for k in n
        KAll[:, k] = u0 .^ (k - 1)
    end
    conds = zeros((length(n), 2))
    rmss = zeros((length(n), 2))

    for k in n
        K = KAll[:, 1:k]
        conds[k, 2] = cond(K)
        conds[k, 1] = cond(K' * K)

        # numerically unstable estimation
        estimator1 = inv(K'K) * K' * y0
        p1 = Polynomial(reverse(estimator1))
        ŷ1 = evalpoly.(u0, p1)
        rmss[k, 1] = sqrt(sum((ŷ1 - y0) .^ 2) / N)

        # numerically stable estimation
        estimator2 = K \ y0
        p2 = Polynomial(reverse(estimator2))
        ŷ2 = evalpoly.(u0, p2)
        rmss[k, 2] = sqrt(sum((y0 - ŷ2) .^ 2) / N)
    end
end

# ╔═╡ 64168108-a6b2-4930-a659-f7738ac45b6b
md"
!!! warning
	what is `KAll`? and why is the RMS figure different from the book?

	What is the purpose of this exercise in the first place
"

# ╔═╡ 02447fa2-5ad9-47a2-93b4-ed0650385a7c
plot(
    conds;
    yaxis = :log,
    ylims = (10e0, 10e20),
    xlabel = "Order",
    labels = reshape(["unstable", "stable"], 1, :),
    ylabel = "Condition number",
    legend = :topleft,
)

# ╔═╡ 7f6f06fa-a857-4554-ba14-6880c52e0b5d
plot(
    rmss;
    yaxis = :log,
    ylims = (10e-10, 10e10),
    xlabel = "Order",
    labels = reshape(["unstable", "stable"], 1, :),
    ylabel = "RMS error",
)

# ╔═╡ Cell order:
# ╠═6c9c5140-969d-11ee-30b0-5d1f501a232b
# ╠═aa954823-05a5-4343-a78d-10f034dddf56
# ╟─8a085c11-88d5-43e1-807a-865a92315440
# ╠═92c8f99e-b3ea-430f-9f60-1ca9f0d4a6a6
# ╟─64168108-a6b2-4930-a659-f7738ac45b6b
# ╟─02447fa2-5ad9-47a2-93b4-ed0650385a7c
# ╟─7f6f06fa-a857-4554-ba14-6880c52e0b5d

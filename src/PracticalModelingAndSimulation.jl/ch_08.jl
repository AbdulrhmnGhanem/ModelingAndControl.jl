### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 04073181-12ac-4d98-a106-7e591c933859
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end


# ╔═╡ 9aea4776-93f7-4873-ae14-179853c41921
using Plots, DifferentialEquations

# ╔═╡ 609d9590-1144-11ee-34dd-6da4e1760243
md"# Boundary Value Problems"

# ╔═╡ 55b6ee4a-fae2-47b7-91f6-4d76405c337d
md"
## Example 1
$ÿ + 2y = e^t$
$y = [y_0,\ y_{end}]$
$y_0 = y(0) = 0, y_{end} = y(π) = 1$
"

# ╔═╡ acc25675-167f-4a6b-8256-5c4dda1d5bc3
begin
    function ode_1(du, u, p, t)
        u₁, u₂ = u
        du[1] = u₂
        du[2] = exp(t) - 2u₁
    end

    function bc_1(residual, u, p, t)
        residual[1] = u[1][1]
        residual[2] = u[end][1] - 1
    end

    u0_1 = [0, pi]
    tspan_1 = (0, pi / 1)
    bvp_1 = BVProblem(ode_1, bc_1, u0_1, tspan_1)
    sol_1 = solve(bvp_1, GeneralMIRK4(), dt = 0.05)
end;

# ╔═╡ 049e5389-7705-4802-ab19-88ac94c12241
plot(sol_1; idxs = [1], label = nothing, ylabel = "u")

# ╔═╡ 64ba4356-b708-418d-8c5b-ea946ae79084
md"
## Example 2
$t^2ÿ - 5tẏ + 8y = 0$
$y(1) = 0, \ y(2) = 24$
"

# ╔═╡ a1581bc1-348b-441c-8535-a0918892989c
begin
    function ode_2(du, u, p, t)
        u₁, u₂ = u
        du[1] = u₂
        du[2] = (5t * u₂ - 8u₁) / t^2
    end

    function bc_2(residual, u, p, t)
        residual[1] = u[1][1]
        residual[2] = u[end][1] - 24
    end

    tspan_2 = (1, 2)
    u0_2 = [0.0, 24.0]
    bvp_2 = BVProblem(ode_2, bc_2, u0_2, tspan_2)
    sol_2 = solve(bvp_2, GeneralMIRK4(), dt = 0.05)
end;

# ╔═╡ 760b3e04-9bec-464f-adb0-c56229becec6
plot(sol_2; idxs = [1], label = nothing, ylabel = "u")

# ╔═╡ 0422a668-c4c3-49d4-b15f-1e3fcdd3fb9c
md"## Example 3
$ÿ + 2y = 0$
$BCs = \begin{cases}
y(0) = 1,\ ẏ(π) = 0 → \text{(Robin or mixed
boundary conditions)}\\
y(0) = 1,\ y(π) = 0 → \text{(Dirichlet
boundary conditions)}
\end{cases}$
"

# ╔═╡ e0ad6e37-8485-44bc-a4f6-29f562db11bf
begin
    function ode_3(du, u, p, t)
        u₁, u₂ = u
        du[1] = u₂
        du[2] = -2u₁
    end

    tspan_3 = (0, π / 1)
    u0_3 = [1.0, 0.0]
end;

# ╔═╡ 405fa717-ba76-4ff9-87f9-af01f1017dfb
md"### Robin conditions"

# ╔═╡ 7c0f6490-e043-4c01-8f96-9a89a26c5c35
begin
    function bc_3_robin(residual, u, p, t)
        residual[1] = u[1][1] - 1
        residual[2] = u[end][2]
    end

    bvp_3_robin = BVProblem(ode_3, bc_3_robin, u0_3, tspan_3)
    sol_3_robin = solve(bvp_3_robin, GeneralMIRK4(), dt = 0.05)
end;

# ╔═╡ 80c5a736-4b1a-4799-b910-c223bf38859a
md"### Dirichlet conditions"

# ╔═╡ eb3dad6c-cdae-4616-bb97-a6e94a4c4fa9
begin
    function bc_3_dirchlet(residual, u, p, t)
        residual[1] = u[1][1] - 1
        residual[2] = u[end][1]
    end

    bvp_3_dirichlet = BVProblem(ode_3, bc_3_dirchlet, u0_3, tspan_3)
    sol_3_dirichlet = solve(bvp_3_dirichlet, GeneralMIRK4(), dt = 0.05)
end;

# ╔═╡ 33ac91e8-43ec-4587-ad70-16158ce5c3ca
begin
    plot(sol_3_dirichlet; idxs = [1], label = "Dirchlet", ylabel = "u")

    plot!(sol_3_robin; idxs = [1], label = "robin")
end

# ╔═╡ bb94f9b7-4738-4d3d-b50e-458ce198d33c
md"## Example 4

$sinh(t) = - ü + ω^2u$
$\omega \in \{0, 3, 7, 13\}$
$u(0) = u(1) = 0$
"

# ╔═╡ cf74bf00-0f77-4276-84a7-39f9ebc7a683
begin
    function ode_4(du, u, ω, t)
        u₁, u₂ = u
        du[1] = u₂
        du[2] = ω^2 * u₁ - sinh(t)
    end

    function bc_4(residual, u, p, t)
        residual[1] = u[1][1]
        residual[2] = u[end][1]
    end

    tspan_4 = (0.0, 1)
    u0_4 = [0.0, 0.0]
    ωs = (0, 3, 7, 13)
    sols_4 = Vector{ODESolution}(undef, length(ωs))
    for (i, ω) in enumerate(ωs)
        bvp_4 = BVProblem(ode_4, bc_4, u0_4, tspan_4, ω)
        sol_4 = solve(bvp_4, GeneralMIRK4(), dt = 0.01)
        sols_4[i] = sol_4
    end
end;

# ╔═╡ ef9abe71-fefb-4b1a-8798-f94980c49dbb
begin
    plot(sols_4[1]; idxs = [1], label = "ω = 0", ylabel = "u")
    plot!(sols_4[2]; idxs = [1], label = "ω = $(ωs[2])")
    plot!(sols_4[3]; idxs = [1], label = "ω = $(ωs[3])")
    plot!(sols_4[4]; idxs = [1], label = "ω = $(ωs[4])")
end

# ╔═╡ 6998636f-6c79-4da4-975f-9042e2c7f597
md"## Example 5
$ÿ - (\alpha + \beta)ẏ + \alpha \beta y = 0$
$y(0) = y(1) = 1$
$\alpha = 100, \ \beta = -100$
"

# ╔═╡ dab27397-76c8-4eee-bab8-e024162494a2
begin
    function ode_5(du, u, p, t)
        u₁, u₂ = u
        α, β = p
        du[1] = u₂
        du[2] = (α + β) * u₂ - α * β * u₁
    end

    function bc_5(residual, u, p, t)
        residual[1] = u[1][1] - 1
        residual[2] = u[end][1] - 1
    end

    tspan_5 = (0.0, 1.0)
    u0_5 = [1.0, 1.0]
    p_5 = [100.0, -100.0]

    bvp_5 = BVProblem(ode_5, bc_5, u0_5, tspan_5, p_5)
    sol_5 = solve(bvp_5, GeneralMIRK4(), dt = 0.005)
end;

# ╔═╡ 5a2c7d8a-d672-483e-b794-bb1e165fa13a
plot(sol_5; idxs = [1], label = nothing, ylabel = "u")

# ╔═╡ Cell order:
# ╟─609d9590-1144-11ee-34dd-6da4e1760243
# ╠═04073181-12ac-4d98-a106-7e591c933859
# ╠═9aea4776-93f7-4873-ae14-179853c41921
# ╟─55b6ee4a-fae2-47b7-91f6-4d76405c337d
# ╠═acc25675-167f-4a6b-8256-5c4dda1d5bc3
# ╟─049e5389-7705-4802-ab19-88ac94c12241
# ╟─64ba4356-b708-418d-8c5b-ea946ae79084
# ╠═a1581bc1-348b-441c-8535-a0918892989c
# ╟─760b3e04-9bec-464f-adb0-c56229becec6
# ╟─0422a668-c4c3-49d4-b15f-1e3fcdd3fb9c
# ╠═e0ad6e37-8485-44bc-a4f6-29f562db11bf
# ╟─405fa717-ba76-4ff9-87f9-af01f1017dfb
# ╠═7c0f6490-e043-4c01-8f96-9a89a26c5c35
# ╟─80c5a736-4b1a-4799-b910-c223bf38859a
# ╠═eb3dad6c-cdae-4616-bb97-a6e94a4c4fa9
# ╟─33ac91e8-43ec-4587-ad70-16158ce5c3ca
# ╟─bb94f9b7-4738-4d3d-b50e-458ce198d33c
# ╠═cf74bf00-0f77-4276-84a7-39f9ebc7a683
# ╟─ef9abe71-fefb-4b1a-8798-f94980c49dbb
# ╟─6998636f-6c79-4da4-975f-9042e2c7f597
# ╠═dab27397-76c8-4eee-bab8-e024162494a2
# ╟─5a2c7d8a-d672-483e-b794-bb1e165fa13a

### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ 8921706e-dd46-11ee-06f4-8de991020711
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
    data_path = joinpath("..", "..", "books", "databook", "DATA")
end;

# ╔═╡ 513a0cbb-6af4-47a7-819f-52f74ded57b3
# ╠═╡ show_logs = false
using LsqFit, Plots, Interpolations, Optim, PlutoUI, StatsBase

# ╔═╡ 5504ba02-95f9-45a0-896d-edf3a1e11585
md"# Chapter 4 - Regression and Model Selection"

# ╔═╡ bbf81c22-a7c4-4a16-ad4f-3563e72abad8
md"## Exercise 4.1

Derive in closed form the 3 × 3 matrix which results from a least-squares
regression to a parabolic fit $f(x) = Ax^2 + Bx + C$.
"

# ╔═╡ 5c834ba4-b3de-44d5-838a-eba0c5c7e6da


# ╔═╡ b01ecc30-b29c-40c8-a59c-d2159440f651
md"## Exercise 4.2

Consider the following temperature data taken over a 24-hour (military time)
cycle

![](https://github.com/AbdulrhmnGhanem/ModelingAndControl.jl/blob/main/src/databook.jl/4.2.png?raw=true)

Fit the data with the parabolic fit

$f(x) = Ax^2 + Bx + C$

and calculate the $E_2$ error. Use both a linear interpolation and a spline to generate an interpolated approximation to the data for x = 1:0.01:24.
Develop a least-squares algorithm and calculate $E_2$ for

$y = Acos(Bx) + C$

Evaluate the resulting fit as a function of the initial guess for the values of A, B, and C.
"

# ╔═╡ c9127143-25df-4d93-aa6c-fd57fd3b029d
begin
	E₂(y, ŷ) = round(sqrt(sum(abs2.(y - ŷ)) / 2); digits=3)
	function solve_two()
		x = 1:24
		y = [75, 77, 76, 73, 69, 68, 63, 59, 57, 55, 54, 52, 50, 50, 49, 49, 49, 50, 54, 56, 59, 63, 67, 72]
		
		model(x, p) = p[1] * x.^2 .+ p[2] * x .+ p[3]
		p0 = [0.0, 50.0, 25.0]
		fit = curve_fit(model, x, y, p0)
		A, B, C = fit.param
		f(x) = A * x^2 + B * x + C
		f, x, y, E₂(f.(x), y)
	end
	
	function plot_two(sol)
		f, x, y, e₂ = sol
		
		scatter(x, y;
			label="y",
		)
		plot!(f, 1, 26;
			label="fitted curve",
		)
		plot!(x, y;
			label="linear interpolation",
			ls=:dash
		)
		plot!(x, interpolate(y, BSpline(Quadratic(Reflect(OnGrid()))));
			label="quadratic interpolation",
			ls=:dashdot,
			leg=:top,
			
		)
	end
end

# ╔═╡ ab98068e-d2d0-40d4-a763-da86a6ffec96
begin
	sol_two = solve_two()
	md"E₂ = $(sol_two[end])"
end

# ╔═╡ 723cdee8-b566-42f0-9107-5af3868396d7
plot_two(sol_two)

# ╔═╡ f1a6de55-221e-413e-bd1d-2cf4715e3b67
md"## Exercise 4.3

For the temperature data of the previous example, consider a polynomial fit
of the form

$f(x)  = \sum^{10}_{k=0} a_k x^k,$

where the loadings $α_k$ are to be determined by four regression techniques: least-squares, LASSO, ridge, and elastic net. Compare the models for each against each other. Randomly pick any time point and corrupt the temperature measurement at that location. For instance, the temperature reading at that location could be zero. Investigate the resulting model and $E_2$ error for the four regression techniques considered. Identify the models that are robust to such an outlier and those that are not. Explicitly calculate the variance of the loading coefficients $α_k$ for each method for a number of random trials with one or more corrupt data points.
"

# ╔═╡ 5cb3fe75-fb66-4862-a769-c4a067933221
begin
	function lasso_objective(p, model, x, y, λ)
	    ŷ = model(x, p)
	    mse_loss = sum((ŷ .- y).^2) / length(y)
	    l1_penalty = λ * sum(abs, p)
	    return mse_loss + l1_penalty
	end

	function ridge_objective(p, model, x, y, λ)
	    ŷ = model(x, p)
	    mse_loss = sum((ŷ .- y).^2) / length(y)
	    l2_penalty = λ * sqrt(sum(p .^ 2))
	    return mse_loss + l2_penalty
	end

	function elastic_objective(p, model, x, y, λ)
		ŷ = model(x, p)
	    mse_loss = sum((ŷ .- y).^2) / length(y)
	    l2_penalty = λ * sqrt(sum(p .^ 2))
	    l1_penalty = λ * sum(abs, p)
	    return mse_loss + l2_penalty + l1_penalty
	end
	
	function solve_three(enable_disturbance_mode)
		x = 1.0:24.0
		y = [75, 77, 76, 73, 69, 68, 63, 59, 57, 55, 54, 52, 50, 50, 49, 49, 49, 50, 54, 56, 59, 63, 67, 72]

		# simulate distrubance by setting one on the readings to 30
		if enable_disturbance_mode
			y[Int(rand(x))] = 30
		end
		
		model(x, p) = sum([p[i] * x .^ (i-1) for i in 1:length(p)])
		p0 = zeros(10)
		lsq_param = curve_fit(model, x, y, p0).param
		lsq_model = model(x, lsq_param)
	
		λ = 0.2
		lasso_param = optimize(
			p -> lasso_objective(p, model, x, y, λ), 
			p0, 
			NelderMead(), 
		   Optim.Options(iterations = 5000),
		) |> Optim.minimizer
		lasso_model = model(x, lasso_param)

		ridge_param = optimize(
			p -> ridge_objective(p, model, x, y, λ), 
			p0, 
			NelderMead(), 
		   Optim.Options(iterations = 5000),
		) |> Optim.minimizer
		ridge_model = model(x, ridge_param)

		elastic_param = optimize(
			p -> elastic_objective(p, model, x, y, λ), 
			p0, 
			NelderMead(), 
		   Optim.Options(iterations = 5000),
		) |> Optim.minimizer
		elastic_model = model(x, elastic_param)

		y, lsq_model, lasso_model, ridge_model, elastic_model
	end

	function variance()
		n = 1000
		lsq_var = zeros(n)
		lasso_var = zeros(n)
		ridge_var = zeros(n)
		elastic_var = zeros(n)
		
		Threads.@threads for i in 1:n
			y, lsq, lasso, ridge, elastic = solve_three(true)
			lsq_var[i] = var(lsq)
			lasso_var[i] = var(lasso)
			ridge_var[i] = var(ridge)
			elastic_var[i] = var(elastic)
		end

		return lsq_var, lasso_var, ridge_var, elastic_var
	end

	function plot_variance(variance)
		lsq_var, lasso_var, ridge_var, elastic_var = variance
		p1 = histogram(lsq_var; title="OLS")
		p2 = histogram(lasso_var; title="lasso")
		p3 = histogram(ridge_var; title="ridge")
		p4 = histogram(elastic_var; title="elastic net")

		plot(p1, p2, p3, p4; layout=(2,2), size=(500,500), legend=false)
	end
	
	function plot_three(sol)
		y, lsq_model, lasso_model, ridge_model, elastic_model = sol
		scatter(y;
			label="y",
			ylims=(20, 80),
		)
		plot!(lsq_model,
			label="OLS",
		)
		plot!(lasso_model,
			label="LASSO",
		)
		plot!(ridge_model,
			label="ridge",
		)
		plot!(elastic_model;
			label="elastic net",
			ls=:dash,
		)
	end
end

# ╔═╡ 9df6a985-45da-4097-911b-1796e81bb1ec
@bind enable_disturbance_mode CheckBox(default=false)

# ╔═╡ e88344c2-ebe7-46e1-8052-dac7eacca08d
@bind distrub Button("Disturb")

# ╔═╡ 09cfbbe8-0254-4a3a-9638-8d0b0e6e6832
begin
	distrub
	sol = solve_three(enable_disturbance_mode)
	d_y, d_lsq_model, d_lasso_model, d_ridge_model, d_elastic_model = sol
	plot_three(sol)
end

# ╔═╡ 5371855f-400a-453f-a70f-715268a68e23
begin
	y, lsq_model, lasso_model, ridge_model, elastic_model = solve_three(false)
md"
- E₂ (OLS) = $(E₂(y, lsq_model))
- E₂ (Disturbed OLS) = $(E₂(y, d_lsq_model))

- E₂ (lasso) = $(E₂(y, lasso_model))
- E₂ (Disturbed lasso) = $(E₂(y, d_lasso_model))

- E₂ (ridge) = $(E₂(y, ridge_model))
- E₂ (Disturbed ridge) = $(E₂(y, d_ridge_model))

- E₂ (elsatic net) = $(E₂(y, elastic_model))
- E₂ (Disturbed elastic) = $(E₂(y, d_elastic_model))
"
end

# ╔═╡ aa68a8ef-a08d-4f2a-b780-faa6743be03f
plot_variance(variance())

# ╔═╡ 4934d12b-f6d1-45a8-8b0a-fc5bfc8162a9
md"## Exercise 4.4
Download the MNIST data set (both training and test sets and labels) from
[http://yann.lecun.com/exdb/mnist](http://yann.lecun.com/exdb/mnist).

The labels will tell you which digit it is: 1, 2, 3, 4, 5, 6, 7, 8, 9, 0. Let each output be denoted by the vector $y_j$.

Now let **B** be the set of output vectors

$B = [y_1 \ y_2 \ y_3 \ ... \ y_n]$

and let the matrix A be the corresponding reshaped (vectorized) MNIST images

$A = [x_1 \ x_2 \ x_3 \ ... \ x_n]$

Thus each vector $x_j \in R^{n^2}$ is a vector reshaped from the $n × n$ image.
Using various **AX = B** solvers, determine a mapping from the image space to the label space.

By promoting sparsity, determine and rank which pixels in the MNIST set are most
informative for correctly labeling the digits. (You will have to come up with your own
heuristics or empirical rules for this. Be sure to visualize the results from **X**.) Apply your most important pixels to the test data set to see how accurate you are with as few pixels as possible. Redo the analysis with each digit individually to find the most important pixels for each digit. Think about the interpretation of what you are doing with this **AX = B** problem.
"

# ╔═╡ 57172e34-146f-4930-ba78-9285f5f5f7f7


# ╔═╡ Cell order:
# ╠═8921706e-dd46-11ee-06f4-8de991020711
# ╠═513a0cbb-6af4-47a7-819f-52f74ded57b3
# ╟─5504ba02-95f9-45a0-896d-edf3a1e11585
# ╟─bbf81c22-a7c4-4a16-ad4f-3563e72abad8
# ╠═5c834ba4-b3de-44d5-838a-eba0c5c7e6da
# ╟─b01ecc30-b29c-40c8-a59c-d2159440f651
# ╠═c9127143-25df-4d93-aa6c-fd57fd3b029d
# ╟─ab98068e-d2d0-40d4-a763-da86a6ffec96
# ╟─723cdee8-b566-42f0-9107-5af3868396d7
# ╟─f1a6de55-221e-413e-bd1d-2cf4715e3b67
# ╠═5cb3fe75-fb66-4862-a769-c4a067933221
# ╟─9df6a985-45da-4097-911b-1796e81bb1ec
# ╟─e88344c2-ebe7-46e1-8052-dac7eacca08d
# ╟─09cfbbe8-0254-4a3a-9638-8d0b0e6e6832
# ╟─5371855f-400a-453f-a70f-715268a68e23
# ╟─aa68a8ef-a08d-4f2a-b780-faa6743be03f
# ╟─4934d12b-f6d1-45a8-8b0a-fc5bfc8162a9
# ╠═57172e34-146f-4930-ba78-9285f5f5f7f7

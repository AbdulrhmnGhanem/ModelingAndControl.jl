### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 8921706e-dd46-11ee-06f4-8de991020711
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
    data_path = joinpath("..", "..", "books", "databook", "DATA")
end;

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


# ╔═╡ f1a6de55-221e-413e-bd1d-2cf4715e3b67
md"## Exercise 4.3

For the temperature data of the previous example, consider a polynomial fit
of the form

$f(x)  = \sum^{10}_{k=0} a_k x^k,$

where the loadings αk are to be determined by four regression techniques: least-squares, LASSO, ridge, and elastic net. Compare the models for each against each other. Randomly pick any time point and corrupt the temperature measurement at that location. For instance, the temperature reading at that location could be zero. Investigate the resulting model and E2 error for the four regression techniques considered. Identify the models that are robust to such an outlier and those that are not. Explicitly calculate the variance of the loading coefficients αk for each method for a number of random trials with one or more corrupt data points.
"

# ╔═╡ 5cb3fe75-fb66-4862-a769-c4a067933221


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
# ╟─5504ba02-95f9-45a0-896d-edf3a1e11585
# ╟─bbf81c22-a7c4-4a16-ad4f-3563e72abad8
# ╠═5c834ba4-b3de-44d5-838a-eba0c5c7e6da
# ╟─b01ecc30-b29c-40c8-a59c-d2159440f651
# ╠═c9127143-25df-4d93-aa6c-fd57fd3b029d
# ╟─f1a6de55-221e-413e-bd1d-2cf4715e3b67
# ╠═5cb3fe75-fb66-4862-a769-c4a067933221
# ╟─4934d12b-f6d1-45a8-8b0a-fc5bfc8162a9
# ╠═57172e34-146f-4930-ba78-9285f5f5f7f7

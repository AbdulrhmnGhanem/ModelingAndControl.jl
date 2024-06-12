### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 890941ea-28a0-11ef-0ae3-73987e387256
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
    data_path = joinpath("..", "..", "books", "databook", "DATA")
end;

# ╔═╡ 747cbaa4-f22e-4190-95e3-a180ab39ae2e
using Flux, CUDA

# ╔═╡ 9138963c-8b3e-4ded-966b-1d67bb84dad7
md"# Chapter 6 - Neural Networks and Deep Learning"

# ╔═╡ ebcc9661-0840-4935-b909-4011fe31ceea
md"## Exercise 6.1

Download the code base for solving (i) a reaction–diffusion system of equa-
tions, and (ii) the Kuramoto–Sivashinsky (KS) equation.
1. Train a NN that can advance the solution from t to t + Δt for the KS equation.
2. Compare your evolution trajectories for your NN against using the ODE time-stepper provided with different initial conditions.
3. For the reaction–diffusion system, first project to a low-dimensional subspace via the SVD and see how forecasting works in the low-rank variables.

For the Lorenz equations, consider the following.

4. Train an NN to advance the solution from t to t + Δt for ρ = 10, 28, and 40. Now see how well your NN works for future state prediction for ρ = 17 and ρ = 35.

5. See if you can train your NN to identify (for ρ = 28) when a transition from one lobe to another is imminent. Determine how far in advance you can make this prediction. (Note: You will have to label the transitions in a test set in order to do this task.)
"

# ╔═╡ d7bf087d-e488-4aeb-bc3a-00fbd7bd2c7c


# ╔═╡ cad54b3c-5c5e-41a7-b49c-caedc53e83fa
md"## Exercise 6.2

Consider time-series data acquired from power grid loads, specifically: T. V.
Jensen and P. Pinson. Re-Europe, a large-scale data set for modeling a highly renewable European electricity system. _Scientific Data_, 4:170175, 2017. Compare the forecasting capabilities of the following neural networks on the power grid data: (i) a feedforward neural network; (ii) an LSTM; (iii) an RNN; and (iv) an echo state network. Consider the performance of each under cross-validation for forecasting ranges of t into the future and Nt into the future (where N ≫ 1).
"

# ╔═╡ 45c64e3c-f803-4f7b-b423-4929f59067e4


# ╔═╡ 47c34e84-0f3f-4aad-bdff-cc7c6e2815fd
md" ## Exercise 6.3

Download the flow around the cylinder data. Using the first _P_% of the tem-
poral snapshots, forecast the remaining (100 − _P_)% future state data. Do this by training a neural network on the high-dimensional data and using: (i) a feedforward neural network; (ii) an LSTM; (iii) an RNN; and (iv) an echo state network. Determine the performance of the algorithms as a function of decreasing data _P_.

Redo the forecasting calculations by training a model in the reduced subspace U from the singular value decomposition. Evaluate the forecasting performance as a function of the percentage of the training data P and the rank of the reduced space r.
"

# ╔═╡ 29740cbc-4360-4f92-9fb0-0dcee0388059


# ╔═╡ 8612a792-127a-43a2-9ad1-1e2f589603a0
md"## Exercise 6.4

Generate simulation data for the Kuramoto–Sivashinsky (KS) equation in
three distinct parameter regimes where non-trivial spatio-temporal dynamics occurs. Using a convolutional neural network, map the high-dimensional snapshots of the system to a classification of the system into one of the three parameter regimes. Evaluate the performance of the classification scheme on test data as a function of different convolutional window sizes and stride lengths. For the best performance, what is the convolutional window size and what spatial length scale is extracted to make the classification decision?"

# ╔═╡ a5400a68-24e3-4210-9093-6bfca7b1c114


# ╔═╡ Cell order:
# ╠═890941ea-28a0-11ef-0ae3-73987e387256
# ╠═747cbaa4-f22e-4190-95e3-a180ab39ae2e
# ╟─9138963c-8b3e-4ded-966b-1d67bb84dad7
# ╟─ebcc9661-0840-4935-b909-4011fe31ceea
# ╠═d7bf087d-e488-4aeb-bc3a-00fbd7bd2c7c
# ╟─cad54b3c-5c5e-41a7-b49c-caedc53e83fa
# ╠═45c64e3c-f803-4f7b-b423-4929f59067e4
# ╟─47c34e84-0f3f-4aad-bdff-cc7c6e2815fd
# ╠═29740cbc-4360-4f92-9fb0-0dcee0388059
# ╟─8612a792-127a-43a2-9ad1-1e2f589603a0
# ╠═a5400a68-24e3-4210-9093-6bfca7b1c114

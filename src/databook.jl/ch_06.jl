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
using Flux, CUDA, FFTW, DifferentialEquations, Plots, JLD2, StatsBase, ProgressLogging

# ╔═╡ 9138963c-8b3e-4ded-966b-1d67bb84dad7
md"# Chapter 6 - Neural Networks and Deep Learning"

# ╔═╡ ebcc9661-0840-4935-b909-4011fe31ceea
md"## Exercise 6.1

Download the code base for solving (i) a reaction–diffusion system of equations, and (ii) the Kuramoto–Sivashinsky (KS) equation.
1. Train a NN that can advance the solution from t to t + Δt for the KS equation.
2. Compare your evolution trajectories for your NN against using the ODE time-stepper provided with different initial conditions.
3. For the reaction–diffusion system, first project to a low-dimensional subspace via the SVD and see how forecasting works in the low-rank variables.

For the Lorenz equations, consider the following.

4. Train an NN to advance the solution from t to t + Δt for ρ = 10, 28, and 40. Now see how well your NN works for future state prediction for ρ = 17 and ρ = 35.

5. See if you can train your NN to identify (for ρ = 28) when a transition from one lobe to another is imminent. Determine how far in advance you can make this prediction. (Note: You will have to label the transitions in a test set in order to do this task.)
"

# ╔═╡ d783d3b1-d528-4e81-9085-65ca7678d67e
md"### The ODE time stepper solution"

# ╔═╡ 65d2550e-5d4e-46ae-a83f-1294f4d82683
function ks(u₀s)
    tspan = (0, 10)
	L = 100.0
    N = 1000
    dx = L / N
	κ = (2π / L) * (-N/2:N/2-1)	
    κ = fftshift(κ, 1) 

    function rhs(dudt, u, p, t)
        û = fft(u)
        dû = im * κ .* û
        ddû = -(κ .^ 2) .* û
        ddddû = (κ .^ 4) .* û
        du = ifft(dû)
        ddu = ifft(ddû)
        ddddu = ifft(ddddû)
        dudt .= real(0.5 .* du .^ 2 - ddu - ddddu)
    end
    t = tspan[1]:tspan[2]/200:tspan[2]
    prob = ODEProblem(rhs, u₀s[1, :], tspan)

	sols = zeros(Float32, size(u₀s, 2), size(u₀s, 1))
	for (i, u₀) in collect(enumerate(eachrow(u₀s)))
		prob = remake(prob; u0=u₀s[1, :])
    	u = solve(prob, AutoTsit5(Rosenbrock23(autodiff = false)))
		@inbounds sols[:, i] = u[:, end]
	end
    sols, t
end

# ╔═╡ 2eefce96-2d4e-4a75-b049-f5dc533c8715
function ks_initial_conditions(num_samples)
	L = 100
	N = 1000
	dx = L / N
	domain = range(-L / 2, stop = L / 2 - dx, length = N)

    initial_conditions = zeros(Float32, num_samples, N)
    for i in 1:num_samples        
        # Apply a random shift to the domain
        shift = rand() * L - L / 2
        shifted_domain = domain .+ shift

        case = rand(1:8)
        if case == 1
            # Sine wave with random frequency and amplitude
            freq = rand(1:5)
            amp = rand()
            u₀ = amp * sin.(2 * π * freq * collect(shifted_domain) / L)
        
        elseif case == 2
            # Gaussian pulse with random width and height
            width = L / (5 + rand(1:5))
            height = rand()
            u₀ = height * exp.(-collect(shifted_domain).^2 / (2 * width^2))
        
        elseif case == 3
            # Random initial condition
            u₀ = rand(N)
        
        elseif case == 4
            # Step function with random step location
            step_loc = rand(-L/2:L/2)
            u₀ = [x < step_loc ? -1.0 : 1.0 for x in collect(shifted_domain)]
        
        elseif case == 5
            # Combination of multiple sine waves
            freq1 = rand(1:5)
            freq2 = rand(6:10)
            amp1 = rand()
            amp2 = rand()
            u₀ = amp1 * sin.(2 * π * freq1 * collect(shifted_domain) / L) + amp2 * sin.(4 * π * freq2 * collect(shifted_domain) / L)
        
        elseif case == 6
            # Parabolic initial condition
            u₀ = 1 .- (collect(shifted_domain) ./ (L / 2)).^2
        
        elseif case == 7
            # Exponential initial condition
            decay_rate = L / (10 + rand(1:10))
            u₀ = exp.(-abs.(collect(shifted_domain)) / decay_rate)

		elseif case == 8
			u₀ = sech.(shifted_domain)
        end
        
        initial_conditions[i, :] = u₀
    end
    
    return initial_conditions
end

# ╔═╡ 21dceb6a-02c8-4df8-a7af-261a490a417d
function generate_training_data(num_samples, save=true)
	f = "ks.jld2" 
	if save && isfile(f)
		fi = jldopen(f)
		if fi["num_samples"] == num_samples
			return fi["sols"], fi["u₀s"], fi["time"]
		end
	close(fi)
	end

	u₀s = ks_initial_conditions(num_samples)
	sols, time = ks(u₀s)

	save && jldsave(f; sols, time, num_samples, u₀s)
	sols, u₀s, time
end

# ╔═╡ a2f3f620-d100-41e0-b82c-60b673baf9b5
function feedforward_ks()
	batch_size = 16
	sols, u₀s, t = generate_training_data(30batch_size)
	input = @view sols[1:end-1, :]
	output = @view sols[2:end, :]

	s = size(input, 1)
	model = Chain(
		Dense(s => 2s, sigmoid_fast),
		Dense(2s => 2s, relu),
		Dense(2s => s)
	) |> gpu

		
	optim = Flux.setup(Flux.Adam(0.01), model)
	folds = Flux.kfolds(Flux.shuffleobs((input, output)), k = 4)
	
	all_losses = []

	i = 1
    for ((train_input, train_output), (val_input, val_output)) in folds
		train_loader = Flux.DataLoader((train_input, train_output) |> gpu, batchsize=batch_size, shuffle=true)
        
		val_loader = Flux.DataLoader((val_input, val_output) |> gpu, batchsize=batch_size)

        losses = []
        @progress for epoch in 1:1_000
            for (x, y) in train_loader
                loss, grads = Flux.withgradient(model) do m
                    # Evaluate model and loss inside gradient context:
                    y_hat = m(x)
                    Flux.mse(y_hat, y)
                end
                Flux.update!(optim, model, grads[1])
                push!(losses, loss)  # logging, outside gradient context
            end
        end

		push!(all_losses, losses)

        # Evaluate on validation fold in batches
        val_losses = []
        for (x, y) in val_loader
            val_loss = Flux.mse(model(x), y)
            push!(val_losses, val_loss)
        end
        avg_val_loss = mean(val_losses)
		@info "Validation loss for fold $i: $avg_val_loss\n"
		i += 1
    end

	all_losses, model
end

# ╔═╡ 90a2bfdd-c817-4927-9111-f30e1e3229c6
losses, model = feedforward_ks()

# ╔═╡ 1be89559-9e0e-4c40-9e12-2aece27d87ab
a = generate_training_data(2, false)

# ╔═╡ d0469011-615c-43b5-b4c4-60bb03b33173
plot(a[2][1,:])

# ╔═╡ 66ed9a87-5bd7-49f9-9d1c-1621725b4aab
plot(model(a[2][1, 1:end-1] |> gpu) |> cpu)

# ╔═╡ 2bcb4c1e-5510-4b0c-9311-0923203cf447
plot(ks(a[2][1,:]')[1])

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
# ╟─d783d3b1-d528-4e81-9085-65ca7678d67e
# ╠═65d2550e-5d4e-46ae-a83f-1294f4d82683
# ╟─2eefce96-2d4e-4a75-b049-f5dc533c8715
# ╠═21dceb6a-02c8-4df8-a7af-261a490a417d
# ╠═a2f3f620-d100-41e0-b82c-60b673baf9b5
# ╠═90a2bfdd-c817-4927-9111-f30e1e3229c6
# ╠═1be89559-9e0e-4c40-9e12-2aece27d87ab
# ╠═d0469011-615c-43b5-b4c4-60bb03b33173
# ╠═66ed9a87-5bd7-49f9-9d1c-1621725b4aab
# ╠═2bcb4c1e-5510-4b0c-9311-0923203cf447
# ╟─cad54b3c-5c5e-41a7-b49c-caedc53e83fa
# ╠═45c64e3c-f803-4f7b-b423-4929f59067e4
# ╟─47c34e84-0f3f-4aad-bdff-cc7c6e2815fd
# ╠═29740cbc-4360-4f92-9fb0-0dcee0388059
# ╟─8612a792-127a-43a2-9ad1-1e2f589603a0
# ╠═a5400a68-24e3-4210-9093-6bfca7b1c114

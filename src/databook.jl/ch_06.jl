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
using Flux, CUDA, FFTW, DifferentialEquations, Plots, JLD2

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

# ╔═╡ ea583da7-a613-4311-a7e9-96ca3654262f
length(0:10/400:10)

# ╔═╡ d783d3b1-d528-4e81-9085-65ca7678d67e
md"### The ODE time stepper solution"

# ╔═╡ 65d2550e-5d4e-46ae-a83f-1294f4d82683
function ks(u₀s, κ)
    tspan = (0, 10)

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

	sols = zeros(Float32, length(u₀s[1, :]), length(t), size(u₀s, 1))
	for (i, u₀) in collect(enumerate(eachrow(u₀s)))
		prob = remake(prob; u0=u₀)
    	u = solve(prob, AutoTsit5(Rosenbrock23(autodiff = false)), saveat = t)
		@inbounds sols[:, :, i] = u[1:end, 1:end]
	end
	
    sols, t
end

# ╔═╡ 2eefce96-2d4e-4a75-b049-f5dc533c8715
function ks_initial_conditions(num_samples)
	L = 100
	N = 1000
	dx = L / N
	domain = range(-L / 2, stop = L / 2 - dx, length = N)

    initial_conditions = zeros(num_samples, N)
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
function generate_training_data(num_samples)
	f = "ks.jld2" 
	if isfile(f)
		f = jldopen(f)
		if f["num_samples"] == num_samples
			return f["sols"], f["time"]
		end
	close(f)
	f = "ks.jld2"
	end
	
	L = 100.0
    N = 1000
    dx = L / N
	κ = (2π / L) * (-N/2:N/2-1)	
    κ = fftshift(κ, 1) 
    domain = range(-L / 2, stop = L / 2 - dx, length = N)
	u₀s = ks_initial_conditions(num_samples)
	sols, time = ks(u₀s, κ)

	jldsave(f; sols, time, num_samples)
	sols, time
end

# ╔═╡ 1733c9f7-f51c-46e7-8f47-8b38a570b794
sols, time = generate_training_data(300)

# ╔═╡ d7bf087d-e488-4aeb-bc3a-00fbd7bd2c7c
# ╠═╡ disabled = true
#=╠═╡
begin
	L = 100.0
    N = 1000
    dx = L / N
	κ = (2π / L) * (-N/2:N/2-1)	
    κ = fftshift(κ, 1)
    domain = range(-L / 2, stop = L / 2 - dx, length = N)
	u₀s = ks_initial_conditions(1)
	sols, time = ks(u₀s, κ)
	sol = sols[:, :, 1]
end
  ╠═╡ =#

# ╔═╡ feb5e7e1-112d-4805-a5ff-b2906e20c90e
 @gif for i = 1:length(eachcol(sol))
	plot(
		domain,
		sol[:, i];
		# ylims = (-0.5, 1.75),
		legend = false,
		title = "t = $(round(time[i] *1000; digits=2)) ms",
	)
end

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
# ╠═ea583da7-a613-4311-a7e9-96ca3654262f
# ╟─d783d3b1-d528-4e81-9085-65ca7678d67e
# ╠═65d2550e-5d4e-46ae-a83f-1294f4d82683
# ╠═2eefce96-2d4e-4a75-b049-f5dc533c8715
# ╠═21dceb6a-02c8-4df8-a7af-261a490a417d
# ╠═1733c9f7-f51c-46e7-8f47-8b38a570b794
# ╠═d7bf087d-e488-4aeb-bc3a-00fbd7bd2c7c
# ╠═feb5e7e1-112d-4805-a5ff-b2906e20c90e
# ╟─cad54b3c-5c5e-41a7-b49c-caedc53e83fa
# ╠═45c64e3c-f803-4f7b-b423-4929f59067e4
# ╟─47c34e84-0f3f-4aad-bdff-cc7c6e2815fd
# ╠═29740cbc-4360-4f92-9fb0-0dcee0388059
# ╟─8612a792-127a-43a2-9ad1-1e2f589603a0
# ╠═a5400a68-24e3-4210-9093-6bfca7b1c114

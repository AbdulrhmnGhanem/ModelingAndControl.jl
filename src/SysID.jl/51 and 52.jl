### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 4b5f1008-b961-11ee-028f-2d0956cc3f1f
# ╠═╡ show_logs = false
begin 
 # If you are running this notebook as a stannalone notebook disable this cell.
 import Pkg 
 Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ fee0201b-075a-4ec8-84e1-0c43be6393f1
using Plots, DSP, FFTW, StatsBase

# ╔═╡ 4e78b344-2745-47da-9f61-ffcaeae36ca3
begin
	function cheby1(n, r, wp)
	    h = digitalfilter(Lowpass(wp), Chebyshev1(n, r))
		h
	end
	
	function freqz(filt::Union{FilterCoefficients,Vector}, fₛ::Int, lower::Int, upper::Int, NPeriod::Int)
	    frequencies = range(lower, upper, length=NPeriod)
		if filt isa Vector
			@assert length(filt) == 2
			filt = PolynomialRatio(filt[1], filt[2])
		end
		H = freqresp(filt, frequencies * 2π / fₛ)
	    return frequencies, H
	end
	
	function random_phase_mutli_sine(excited_harm, N, total_size)
	    u = Vector{Float64}(undef, total_size)
	    S = zeros(ComplexF64, N)
	    
		for i in 1:total_size÷N
	        S[excited_harm .+ 1] .= exp.(im .* 2 * π * rand(size(excited_harm)))
	        r = 2 * real(ifft(S))
	        r /= std(r)  # rms value = 1
	        u[(i-1)*N + 1:i*N] .= r
	    end
	    u
	end
	md"### Matlab compatability utils"
end

# ╔═╡ ef83595b-ca54-4663-8e20-f6459535306d
md"## Exercise 51 and 52: Estimating the FRF and noise power spectrum using the local ploynomial method
!!! purpose
	- We have seen that the best selection for the exictation signal is an ergodic multisine, and wait for the transient to vanish. But it has two issues:
	1. Wating for the transient is wasted measurement time.
	2. In some cases we can't apply multisine excitation signal to the system.

	In this exercise the local polynomial method is used to estimate the system using random excitation with the a tradeoff between execution time and accuracy.
"

# ╔═╡ 706550a2-d960-4bdc-bff4-2a8477613245
begin
	function local_polyn_method(data_poly, method)
	    w = method.Bandwidth
	    U = data_poly.U
	    Y = data_poly.Y
	    N = length(U)
	
	    # Window lines
	    w1 = ceil(w / 2)  # calculation of the window width
	    w_lines = -w1:w1
	
	    G = zeros(ComplexF64, N)
	    varY = zeros(N)
	
	    for k = N:-1:1
	        select_lines = collect(k .+ w_lines)
	        k_freq = k  # frequency at which an estimate is made
	        f_fit = collect(w_lines)
	
	        # Edges setting
	        if select_lines[end] > N  # keep window within upper limits of the measurement band
	            f_shift = select_lines[end] - N
	            select_lines .= select_lines .- f_shift
	            f_fit .= w_lines .- f_shift
	        end
	
	        if select_lines[1] < 1  # keep window within lower limits of the measurement band
	            f_shift = abs(select_lines[1]) + 1
	            select_lines .= select_lines .+ f_shift
	            f_fit .= w_lines .+ f_shift
	        end
	
	        # The polynomial fit
	        U_B = U[select_lines .|> Int]  # select lines to be fitted
	        Y_B = Y[select_lines .|> Int]
	        f_fit = select_lines .- k  # select the corresponding frequencies
	
	        f_fit_scaled = f_fit / w1  # normalize the frequencies between [-1, 1]
	
	        G_fit, Y_var = local_lpm_local(Y_B, U_B, f_fit_scaled, method)
	
	        # Evaluate the results
	        G[k] = G_fit[1]
	        varY[k] = Y_var
	    end
	
	    return G, varY
	end
	
	function local_lpm_local(Y, U, f, method)
	    Y = Y
	    U = U
	    f = f
	
	    n_order = method.Order
	
	    K = zeros(ComplexF64, length(f), (n_order + 1) * (1 + method.Transient))
	
	    K[:, 1] .= U
	    counter = 1
	    for k in 1:n_order  # U .* f^k    input contributions
	        counter += 1
	        K[:, counter] .= U .* (f.^k)
	    end
	
	    if method.Transient
	        for k in 0:n_order  # 1 .* f^k    transient contributions
	            counter += 1
	            K[:, counter] .= f.^k
	        end
	    end  # end transient contribution
	
	    # solve equations
	    Theta = K \ Y
	
	    G = Theta[1]
	    p = Theta[1:n_order+1]
	
	    # estimate the variance
	    Y_error = Y - K * Theta  # calculate residual
	    q = length(Y) - length(Theta)  # number of (complex) degrees of freedom
	    Y_var = sum(abs2, Y_error) / q  # estimate complex variance
	
	    return G, Y_var
	end

	
	struct DataPoly
	    U::Vector{ComplexF64}
	    Y::Vector{ComplexF64}
		Freq::Vector{Float64}
	end

	struct MethodPoly
	   	Order::Int8
		Bandwidth::Int8
		Transient::Bool
	end
end;

# ╔═╡ f335a531-ea16-4152-80af-c01915760985
begin
	M = 16
	NPeriod = 128
	N = NPeriod * M
	fs = NPeriod
	Ampl = 1
	NTrans = 1024
		
	Lines = 1:NPeriod ÷ 2
	f = (Lines .- 1) / NPeriod * fs
	
	h = cheby1(2, 10, 0.2)
	
	u = randn(N + NTrans) * Ampl
	y = filt(h, u)
	deleteat!(y, 1:NTrans)
	deleteat!(u, 1:NTrans)
	
	y = reshape(y, NPeriod, M)
	u = reshape(u, NPeriod, M)
	wind = kron(hanning(NPeriod), ones(1, M))
	Y = fft(y .* wind, 1) / sqrt(NPeriod / 2)
	Y = Y[Lines, :]
	U = fft(u .* wind, 1) / sqrt(NPeriod / 2)
	U = U[Lines, :]

	UU = mean(abs2.(U[:, 1:M]), dims=2)
	YU = mean(Y[:, 1:M] .* conj.(U[:, 1:M]), dims=2)
	G = YU ./ UU
	
	# true FRF of the system at Hann frequencies
	_, G0Hann = freqz(h, fs, 0, 63, NPeriod)
	
	# Processing with the local polynomial method
	U = fft(u[:], 1) / sqrt(N)
	Y = fft(y[:], 1) / sqrt(N)
	# zero out the negative frequencies
	U[N ÷ 2 + 1:end] .= 0.0            
	Y[N ÷ 2 + 1:end] .= 0.0
	
	# # Prepare the data
	data = DataPoly(U, Y, (0:N ÷ 2) / N)
	method = MethodPoly(2, 6, true)
	GLocPol, varY = local_polyn_method(data, method)
	
	# # Create frequency vector for local polynomial method
	fLocPol = (0:N - 1) / N * fs
	
	# # True FRF of the system at local polynomial frequencies
	_, G0LocPol = freqz(h, fs, 0, 63, NPeriod)
end;

# ╔═╡ e8a353ea-26c2-4051-8cab-53f00bb86b7c
begin
	plot(f, G .|> abs .|> amp2db;
		xlabel="Frequency (Hz)",
		ylabel="Amplitude (dB)",
		ylims=(-150, 5),
		xlims=(0, 63),
		ls=:dash,
		label="G₀"
	)
	scatter!(f, amp2db.(abs.(G0Hann[1:2:end])); ms=2,label="G₀ Hann")
	plot!(fLocPol, amp2db.(abs.(GLocPol)); label="G local poly", ls=:dot)
end

# ╔═╡ 243a69b3-6e32-4ffd-967c-83fb93936d6e
scatter(fLocPol, [Y, varY / 2] .|> x -> x .|> abs .|> amp2db;
	xlabel="Frequency (Hz)",
	ylabel="Amplitude (dB)",
	ylims=(-300, 5),
	xlims=(0, 63),
	labels=reshape(["Y", "σ²"], 1, :),
	ms=2
)

# ╔═╡ Cell order:
# ╠═4b5f1008-b961-11ee-028f-2d0956cc3f1f
# ╠═fee0201b-075a-4ec8-84e1-0c43be6393f1
# ╠═4e78b344-2745-47da-9f61-ffcaeae36ca3
# ╟─ef83595b-ca54-4663-8e20-f6459535306d
# ╠═706550a2-d960-4bdc-bff4-2a8477613245
# ╠═f335a531-ea16-4152-80af-c01915760985
# ╟─e8a353ea-26c2-4051-8cab-53f00bb86b7c
# ╟─243a69b3-6e32-4ffd-967c-83fb93936d6e

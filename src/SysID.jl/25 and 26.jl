### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 89c57041-7bf5-4b60-8263-92177041ca22
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 62c7805d-f083-404d-8f4a-1a732d73920b
using DataStructures, Plots, FFTW

# ╔═╡ bcc0c1b6-196c-452c-9a3d-d73cbd6fc93b
md"## Exercise 25 and 26: Generation of maximum length binary sequence
!!! outcomes
	 -  An MLBS injects a maximum power for a given amplitude, its crest factor equals one.
	- The amplitude spectrum of the discrete time sequence is flat.
	- The amplitude spectrum of the ZOH sequence drops as a sine-function with its first zero crossing at the clock frequency
	- Within the frequency band of interest, the amplitude spectrum rolls off slowly **(but still finite)**.
"

# ╔═╡ 4399fb37-636a-4a01-b4af-69fc07769c1c
"""
Feedback configuration for generatinbg maximum length binary sequence
Keith Godfrey

Refernce: K. Godfrey, Perturbation Signals for System Identification. Prentice Hall, 1993.
"""
function mlbs_fedback(n)
    forks = if n == 1
        (1,)
    elseif n == 2
        (1, 2)
    elseif n == 3
        (2, 3)
    elseif n == 4
        (3, 4)
    elseif n == 5
        (3, 5)
    elseif n == 6
        (5, 6)
    elseif n == 7
        (4, 7)
    elseif n == 8
        (4, 5, 6, 8)
    elseif n == 9
        (5, 9)
    elseif n == 10
        (7, 10)
    elseif n == 11
        (9, 11)
    elseif n == 12
        (6, 8, 11, 12)
    elseif n == 13
        (9, 10, 12, 13)
    elseif n == 14
        (4, 8, 13, 14)
    elseif n == 15
        (14, 15)
    elseif n == 16
        (4, 13, 15, 16)
    elseif n == 17
        (14, 17)
    elseif n == 18
        (11, 12)
    elseif n == 19
        (14, 17, 18, 19)
    elseif n == 20
        (17, 20)
    elseif n == 21
        (19, 21)
    elseif n == 22
        (21, 22)
    elseif n == 23
        (18, 23)
    elseif n == 24
        (17, 22, 23, 24)
    elseif n == 25
        (22, 25)
    elseif n == 26
        (20, 24, 25, 26)
    elseif n == 27
        (22, 25, 26, 27)
    elseif n == 28
        (25, 28)
    elseif n == 29
        (27, 29)
    elseif n == 30
        (7, 28, 29, 30)
    elseif n == 31
        (28, 31)
    elseif n == 32
        (10, 30, 31, 32)
    elseif n == 33
        (20, 33)
    elseif n == 34
        (7, 32, 33, 34)
    elseif n == 35
        (2, 35)
    elseif n == 36
        (1, 2, 4, 5, 6, 36)
    elseif n == 37
        (1, 2, 3, 4, 5, 37)
    elseif n == 38
        (1, 5, 6, 38)
    elseif n == 39
        (4, 39)
    elseif n == 40
        (3, 4, 5, 40)
    elseif n == 41
        (3, 41)
    elseif n == 42
        (1, 2, 3, 4, 5, 42)
    elseif n == 43
        (3, 4, 6, 43)
    elseif n == 44
        (2, 5, 6, 44)
    elseif n == 45
        (1, 3, 4, 45)
    elseif n == 46
        (1, 2, 3, 5, 8, 46)
    elseif n == 47
        (5, 47)
    elseif n == 48
        (1, 2, 3, 5, 8, 46)
    elseif n == 49
        (4, 5, 6, 49)
    elseif n == 50
        (2, 3, 4, 50)
    elseif n == 51
        (1, 3, 6, 51)
    elseif n == 52
        (3, 52)
    elseif n == 53
        (1, 2, 6, 53)
    elseif n == 54
        (2, 3, 4, 5, 6, 54)
    elseif n == 55
        (1, 2, 6, 55)
    elseif n == 56
        (2, 4, 7, 56)
    elseif n == 57
        (2, 3, 5, 57)
    elseif n == 58
        (1, 5, 6, 56)
    elseif n == 59
        (1, 3, 4, 5, 6, 59)
    elseif n == 60
        (1, 60)
    elseif n == 61
        (1, 2, 5, 61)
    elseif n == 62
        (3, 5, 6, 62)
    else
        throw(error("n has to be between 1 and 62, got $n instead!"))
    end

    state -> reduce(⊻, [state[s] for s in forks])
end

# ╔═╡ c3269876-460c-469b-8cc7-bda1604657be
"""
Pseudo-random binary sequnce base on maximum length binary sequence algorithm.
Maximum allowed order `n=62`.
The signal switched between `dc` and `-dc`, defualt is binary levels: 1, 0.
"""
function mlbs(n, dc = nothing)
    maximum_length = begin
        state = CircularBuffer{Bool}(n)
        # intialize the buffer to all zeros expect for the last one
        append!(state, zeros(n))
        push!(state, 1)
    end

    maximum_length = (2^capacity(state)) - 1
    signal = Vector{Int8}(undef, maximum_length)
    feedback = mlbs_fedback(n)

    if isnothing(dc)
        dc = 1
        ndc = 0
    else
        ndc = -dc
    end

    for i = 1:maximum_length
        signal[i] = state[end] == 0 ? ndc : dc
        pushfirst!(state, feedback(state))
    end
    signal
end

# ╔═╡ 3d57368d-bf24-4fd0-9d8a-b0824d2f1616
begin
    fₛ = 100
    u = mlbs(5, 1)
    U = fft(u) / length(u)

    over_sampling_rate = 16
    oversampled = kron(u, ones(over_sampling_rate))
    Oversampled = fft(oversampled) / length(oversampled)

    N = length(u)
    freqs = (0:N-1) * fₛ / N
    oversampled_freqs = (0:over_sampling_rate*N-1) * fₛ / N
end;

# ╔═╡ 0d9edf18-15bc-4622-994c-62638b9bdefb
# ╠═╡ show_logs = false
begin
    p1 = plot(
        LinRange(0, 1, length(oversampled)),
        oversampled;
        xlabel = "t (s)",
        ylabel = "u(t)",
        title = "ZOH-MLBS",
    )
    p2 = plot(
        freqs,
        abs.(U[1:N]);
        line = :stem,
        marker = :dot,
        xlim = (1, 100),
        # ylim=(0, 0.21),
        xlabel = "Frequncey (Hz)",
        ylabel = "Amplitude (linear)",
        title = "ZOH-MLBS FFT",
    )

    p3 = plot(
        freqs,
        abs.(Oversampled);
        line = :stem,
        marker = :dot,
        xlim = (1, 100),
        ylim = (0, 0.21),
        xlabel = "Frequncey (Hz)",
        ylabel = "Amplitude (linear)",
        title = "Continous MLBS FFT",
    )

    p4 = plot(
        oversampled_freqs[1:N*over_sampling_rate÷2],
        abs.(Oversampled[1:N*over_sampling_rate÷2]);
        xlim = (1, 900),
        ylim = (0, 0.21),
        xlabel = "Frequncey (Hz)",
        ylabel = "Amplitude (linear)",
        title = "Continous MLBS FFT",
    )

    plot(p1, p2, p3, p4; size = (800, 500), legend = false)
end

# ╔═╡ Cell order:
# ╠═89c57041-7bf5-4b60-8263-92177041ca22
# ╠═62c7805d-f083-404d-8f4a-1a732d73920b
# ╟─bcc0c1b6-196c-452c-9a3d-d73cbd6fc93b
# ╠═3d57368d-bf24-4fd0-9d8a-b0824d2f1616
# ╟─0d9edf18-15bc-4622-994c-62638b9bdefb
# ╠═4399fb37-636a-4a01-b4af-69fc07769c1c
# ╠═c3269876-460c-469b-8cc7-bda1604657be

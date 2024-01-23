### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 92cfc11a-9fa0-11ee-0dae-55f4f9e89892
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ fa140b32-2192-491e-b044-bb55f65242b0
using Plots, FFTW

# ╔═╡ a053dcf7-09e3-491c-ab77-3eb2ecf6ba3a
md"# Exercise 14: Discretization in time (Choice of sampling frequency)

!!! purpose
	Show the effect of choosing low sampling frequency.
	
	- The sampling frequency (fₛ) should be at least two times the frequency of the sampled data otherwise the signal gets `aliased`.

	- In order to avoid aliasing signals with bandwidth larger than (fₛ/2), it's advised to filter the signal first with a lowpass filter (anti-aliasing filter).
"

# ╔═╡ ec888e0d-d19c-4588-9bcb-aa81a0482418
begin
    f₁ = 1
    f₂ = 9
    fₛ = 10
    A = 1
    T = 1

    # the time interval is [0, 1[
    interval = 0:T/1000:T-T/1000
    sampling_points = 0:1/fₛ:T-1/fₛ  # 1/fₛ sample in the second
    N = length(sampling_points)

    u₁ = map(t -> A * sin(2π * f₁ * t), interval)
    discretized_u₁ = map(t -> A * sin(2π * f₁ * t), sampling_points)
    U₁ = fft(discretized_u₁) / N

    u₂ = map(t -> -A * sin(2π * f₂ * t), interval)
    discretized_u₂ = map(t -> -A * sin(2π * f₂ * t), sampling_points)
    U₂ = fft(discretized_u₂) / N
end;

# ╔═╡ 2b46ca93-a8d2-4375-9bcf-df9b67f58a42
begin
    p1 = plot(interval, u₁; xlabel = "Time (s)", ylabel = "Amplitude", label = "Original")
    p1 = scatter!(sampling_points, discretized_u₁; label = "Sampled data")

    p2 = plot(
        [-fₛ, fₛ],
        [1, 1];
        line = :stem,
        marker = :circle,
        legend = false,
        xlims = (-fₛ - 1, fₛ + 1),
        xlabel = "Freq (Hz)",
    )
    p2 = plot!(
        [-f₁, f₁],
        [A / 2, A / 2];
        line = :stem,
        marker = :circle,
        legend = false,
        xlims = (-fₛ - 1, fₛ + 1),
    )
    plot(p1, p2)
end

# ╔═╡ e03abf74-d626-4650-a58a-548fab0f9dc1
plot(
    [
        scatter(
            sampling_points,
            [discretized_u₁, discretized_u₂];
            markershapes = reshape([:+, :x], 1, :),
            xlims = (-0.02, 1),
            labels = reshape(["sampling points for u₁", "sampling points for u₂"], 1, :),
            xlabel = "Time (s)",
            ylabel = "Amplitude",
        ),
        plot(
            fₛ * sampling_points,
            abs.(U₁);
            line = :stem,
            marker = :circle,
            legend = false,
            xlims = (-fₛ - 1, fₛ + 1),
            xlabel = "Freq (Hz)",
        ),
    ]...,
)

# ╔═╡ 835a5dc6-1d02-4dd5-868e-9fe002f537fc
begin
    p3 = plot(interval, u₂; xlabel = "Time (s)", label = "Original", ylabel = "Amplitude")
    p3 = scatter!(sampling_points, discretized_u₂; label = "Sampled data")

    p4 = plot(
        [-fₛ, fₛ],
        [1, 1];
        line = :stem,
        marker = :circle,
        legend = false,
        xlims = (-fₛ - 1, fₛ + 1),
        xlabel = "Freq (Hz)",
    )
    p4 = plot!(
        [-f₂, f₂],
        [A / 2, A / 2];
        line = :stem,
        marker = :circle,
        legend = false,
        xlims = (-fₛ - 1, fₛ + 1),
    )
    plot(p3, p4)
end

# ╔═╡ Cell order:
# ╠═92cfc11a-9fa0-11ee-0dae-55f4f9e89892
# ╠═fa140b32-2192-491e-b044-bb55f65242b0
# ╟─a053dcf7-09e3-491c-ab77-3eb2ecf6ba3a
# ╠═ec888e0d-d19c-4588-9bcb-aa81a0482418
# ╟─2b46ca93-a8d2-4375-9bcf-df9b67f58a42
# ╟─e03abf74-d626-4650-a58a-548fab0f9dc1
# ╟─835a5dc6-1d02-4dd5-868e-9fe002f537fc

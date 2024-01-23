### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 3b4420fe-a053-11ee-1fdf-a9dd111cdac1
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ 146a58f9-2848-414b-a87c-9eb15d1e3a08
using Plots, FFTW, DSP, StatsBase, Distributions

# ╔═╡ fe2809da-37d0-4cd2-9143-28ec17b4b786
md"# Exercise 20: Generation of a multisine with flat amplitude spectrum

!!! purpose
	Same as Exercise 16: shoow the leakge effect of measuring non integer number of periods (but this one for multi sine signal).
"

# ╔═╡ 1afece58-a648-4285-95e8-1d7a68bec830
begin
    N = 128
    interval = 0:N-1

    F = length(0.15N:0.35N)
    excitation_lines = 1 .+ floor.(0.15N:0.35N) .|> Int

    ϕ_rand = rand(Uniform(), N)
    Uᵣ = zeros(Complex{Float64}, N, 1)
    Uᵣ[excitation_lines] = exp.(im * ϕ_rand[excitation_lines])
    uᵣ = 2real(ifft(Uᵣ))
    uᵣ = uᵣ / std(uᵣ)
    M = 4
    extra = div(length(uᵣ), 4)
    u_long = vcat(repeat(uᵣ, M), uᵣ[begin:extra])
    @assert length(u_long) == length(uᵣ) * (M + 0.25)

    # windows
    Z(u, w) = begin
        z = w .* u
        abs.(fft(z)) / sqrt(length(u))
    end
    # rectangular window
    u_rect_window = Z(u_long, 1)
    u_hanning_window = Z(u_long, u_long |> length |> hanning)
end;

# ╔═╡ 691361e0-038b-41fc-9e75-2f496fab758f
plot(
    u_long;
    ylims = (-2.8, 2.8),
    xlabel = "Time (ms)",
    ylabel = "Amplitude",
    legend = false,
)

# ╔═╡ 1c238747-6435-415f-b6b3-903825aef7f0
md"""!!! danger
	I got different results compared to the book I need to get back to this exercise! 
	Mine has both sides of the DFT components.
![book result](https://i.imgur.com/Xsl1zdK.png)
"""

# ╔═╡ 7915f87e-dcd4-4800-a8c3-3b800bdc7cd2
begin
    plot(
        amp2db.(real(u_rect_window));
        seriestype = :scatter,
        # xlims=(0, 500),
        ylims = (-80, 20),
        xlabel = "Frequency (Hz)",
        ylabel = "Amplitude (dB)",
    )
    plot!(amp2db.(u_hanning_window); seriestype = :scatter)
end

# ╔═╡ Cell order:
# ╠═3b4420fe-a053-11ee-1fdf-a9dd111cdac1
# ╠═146a58f9-2848-414b-a87c-9eb15d1e3a08
# ╟─fe2809da-37d0-4cd2-9143-28ec17b4b786
# ╠═1afece58-a648-4285-95e8-1d7a68bec830
# ╟─691361e0-038b-41fc-9e75-2f496fab758f
# ╟─1c238747-6435-415f-b6b3-903825aef7f0
# ╟─7915f87e-dcd4-4800-a8c3-3b800bdc7cd2

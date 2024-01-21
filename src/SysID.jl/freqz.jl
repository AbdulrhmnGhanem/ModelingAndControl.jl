### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 2bd75bfc-b880-11ee-1011-0d2b678384d8
# ╠═╡ show_logs = false
begin
    # If you are running this notebook as a stannalone notebook disable this cell.
    import Pkg
    Pkg.activate(joinpath("..", ".."))
end

# ╔═╡ bfeca572-9a85-49f3-9fc3-e50ffd781f27
using DSP, Plots

# ╔═╡ a61c53cc-a192-4e1d-9827-a5e599ea7263
"""
Similar to [MATLAB's Freqz](https://www.mathworks.com/help/signal/ref/freqz.html)

- filt: Filter coofficients
- fₛ: sampling frequency
- lower: lower limit for frequncy band of interest
- upper: upper limit for frequecny band of interest
- NPeriod: number of data points per block
"""
function freqz(filt::FilterCoefficients, fₛ::Int, lower::Int, upper::Int, NPeriod::Int)
    frequencies = range(lower, upper, length=NPeriod)
    H = freqresp(filt, frequencies * 2π / fₛ)
    return frequencies, H
end

# ╔═╡ f285708f-916a-489a-a24e-0b5a7da1d0bb
function cheby1(n, r, wp)
    digitalfilter(Lowpass(wp), Chebyshev1(n, r))
end

# ╔═╡ ee226799-88ab-48a1-881e-b76a20f373d7
begin
    fₛ = 128
    NPeriod = 1024
    filter_order = 2
    filter_resonance = 10
    filter_cutoff = 0.1fₛ

    h = cheby1(filter_order, filter_resonance, 2filter_cutoff / fₛ)
    w, H = freqz(h, fₛ, 1, 60, NPeriod)
end;

# ╔═╡ f4a45d38-eefa-420a-acf9-952fd58675a0
plot(w, amp2db.(abs.(H));
    ylims=(-60, 1)
)

# ╔═╡ Cell order:
# ╠═2bd75bfc-b880-11ee-1011-0d2b678384d8
# ╠═bfeca572-9a85-49f3-9fc3-e50ffd781f27
# ╠═a61c53cc-a192-4e1d-9827-a5e599ea7263
# ╠═f285708f-916a-489a-a24e-0b5a7da1d0bb
# ╠═ee226799-88ab-48a1-881e-b76a20f373d7
# ╠═f4a45d38-eefa-420a-acf9-952fd58675a0

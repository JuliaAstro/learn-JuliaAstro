### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 7fe5fdc0-d4c6-11f0-a48a-035218fa12ed
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using CairoMakie, DynamicQuantities
end

# ╔═╡ 94c34cd5-e77c-4e0e-8546-5c2ddd9d6ed2
lines(rand(10)u"m")

# ╔═╡ Cell order:
# ╠═7fe5fdc0-d4c6-11f0-a48a-035218fa12ed
# ╠═94c34cd5-e77c-4e0e-8546-5c2ddd9d6ed2

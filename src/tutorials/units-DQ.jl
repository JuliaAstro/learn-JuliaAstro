### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "Using units in astrophysical calculations"
#> layout = "layout.jlhtml"
#> date = "2025-11-19"
#> description = "Work with units in astrophysical calculations."
#> tags = ["units", "plots"]
#> 
#>     [[frontmatter.author]]
#>     name = "Ian Weaver"
#>     url = "https://github.com/icweaver"

using Markdown
using InteractiveUtils

# ╔═╡ fd88a6c1-0abe-4a5a-9414-bb15730c9d18
begin
	import Pkg
	Pkg.add(; url = "https://github.com/MakieOrg/Makie.jl", subdir = "Makie", rev = "ff/dim-converts")
	Pkg.add(["PlutoUI", "DynamicQuantities", "StatsBase", "Distributions", "CairoMakie", "Makie", "DimensionalData"])

	using DynamicQuantities: @u_str, @us_str, dimension, uconvert, ustrip
	using DynamicQuantities.Constants: pc, G
	using StatsBase: mean
	using Distributions: Normal
	using CairoMakie: Colorbar, stephist, heatmap
	using Makie: DQConversion
	using DimensionalData: DimArray, val, dims
end

# ╔═╡ 17c6b7df-a8b3-45d1-9491-526afce11318
using PlutoUI: TableOfContents, details

# ╔═╡ c6ad0267-65d1-4372-a538-22acd9b5d02b
md"""
# Using units in astrophysical calculations

*Authors: Ian Weaver*

This notebook is modified from <https://learn.astropy.org/tutorials/quantities.html>

!!! tip ""
	## Learning Goals
	
	* 

!!! note ""
	## Keywords

	units, plots

!!! warning ""
	## Summary
	
	
"""

# ╔═╡ 05b485e7-115a-4dbb-aa73-0ca6ace2f5c0
md"""
### Imports
"""

# ╔═╡ 60c89d86-942e-4c97-bd7a-ad2f792b1155
md"""
## 1. Galaxy mass
"""

# ╔═╡ 1d2293bd-a236-4b41-a9d7-9c27463b5062
Reff = 29 * u"Constants.pc" # Or 29 * pc

# ╔═╡ 9b6e871f-0a49-4a24-a9b0-51e3008d6db4
ustrip(Reff)

# ╔═╡ 1b239188-bba1-44d1-bc9c-10b60f762e0d
ustrip(pc, Reff)

# ╔═╡ cecd879b-2d95-4150-8caf-c327165ddec6
dimension(Reff)

# ╔═╡ cc736aec-c48a-411a-af83-4255309d77e9
md"""
Furthermore, we can convert the radius to any other unit of length. Here, we convert it to kilometers:
"""

# ╔═╡ 2c7a71f6-e618-40db-80b7-1ad5f09e59d5
Reff |> us"km" # Or uconvert(us"km", Reff)

# ╔═╡ c5b4f340-c774-4f09-af4c-f326afce5de3
md"""
Synth. dataset of rad vels
"""

# ╔═╡ 892f19a6-af21-4353-9987-9de795ba7ad7
v = rand(Normal(206, 4.3), 500)u"km/s"

# ╔═╡ 808f338d-7faa-4d97-ad7b-a7b4066c83f8
first(v, 10) .|> us"km/s"

# ╔═╡ 59af301c-0713-4c93-a824-7375f8c4f761
stephist(v |> us"km/s")

# ╔═╡ d8a2219e-1094-4eb7-be34-2ce58b6bd462
md"""
!!! todo
	Units + histograms (and more!) working in <https://github.com/MakieOrg/Makie.jl/pull/5323>.
"""

# ╔═╡ 2d236d32-faf4-41ef-84de-1d40c02fb238
sigma = (sqrt ∘ sum)((v .- mean(v)).^2 / length(v))

# ╔═╡ 64bfb535-ab84-411c-9a1b-3ed8778e8516
sigma |> us"km/s"

# ╔═╡ a3baf049-5419-4773-8b7f-0ba07a9f1728
M = 4 * sigma^2 * Reff / G

# ╔═╡ bb2622f9-98a5-48ca-92cb-bbce84d797f1
M .|> (us"Constants.M_sun", us"g")

# ╔═╡ 9e13bde0-eccf-4d21-b603-fd186e87b7d0
(log10 ∘ ustrip)(us"Constants.M_sun", M)

# ╔═╡ 5a9dfdf7-aaf8-421d-a451-6a6a5c35e1cc
md"""
Note that this is different than:
"""

# ╔═╡ a13de611-145a-40da-9414-8f3f8f85ad98
(log10 ∘ ustrip)(M)

# ╔═╡ a3f42e7e-af38-46b0-b1fb-bdbf77312e16
md"""
emphasizing the importance of being explicit with our units. Similarly, taking the logarithm of something with units is not mathematically well defined, so this will sensibly error as well:
"""

# ╔═╡ 095e7efe-f4c6-4a8a-a029-03ce97cd15bd
log10(M)

# ╔═╡ 2648b454-2762-4941-b13c-2303ffcd6521
let
	sol = details("Solution",
	md"""
	```julia
	using DynamicQuantities: Constants as C
	
	# Speed from Kepler's law
	v_kep = sqrt(C.G * C.M_sun / C.au)
	
	# View in units of km/s (≈ 29.7847 km/s)
	v_kep |> us"km/s"
	
	# Speed from kinematics (≈ 29.7853 km/s)
	v_kin = 2π * C.au / u"yr"
	
	# Percent difference (≈ 0.002%)
	percent_diff = 100 * (v_kin - v_kep) / v_kep
	```
	""")
	
	md"""
	!!! warn "Exercise"
	
		Use DynamicQuantities.jl and Kepler's law in the form given below to determine the (circular) orbital speed of the Earth around the Sun in km/s:
		
		```math
		v = \sqrt{\frac{G M_⊙}{r}}
		```
	
		No need to look up constants or conversion factors to do this calculation -- it's all in `DynamicQuantities.Units` and `DynamicQuantities.Constants`.
	
		There's a much easier way to figure out the velocity of the Earth using just two units or quantities. Do that and then compare to the Kepler's law answer (the easiest way is probably to compute the percentage difference, if any).

		$(sol)

		Completely optional, but a good way to convince ourselves of the value of DynamicQuantities.jl: Do the above calculations by hand. Look up all the appropriate conversion factors and use paper-and-pencil / basic approaches for keeping track of them all. Using Julia as a basic calculator is also fine. Which one took longer?
	"""
end

# ╔═╡ 48fac714-ab96-4475-ad0b-0c61432bf849
# Enter here

# ╔═╡ 9e8eabc7-f1ce-4dbb-bfd5-e6d28793f9d6
md"""
## 2. Molecular cloud mass

In this second example, we will demonstrate how using units can facilitate a full derivation of the total mass of a molecular cloud using radio observations of isotopes of Carbon Monoxide (CO).
"""

# ╔═╡ 78740b86-45e9-45f5-b4c2-be6cb65f1368
md"""
### Setting up the data cube

Let's assume that we've mapped the inner part of a molecular cloud in the ``J = 1 - 0`` rotational transition of ``\text{C}^{18}\text{O}`` and are interested in measuring its total mass. The measurement produced a data cube with RA and Dec as spatial coordiates and velocity as the third axis. Each voxel in this data cube represents the brightness temperature of the emission at that position and velocity. Furthermore, we'll assume that we have an independent measurement of distance to the cloud ``d = 250\text{ pc}`` and that the excitation temperature is known and constant throughout the cloud: ``T_\text{ex} = 25\text{ K}``:
"""

# ╔═╡ 2c76c406-d15c-4271-baa4-ebebf5299429
d = 250u"Constants.pc"

# ╔═╡ cb5b1020-b1f8-4ea1-ad03-a91b5c3ab0c2
T_ex = 25u"K"

# ╔═╡ b86a31df-9de1-4e27-8ca7-ce7f3d2576fa
md"""
We'll generate a synthetic dataset, assuming the cloud follows a Gaussian distribution in each of RA, Dec, and velocity. We start by creating a 100×100×300 array, such that the first coordinate is right ascension, the second is declination, and the third is velocity. We use the numpy.meshgrid function to create data cubes for each of the three coordinates, and then use them in the formula for a Gaussian to generate an array with the synthetic data cube. In this cube, the cloud is positioned at the center of the cube, with _σ_ and the center in each dimension shown below. Note in particular that the _σ_ for RA and Dec have different units from the center, but astropy automatically does the relevant conversions before computing the exponential.
"""

# ╔═╡ aab4113d-d7ad-4521-a208-e192ad218cf0
begin
	# 1D coordinate quantities
	ras = range(52, 52.5; length = 100)u"deg"
	decs = range(0, 0.5; length = 100)u"deg"
	vs = range(0, 30; length = 300)u"km/s"
end;

# ╔═╡ 539290af-787d-4deb-928c-50e4e9f28173
data = let
    # Cloud's center
    cen_ra = 52.25u"deg"
    cen_dec = 0.25u"deg"
    cen_v = 15u"km/s"

    # Cloud's size
    sig_ra = 3u"arcmin"
    sig_dec = 4u"arcmin"
    sig_v = 3u"km/s"

	A = [
	    exp(
	        -0.5 * ((ra - cen_ra) / sig_ra)^2
	        -0.5 * ((dec - cen_dec) / sig_dec)^2
	        -0.5 * ((v - cen_v) / sig_v)^2
	    )
	    for ra in ras, dec in decs, v in vs
	]

	DimArray(A * u"K", (RA = ras, Dec = decs, Vel = vs))
	# DimArray(A * u"K", (:RA, :Dec, :Vel))
end

# ╔═╡ e86d6cf7-2274-403a-b1db-017973f33fb7
md"""
!!! note
	The units of the exponential are dimensionless, so we multiplied the data cube by K to get brightness temperature units. As an aside for experts, we're setting up our artificial cube on the main-beam temperature scale ``\left(T_\text{MB}\right)`` which is the closest we can normally get to the actual brightness temperature of our source.
"""

# ╔═╡ 00b9d37f-ce7b-4491-b1df-f0963f2598a8
md"""
We will also need to know the width of each velocity bin and the size of each pixel, so let's calculate that now:
"""

# ╔═╡ 1f6a0d6c-832b-474d-8fbf-a9de6e821d80
# Average pixel size
# This is only right if dec ~ 0, because of the cos(dec) factor.
Δra = (maximum(ras) - minimum(ras)) / length(ras) # Typed |Delta<TAB>

# ╔═╡ f80f8b06-976f-4b43-8074-8d0bf725ba48
Δdec = (maximum(decs) - minimum(decs)) / length(decs)

# ╔═╡ 6fa07efa-fd25-4bca-bd18-8ba5c557d146
# Average velocity bin width
Δv = (maximum(vs) - minimum(vs)) / length(vs)

# ╔═╡ 3211ba19-0140-450d-95c6-698a4f0ccbff
md"""
Note that DynamicQuantities.jl uses the unitless radian by default. We can easily display this in our desired unit system:
"""

# ╔═╡ 175854bb-4385-4199-ada9-358def53a822
(Δra, Δdec) .|> us"arcsec" # Display in arcseconds

# ╔═╡ 7ab2fe80-812f-4a21-933a-e169a17f6c32
md"""
We're interested in the integrated intensity over all of the velocity channels, so let's create a 2D quantity array by summing our data cube along the velocity axis (multiplying by the velocity width of a pixel):
"""

# ╔═╡ 6ddeb42d-37a5-48a3-ac82-47e4ec2ca541
intcloud = reduce(+, eachslice(data * Δv; dims = :Vel))

# ╔═╡ 4d4794fb-804d-4b4a-8c32-a32d88e43e30
md"""
!!! todo
	Get `Base.sum` support for DQ. Discussion here: <https://github.com/JuliaPhysics/DynamicQuantities.jl/issues/76#issuecomment-3614719247>
"""

# ╔═╡ 1fa8010d-9a5d-4297-9665-9fa8795ef5f7
md"""
!!! note
	Radio astronomers use a rather odd set of units [K km/s] for integrated intensity (that is, summing all the emission from a line over velocity).
"""

# ╔═╡ f2b222ec-0783-487e-9c52-835976a555b6
let
	A = intcloud
	x, y = dims(A)
	u_A = us"K*km/s"
	
	fig, ax, p =  heatmap(
		val(x) .|> us"deg",
		val(y) .|> us"deg",
		ustrip.(u_A, parent(A))
	)

	ax.xlabel = "RA"
	ax.ylabel = "Dec"
	
	Colorbar(fig[1, 2], p; label = string("Intensity [", dimension(u_A), " ]"))
	
	fig
end

# ╔═╡ e4232b15-3369-438e-994b-042aab477a7f
md"""
!!! todo
	See if something like this can be upstreamed to DD. Would be nice to be able to just do:

	```julia
	heatmap(intcloud)
	```
"""

# ╔═╡ 59b4d441-9a74-468f-ad8c-882516a09049
md"""
# Notebook setup 🔧
"""

# ╔═╡ bedc8ccd-e6f6-4dd1-a0b6-1889f4b5b658
TableOfContents()

# ╔═╡ Cell order:
# ╟─c6ad0267-65d1-4372-a538-22acd9b5d02b
# ╟─05b485e7-115a-4dbb-aa73-0ca6ace2f5c0
# ╠═fd88a6c1-0abe-4a5a-9414-bb15730c9d18
# ╟─60c89d86-942e-4c97-bd7a-ad2f792b1155
# ╠═1d2293bd-a236-4b41-a9d7-9c27463b5062
# ╠═9b6e871f-0a49-4a24-a9b0-51e3008d6db4
# ╠═1b239188-bba1-44d1-bc9c-10b60f762e0d
# ╠═cecd879b-2d95-4150-8caf-c327165ddec6
# ╟─cc736aec-c48a-411a-af83-4255309d77e9
# ╠═2c7a71f6-e618-40db-80b7-1ad5f09e59d5
# ╟─c5b4f340-c774-4f09-af4c-f326afce5de3
# ╠═892f19a6-af21-4353-9987-9de795ba7ad7
# ╠═808f338d-7faa-4d97-ad7b-a7b4066c83f8
# ╠═59af301c-0713-4c93-a824-7375f8c4f761
# ╟─d8a2219e-1094-4eb7-be34-2ce58b6bd462
# ╠═2d236d32-faf4-41ef-84de-1d40c02fb238
# ╠═64bfb535-ab84-411c-9a1b-3ed8778e8516
# ╠═a3baf049-5419-4773-8b7f-0ba07a9f1728
# ╠═bb2622f9-98a5-48ca-92cb-bbce84d797f1
# ╠═9e13bde0-eccf-4d21-b603-fd186e87b7d0
# ╟─5a9dfdf7-aaf8-421d-a451-6a6a5c35e1cc
# ╠═a13de611-145a-40da-9414-8f3f8f85ad98
# ╟─a3f42e7e-af38-46b0-b1fb-bdbf77312e16
# ╠═095e7efe-f4c6-4a8a-a029-03ce97cd15bd
# ╟─2648b454-2762-4941-b13c-2303ffcd6521
# ╠═48fac714-ab96-4475-ad0b-0c61432bf849
# ╟─9e8eabc7-f1ce-4dbb-bfd5-e6d28793f9d6
# ╟─78740b86-45e9-45f5-b4c2-be6cb65f1368
# ╠═2c76c406-d15c-4271-baa4-ebebf5299429
# ╠═cb5b1020-b1f8-4ea1-ad03-a91b5c3ab0c2
# ╟─b86a31df-9de1-4e27-8ca7-ce7f3d2576fa
# ╠═aab4113d-d7ad-4521-a208-e192ad218cf0
# ╠═539290af-787d-4deb-928c-50e4e9f28173
# ╟─e86d6cf7-2274-403a-b1db-017973f33fb7
# ╟─00b9d37f-ce7b-4491-b1df-f0963f2598a8
# ╠═1f6a0d6c-832b-474d-8fbf-a9de6e821d80
# ╠═f80f8b06-976f-4b43-8074-8d0bf725ba48
# ╠═6fa07efa-fd25-4bca-bd18-8ba5c557d146
# ╟─3211ba19-0140-450d-95c6-698a4f0ccbff
# ╠═175854bb-4385-4199-ada9-358def53a822
# ╟─7ab2fe80-812f-4a21-933a-e169a17f6c32
# ╠═6ddeb42d-37a5-48a3-ac82-47e4ec2ca541
# ╟─4d4794fb-804d-4b4a-8c32-a32d88e43e30
# ╟─1fa8010d-9a5d-4297-9665-9fa8795ef5f7
# ╠═f2b222ec-0783-487e-9c52-835976a555b6
# ╟─e4232b15-3369-438e-994b-042aab477a7f
# ╟─59b4d441-9a74-468f-ad8c-882516a09049
# ╠═17c6b7df-a8b3-45d1-9491-526afce11318
# ╠═bedc8ccd-e6f6-4dd1-a0b6-1889f4b5b658

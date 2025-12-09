### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "Using units in astrophysical calculations"
#> layout = "layout.jlhtml"
#> date = "2025-11-19"
#> description = "Work with units in astrophysical calculations."
#> tags = ["units", "plots", "radio astronomy", "data cubes"]
#> 
#>     [[frontmatter.author]]
#>     name = "Ian Weaver"
#>     url = "https://github.com/icweaver"

using Markdown
using InteractiveUtils

# ╔═╡ fd88a6c1-0abe-4a5a-9414-bb15730c9d18
begin
	import Pkg
	Pkg.activate(Base.current_project())
    Pkg.instantiate()
	
	using DynamicQuantities: @u_str, @us_str, dimension, uconvert, ustrip
	using DynamicQuantities.Constants: G, pc, h, k_B, c as c_0
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

This notebook is modified from <https://learn.astropy.org/tutorials/quantities.html>

!!! tip ""
	## Learning Goals
	
	- Use `Quantity` objects to estimate a hypothetical galaxy's mass
	- Take advantage of constants in the DynamicalQuantities.jl package
	- Print formatted unit strings
	- Plot `Quantity` objects with unit labels, using Makie.jl
	- Do math with `Quantity` objects
	- Convert quantities
	- Convert between wavelength and energy
	- Write functions that take `Quantity` objects instead of plain arrays
	- Make synthetic radio observations
	- Use `Quantity` objects such as data cubes to facilitate a full derivation of the total mass of a molecular cloud

!!! note ""
	## Keywords

	units, plots, radio astronomy, data cubes

!!! warning ""
	## Summary
	
	In this tutorial we present some examples showing how `Quantity` objects can make astrophysics calculations easier. The examples include calculating the mass of a galaxy from its velocity dispersion and determining masses of molecular clouds from ``\mathrm{CO}`` intensity maps. We end with an example of good practices for using quantities in functions you might distribute to other people.
"""

# ╔═╡ 05b485e7-115a-4dbb-aa73-0ca6ace2f5c0
md"""
### Imports
"""

# ╔═╡ 60c89d86-942e-4c97-bd7a-ad2f792b1155
md"""
## 1. Galaxy mass

In this first example, we will use Quantity objects to estimate a hypothetical galaxy's mass, given its half-light radius and radial velocities of stars in the galaxy.

Let's assume that we measured the half-light radius of the galaxy to be 29 pc projected on the sky at the distance of the galaxy. This radius is often called the "effective radius", so we'll store it as a `Quantity` object with the name `Reff`. The easiest way to create a `Quantity` object is by multiplying the value with its unit:
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
	sol = details("Example solution",
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

Let's assume that we've mapped the inner part of a molecular cloud in the ``J = 1 - 0`` rotational transition of ``\mathrm{C}^{18}\mathrm{O}`` and are interested in measuring its total mass. The measurement produced a data cube with RA and Dec as spatial coordiates and velocity as the third axis. Each voxel in this data cube represents the brightness temperature of the emission at that position and velocity. Furthermore, we'll assume that we have an independent measurement of distance to the cloud ``d = 250\, \mathrm{pc}`` and that the excitation temperature is known and constant throughout the cloud: ``T_\text{ex} = 25\, \mathrm{K}``:
"""

# ╔═╡ 2c76c406-d15c-4271-baa4-ebebf5299429
d = 250u"Constants.pc"

# ╔═╡ cb5b1020-b1f8-4ea1-ad03-a91b5c3ab0c2
T_ex = 25u"K"

# ╔═╡ b86a31df-9de1-4e27-8ca7-ce7f3d2576fa
md"""
We'll generate a synthetic dataset, assuming the cloud follows a Gaussian distribution in each of RA, Dec, and velocity. We start by creating a 100×100×300 array, such that the first coordinate is right ascension, the second is declination, and the third is velocity. In this data cube, the cloud is positioned at the center, with ``\sigma`` and the center in each dimension shown below. Note in particular that the ``\sigma`` for RA and Dec have different units from the center, but DynamicQuantities.jl automatically does the relevant conversions before computing the exponential.
"""

# ╔═╡ 9a002059-650b-4dc1-9395-6f266ed35500
md"""
!!! todo
	Explain. What does this part mean?

	> Note in particular that the ``\sigma`` for RA and Dec have different units from the center, but DynamicQuantities.jl automatically does the relevant conversions before computing the exponential.
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
end

# ╔═╡ e86d6cf7-2274-403a-b1db-017973f33fb7
md"""
!!! note
	The units of the exponential are dimensionless, so we multiplied the data cube by ``\mathrm{K}`` to get brightness temperature units. As an aside for experts, we're setting up our artificial cube on the main-beam temperature scale ``\left(T_\text{MB}\right)``, which is the closest we can normally get to the actual brightness temperature of our source.
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
Note that DynamicQuantities.jl uses the unitless radian by default. We can still easily display this in our desired unit system:
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
	Radio astronomers use a rather odd set of units ``[\mathrm{K\, km/s}]`` for integrated intensity (that is, summing all the emission from a line over velocity).
"""

# ╔═╡ 4ee79c7c-5c4e-457b-b713-103937e50355
md"""
We can plot the 2D quantity using Makie's `heatmap` function:
"""

# ╔═╡ f2b222ec-0783-487e-9c52-835976a555b6
let
	A = intcloud
	x, y = dims(A)
	u_A = us"K*km/s"
	
	fig, ax, p =  heatmap(
		val(x) .|> us"deg",
		val(y) .|> us"rad",
		ustrip.(u_A, parent(A));
		colormap = :cividis,
	)

	ax.xreversed = true
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

# ╔═╡ c749ce2d-17ae-45f4-b721-3f486b1cbc23
md"""
### Measuring The Column Density of ``\mathrm{CO}``

In order to calculate the mass of the molecular cloud, we need to measure its column density. A number of assumptions are required for the following calculation; the most important are that the emission is optically thin (typically true for ``\mathrm{C}^{18}\mathrm{O}``), and that conditions of local thermodynamic equilibrium hold along the line of sight. In the case where the temperature is large compared to the separation in energy levels for a molecule and the source fills the main beam of the telescope, the total column density for ``\mathrm{C}^{13}\mathrm{O}`` is:

```math
N = C \frac{\int T_\text{B}(V) \, dV}{1 - e^{-B}}\ ,
```

where ``T_\text{B}`` is the brightness temperature, and the constants ``C`` and ``B`` are given by:

```math
\begin{align*}
C &= 3.0 \times {10}^{14}\ \mathrm{K^{-1}\, cm^{-2}\, km^{-1}\, s}
	\left(\frac{\nu}{\nu_{13}}\right)^2 \frac{A_{13}}{A} \\
B &= \frac{h\nu}{k_\text{B} T}
\end{align*}
```

(Rohlfs & Wilson [Tools for Radio Astronomy](https://www.springer.com/gp/book/9783662053942)).
"""

# ╔═╡ 3a3328d5-31de-4deb-bff7-d25c1fcbc4ef
md"""
Here we have given an expression for ``C`` scaled to the values for ``\mathrm{C}^{13}\mathrm{O}`` (``\nu_{13}`` and ``A_{13}``). In order to use this relation for ``\mathrm{C}^{18}\mathrm{O}``, we need to rescale the frequencies ``\nu`` and the Einstein coefficients (``A``). Lastly, ``C`` is in funny mixed units, but that's okay. We'll be able to do our unit handling in the usual way.
"""

# ╔═╡ fb0f2941-2e56-45db-91fe-e8d6744a00e0
md"""
First, we look up the wavelength for these emission lines and store them as quantities:
"""

# ╔═╡ bf618161-ef96-4445-8fc8-25dc5f662242
const λ_13 = 2.60076u"mm"

# ╔═╡ 1a16d2e1-7898-415c-9a98-a1f2f87e08be
const λ_18 = 2.73079u"mm"

# ╔═╡ e0a2c745-41cf-4359-bb99-3117fbb507cc
md"""
And compute their corresponding frequencies:
"""

# ╔═╡ d45bad3e-4c82-436b-8e64-c61e5cf65c2f
λ_to_ν(λ) = c_0 / λ

# ╔═╡ ee46b865-dbd4-4939-b241-514941dd138d
const ν_13 = λ_13 |> λ_to_ν

# ╔═╡ b183c219-a39b-427a-b740-5675fcc175ca
const ν_18 = λ_18 |> λ_to_ν

# ╔═╡ 323d1169-c59f-456c-9b19-093e26f84214
md"""
!!! todo
	See how helpful UnitfulEquivalencies.jl-like functionality for DQ would be.
"""

# ╔═╡ d0215082-52df-4eae-a409-7c3c7cfc69ed
md"""
Next, we look up Einstein coefficients (in units of s⁻¹), and calculate the ratios in constant ``C``. Note how the ratios of frequency and Einstein coefficient units are dimensionless, so the unit of _C_ is unchanged.
"""

# ╔═╡ 3a24f4aa-b074-4beb-b77d-6778e2fe580a
const A_13 = 7.4e-8 / u"s"

# ╔═╡ 1e203abc-4a02-4d01-bf9e-9df636e70ebd
const A_18 =  8.8e-8 / u"s"

# ╔═╡ 98869c25-2644-47d1-b8cc-05699292f2a8
C = 3e14u"s/(K*cm^2*km)" * (ν_18/ν_13)^3 * (A_13/A_18)

# ╔═╡ 18985bb2-6cce-4417-97a2-c2da7b9e0428
C |> us"s / K / km / cm^2"

# ╔═╡ d9e2b828-9c32-4638-bfcd-0c16b221aa43
md"""
Now we move on to calculate the constant ``B``. This is given by the ratio of ``\dfrac{hν}{k_\text{B}T}``, where ``h`` is Planck's constant, ``k_\text{B}`` is the Boltzmann's constant, ``ν`` is the emission frequency, and ``T`` is the excitation temperature. The constants were imported from `DynamicQuantities.Constants`, and the other two values are already calculated, so here we just take the ratio:
"""

# ╔═╡ 4092a893-818e-49f4-93d0-be7bd652dddc
B = h * ν_18 / (k_B * T_ex)

# ╔═╡ c5c907d0-a669-43f3-8030-0bd67452f0a1
k_B

# ╔═╡ a0b64ba9-cbec-404e-bf39-aee03ae407ae
md"""
Note how DynamicQuantities.jl intelligently cancelled the units for us, while still keeping this as a Quantity object:
"""

# ╔═╡ 82f11626-722b-46fa-9173-8b5d7a80a190
typeof(B)

# ╔═╡ 6aff587f-be6f-4fd8-96df-9c12f3769f32
md"""
At this point we have all the ingredients to calculate the number density of ``\mathrm{CO}`` molecules in this cloud. We already integrated (summed) over the velocity channels above to show the integrated intensity map, but we'll do it again here for clarity. This gives us the column density of ``\mathrm{CO}`` for each spatial pixel in our map. We can then print out the peak column column density:
"""

# ╔═╡ ef498536-8fb8-46e3-9d8d-f7eb3994a3f5
NCO = C * reduce(+, eachslice(data * Δv; dims = :Vel)) / (1 - exp(-B))

# ╔═╡ 25598fcd-761d-40b8-95ac-8b68a00026da
md"""
!!! note ""
	**Peak CO Column density:** $(maximum(NCO) |> us"cm^-2")
"""

# ╔═╡ eba2f06f-ebe9-492d-81d2-1cc4fccd5b0a
md"""
### ``\mathrm{CO}`` to Total Mass

We are using ``\mathrm{CO}`` as a tracer for the much more numerous ``\mathrm{H}_2``, the quantity we are actually trying to infer. Since most of the mass is in ``\mathrm{H}_2``, we calculate its column density by multiplying the ``\mathrm{CO}`` column density with the (known/assumed) ``\mathrm{H}_2 / \mathrm{CO}`` ratio:
"""

# ╔═╡ de22c778-1529-4177-85e6-a0178f437a8c
H₂_CO_ratio = 5.9e6

# ╔═╡ c262aa56-b9b6-4da4-8810-d6120f0724c6
NH₂ = NCO * H₂_CO_ratio

# ╔═╡ 1e8ec5c5-3d41-44a7-8dfa-81d981750d9e
md"""
!!! note ""
	**Peak ``\mathrm{H}_2`` column density:** $(maximum(NH₂) |> us"cm^-2")
"""

# ╔═╡ c6b35992-973b-49d4-bfa7-855e8fe10924
md"""
That's a peak column density of roughly 50 magnitudes of visual extinction (assuming the conversion between ``N_{\mathrm{H}_2}`` and ``A_V`` from Bohlin et al. 1978), which seems reasonable for a molecular cloud.

We obtain the mass column density by multiplying the number column density by the mass of an individual ``\mathrm{H_2}`` molecule:
"""

# ╔═╡ 85f84c90-46e6-4d3f-bba2-6fc190a05533
mH₂ = 2*1.008u"Constants.u"

# ╔═╡ 16cdd7e8-0ec2-4953-ad8b-f53669d05eda
mH₂ |> us"Constants.u"

# ╔═╡ 522d88d6-5f60-401b-8786-0236c0859eda
ρ = NH₂ * mH₂

# ╔═╡ fbd0d75b-b3af-4eb5-a249-423317d656e5
md"""
A final step in going from the column density to mass is summing up over the area. If we do this in the straightforward way of length × width of a pixel, this area is then in units of deg²:
"""

# ╔═╡ 0eefe540-1669-4510-bef7-8fb7f8be68f1
Δap = Δra * Δdec

# ╔═╡ 4f59f4e6-c735-4a40-b171-744b17eb37fd
Δap |> us"deg^2"

# ╔═╡ b057812d-1a8e-479c-8192-1c1acdb97d23
Δa = Δap * d^2

# ╔═╡ 586b1993-be11-4437-9d3b-789511ea26bf
Δa |> us"cm^2"

# ╔═╡ e81f3b44-bc04-4768-81c2-c2733958f5a9
md"""
Finally, multiplying the column density with the pixel area and summing over all the pixels gives us the cloud mass:
"""

# ╔═╡ 88966b7e-b261-45a6-9f12-6b74f24c03e4
M_cloud = sum(ρ * Δa)

# ╔═╡ ff7def67-eb96-45b3-959a-2b01c7804e07
M_cloud |> us"Constants.M_sun"

# ╔═╡ 965a8e44-e96e-434e-8454-cced531ae9d2
md"""
!!! note ""
	**Total cloud mass:** $(M_cloud |> us"Constants.M_sun")
"""

# ╔═╡ 802f5cad-a0c0-4426-94f2-426f89dea7e1
md"""
!!! tip "Exercises"
	The astro material was pretty heavy on that one, so let's focus on some associated statistics using DynamicQuantities.jl's array capabililities. Compute the median and mean of the data with the `mean` and `median` functions. Why are their values so different?

	Similarly, compute the standard deviation and variance. Do they have the units you expect?
"""

# ╔═╡ 2eed0bc8-2306-42bb-9803-33ae131021e6
md"""
## 3. Using `Quantities` with functions

`Quantity` is also a useful tool if you plan to share some of your code, either with collaborators or the wider community. By writing functions that take `Quantity` objects instead of raw numbers or arrays, you can write code that is agnostic to the input unit. In this way, you may even be able to prevent [the destruction of Mars orbiters](http://en.wikipedia.org/wiki/Mars_Climate_Orbiter#Cause_of_failure). Below, we provide a simple example.

Suppose you are working on an instrument, and the person funding it asks for a function to give an analytic estimate of the response function. You determine from some tests it's basically a Lorentzian, but with a different scale along the two axes. Your first thought might be to do this:
"""

# ╔═╡ 9922bdbe-cf9e-486a-8e8a-28a7ded12d04
function response_func_bad(xinarcsec, yinarcsec)
    xscale = 0.9
    yscale = 0.85
    xfactor = 1 / (1 + xinarcsec / xscale)
    yfactor = 1 / (1 + yinarcsec / yscale)

    return xfactor * yfactor
end

# ╔═╡ 118aff27-dc6d-4686-b64a-ff840b731030
md"""
You meant the inputs to be in arcsec, but alas, you send that to your collaborator and they don't look closely and think the inputs are instead supposed to be in arcmin. So they do:
"""

# ╔═╡ ffc4bdc3-16e3-4b3e-b448-cfdc16428378
response_func_bad(1.0, 1.2)

# ╔═╡ 71cca345-cffc-4aee-a671-08901a92babe
md"""
And now they tell all their friends how terrible the instrument is, because it's supposed to have arcsecond resolution, but your function clearly shows it can only resolve an arcmin at best. But you can solve this by requiring they pass in Quantity objects. The new function could simply be:
"""

# ╔═╡ 1c3faf9f-adcf-4232-acf4-655d828623e3
function response_func_good(x, y)
    xscale = 0.9us"arcsec" # We use symbolic dimensions here
    yscale = 0.85us"arcsec" # to treat radians as unitful quantities
    xfactor = 1 / (1 + x / xscale)
    yfactor = 1 / (1 + y / yscale)

    return xfactor * yfactor
end

# ╔═╡ 7cfb17fc-81a4-4522-a88b-2fd10001baaa
md"""
And your collaborator now has to pay attention. If they just blindly put in a number, they get an error:
"""

# ╔═╡ ce1804da-df61-4e4c-aa3b-d950e83b7c13
response_func_good(1.0, 1.2)

# ╔═╡ 5590c909-6103-4427-b1a2-4cf35265dc07
md"""
Which is their cue to provide the units explicitly:
"""

# ╔═╡ f7486c08-b74b-471f-bbc2-690439a487f8
response_func_good(1.0u"arcmin", 1.2u"arcmin")

# ╔═╡ c143ab62-b3f0-4b6a-83aa-76e5dc3548ac
md"""
The funding agency is impressed at the resolution you achieved, and your instrument is saved! You now go on to win the Nobel Prize due to discoveries the instrument makes. And it was all because you used `Quantity` as the input of code you shared.
"""

# ╔═╡ 7d4f9a02-81a7-4bdc-a22d-3435192f9f15
let
	sol = details("Example solution",
	md"""
	```julia
	v_orb(M, r) = sqrt(u"Constants.G" * M / r)
	```
	""")
	
	md"""
	!!! tip "Exercises"
		Write a function that computes the Keplerian velocity you worked out in section 1 (using `Quantity` input and outputs, of course), but allowing for an arbitrary mass and orbital radius. Try it with some reasonable numbers for satellites orbiting the Earth, a moon of Jupiter, or an extrasolar planet. Feel free to use wikipedia or similar for the masses and distances.
	
		$(sol)
	"""
end

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
# ╟─9a002059-650b-4dc1-9395-6f266ed35500
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
# ╟─4ee79c7c-5c4e-457b-b713-103937e50355
# ╠═f2b222ec-0783-487e-9c52-835976a555b6
# ╟─e4232b15-3369-438e-994b-042aab477a7f
# ╟─c749ce2d-17ae-45f4-b721-3f486b1cbc23
# ╟─3a3328d5-31de-4deb-bff7-d25c1fcbc4ef
# ╟─fb0f2941-2e56-45db-91fe-e8d6744a00e0
# ╠═bf618161-ef96-4445-8fc8-25dc5f662242
# ╠═1a16d2e1-7898-415c-9a98-a1f2f87e08be
# ╟─e0a2c745-41cf-4359-bb99-3117fbb507cc
# ╠═d45bad3e-4c82-436b-8e64-c61e5cf65c2f
# ╠═ee46b865-dbd4-4939-b241-514941dd138d
# ╠═b183c219-a39b-427a-b740-5675fcc175ca
# ╟─323d1169-c59f-456c-9b19-093e26f84214
# ╟─d0215082-52df-4eae-a409-7c3c7cfc69ed
# ╠═3a24f4aa-b074-4beb-b77d-6778e2fe580a
# ╠═1e203abc-4a02-4d01-bf9e-9df636e70ebd
# ╠═98869c25-2644-47d1-b8cc-05699292f2a8
# ╠═18985bb2-6cce-4417-97a2-c2da7b9e0428
# ╟─d9e2b828-9c32-4638-bfcd-0c16b221aa43
# ╠═4092a893-818e-49f4-93d0-be7bd652dddc
# ╠═c5c907d0-a669-43f3-8030-0bd67452f0a1
# ╟─a0b64ba9-cbec-404e-bf39-aee03ae407ae
# ╠═82f11626-722b-46fa-9173-8b5d7a80a190
# ╟─6aff587f-be6f-4fd8-96df-9c12f3769f32
# ╠═ef498536-8fb8-46e3-9d8d-f7eb3994a3f5
# ╟─25598fcd-761d-40b8-95ac-8b68a00026da
# ╟─eba2f06f-ebe9-492d-81d2-1cc4fccd5b0a
# ╠═de22c778-1529-4177-85e6-a0178f437a8c
# ╠═c262aa56-b9b6-4da4-8810-d6120f0724c6
# ╟─1e8ec5c5-3d41-44a7-8dfa-81d981750d9e
# ╟─c6b35992-973b-49d4-bfa7-855e8fe10924
# ╠═85f84c90-46e6-4d3f-bba2-6fc190a05533
# ╠═16cdd7e8-0ec2-4953-ad8b-f53669d05eda
# ╠═522d88d6-5f60-401b-8786-0236c0859eda
# ╟─fbd0d75b-b3af-4eb5-a249-423317d656e5
# ╠═0eefe540-1669-4510-bef7-8fb7f8be68f1
# ╠═4f59f4e6-c735-4a40-b171-744b17eb37fd
# ╠═b057812d-1a8e-479c-8192-1c1acdb97d23
# ╠═586b1993-be11-4437-9d3b-789511ea26bf
# ╟─e81f3b44-bc04-4768-81c2-c2733958f5a9
# ╠═88966b7e-b261-45a6-9f12-6b74f24c03e4
# ╠═ff7def67-eb96-45b3-959a-2b01c7804e07
# ╟─965a8e44-e96e-434e-8454-cced531ae9d2
# ╟─802f5cad-a0c0-4426-94f2-426f89dea7e1
# ╟─2eed0bc8-2306-42bb-9803-33ae131021e6
# ╠═9922bdbe-cf9e-486a-8e8a-28a7ded12d04
# ╟─118aff27-dc6d-4686-b64a-ff840b731030
# ╠═ffc4bdc3-16e3-4b3e-b448-cfdc16428378
# ╟─71cca345-cffc-4aee-a671-08901a92babe
# ╠═1c3faf9f-adcf-4232-acf4-655d828623e3
# ╟─7cfb17fc-81a4-4522-a88b-2fd10001baaa
# ╠═ce1804da-df61-4e4c-aa3b-d950e83b7c13
# ╟─5590c909-6103-4427-b1a2-4cf35265dc07
# ╠═f7486c08-b74b-471f-bbc2-690439a487f8
# ╟─c143ab62-b3f0-4b6a-83aa-76e5dc3548ac
# ╟─7d4f9a02-81a7-4bdc-a22d-3435192f9f15
# ╟─59b4d441-9a74-468f-ad8c-882516a09049
# ╠═17c6b7df-a8b3-45d1-9491-526afce11318
# ╠═bedc8ccd-e6f6-4dd1-a0b6-1889f4b5b658

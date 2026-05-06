### A Pluto.jl notebook ###
# v0.20.24

#> [frontmatter]
#> title = "Using Astropy Quantities and Units for astrophysical calculations"
#> layout = "layout.jlhtml"
#> date = "2025-12-31"
#> description = "Work with units in astrophysical calculations using Unitful.jl."
#> tags = ["units", "Unitful", "radio astronomy", "data cubes", "Makie", "plotting"]

using Markdown
using InteractiveUtils

# ╔═╡ 190d57d0-a588-48b5-9345-ada852ffb6b0
using Printf

# ╔═╡ f4b328d3-fa70-472f-9f38-fd740fe455f2
begin
	using Random
	using CairoMakie
	using Unitful, UnitfulAstro, UnitfulEquivalences
	using PhysicalConstants.CODATA2018: G, h, k_B
	using Distributions
end

# ╔═╡ 98b45be6-2e62-4aad-92e5-ef5a33f5af84
using Unitful: °

# ╔═╡ 713d2246-e45d-478f-9e4b-eceb25f31467
using PlutoUI

# ╔═╡ 053db95a-6175-11ef-2347-21d349323cd9
md"""
# Using Unitful quantities for astrophysical calculations
"""

# ╔═╡ c21af58f-45c8-4221-8832-ec9e9b2aa652
md"""
This notebook is modified from <https://learn.astropy.org/tutorials/quantities.html>

Original authors: Ana Bonaca, Erik Tollerud, Jonathan Foster, Lia Corrales, Kris Stern, Stephanie T. Douglas
"""

# ╔═╡ c79f86dd-59d7-4a23-bf9e-534cdb68d4bf
md"""
!!! tip "Learning Goals"
	
	- Use `Quantity` objects to estimate a hypothetical galaxy's mass
	- Take advantage of constants in the PhysicalConstants.jl package
	- Print formatted unit strings
	- Plot `Quantity` objects with unit labels
	- Do math with `Quantity` objects
	- Convert quantities with `astropy.units`
	- Convert between wavelength and energy with `Spectral` equivalencies
	- Use the small angle approximation with `astropy.units.dimensionless_angles` equivalencies
	- Write functions that take `Quantity` objects instead of plain floats
	- Make synthetic radio observations
	- Use `Quantity` objects such as data cubes to facilitate a full derivation of the total mass of a molecular cloud
"""

# ╔═╡ 89a6629b-052f-4b35-8fe8-91e5cfc30115
md"""
!!! note "Keywords"
	units, radio astronomy, data cubes, Makie, plotting
"""

# ╔═╡ 77233996-43ba-40c6-a582-f7e0d4681d5d
md"""
!!! warning "Summary"
	In this tutorial we present some examples showing how Astropy's `Quantity` object can make astrophysics calculations easier. The examples include calculating the mass of a galaxy from its velocity dispersion and determining masses of molecular clouds from CO intensity maps. We end with an example of good practices for using quantities in functions you might distribute to other people.
"""

# ╔═╡ 18c6c7a5-a7d9-42ff-b7e0-330cbaf84b5b
md"""
## Companion Content
"""

# ╔═╡ 6a2f1337-6535-49bb-9da2-8cfca623085b
md"""
[Tools for Radio Astronomy](https://www.springer.com/gp/book/9783662053942) by Rohlfs & Wilson
"""

# ╔═╡ 76ed032d-6c8d-419d-8cf2-ecab5eab9652
md"""
For an in-depth discussion of `Quantity` objects, see the astropy documentation section.
"""

# ╔═╡ 3f0d31af-d6f1-44e3-bf93-9205f54e29ab
md"""
## Imports
"""

# ╔═╡ 022039cc-5b17-4a33-9e43-639bd79e1365
md"""
## Preliminaries
"""

# ╔═╡ cb2c8246-e576-4476-8640-6c7602245306
# The `Random.seed!` function makes the "random" number generator always give the same numbers wherever you run it, so it's reproducible randomness.
rng = Random.seed!(0)

# ╔═╡ 639c4861-cb16-417f-8de3-9aa3c27d0db9
md"""
The Unitful.jl package gives us a handy macro `u"..."` demonstrated below. This will make working with `Quantity` objects much easier.
"""

# ╔═╡ 7e302677-45cd-401d-945c-dcdea7349bb1
md"""
We can also use the PhysicalConstants.jl module where typical physical constants are available. The constants are stored as objects of a subtype of `Quantity`, so they behave just like a `Quantity`. Here, we'll only need the gravitational constant ``G``, Planck's constant ``h``, and Boltzmann's constant, ``k_B``.
"""

# ╔═╡ 4ece2ccf-5113-4eb0-bb12-b1574738d649
md"""
## 1. Galaxy mass
"""

# ╔═╡ 8bb28f7a-a625-4fec-8677-0f0b41c1c432
md"""
In this first example, we will use `Quantity` objects to estimate a hypothetical galaxy's mass, given its half-light radius and radial velocities of stars in the galaxy.
"""

# ╔═╡ fc73c5b7-8fde-42f2-9ffd-d99f75dab037
md"""
Let's assume that we measured the half-light radius of the galaxy to be 29 pc projected on the sky at the distance of the galaxy. This radius is often called the "effective radius", so we'll store it as a `Quantity` object with the name `Reff`. The easiest way to create a `Quantity` object is by multiplying the value with its unit. Units are accessed as `u"unit"`, in this case `u"pc"`.
"""

# ╔═╡ b291375e-c68a-4c92-a8ef-2528cf9be177
Reff = 29u"pc"

# ╔═╡ 30c66d47-c82c-487e-8e23-3786fbfe2f0e
md"""
We can access the value and unit of a `Quantity` using the `ustrip()` and `unit()` functions.
"""

# ╔═╡ dc3bdf39-c7b4-402d-851a-d1971d79d7bc
Print("""Half light radius
Value: $(ustrip(Reff))
unit: $(unit(Reff))
""")

# ╔═╡ 6eb5071d-327c-4ebb-8167-f682c6900faf
md"""
Furthermore, we can convert the radius in parsecs to any other unit of length using the `uconvert()` function. Here, we convert it to metres.
"""

# ╔═╡ 91eb8613-7200-4176-9fb3-309450a9a052
uconvert(u"m", Reff)

# ╔═╡ 9bacca3a-ea0a-4de2-a078-39ea6767639a
md"""
Units are also "callable," meaning we can treat them like functions. Even though the first example with calling a unit is a bit awkward, the `|>` syntax makes it much more natural.
"""

# ╔═╡ beb01856-9d91-4fd8-8253-f73b87ee47d1
u"m"(Reff)

# ╔═╡ 39fa4bb9-312a-491d-bc22-1c8a4488cb47
Reff |> u"m"

# ╔═╡ 38f72bb4-8fac-456d-825d-205a47f46042
md"""
Next, we'll first create a synthetic dataset of radial velocity measurements, assuming a normal distribution with a mean velocity of 206 km/s and a velocity dispersion of 4.3 km/s. We can do this by creating a `Normal` distribution from the Distributions.jl package.
"""

# ╔═╡ 7173b8c2-0599-4170-a5f0-c2258dcee259
v̄, σin = 206, 4.3;

# ╔═╡ e39d8e4c-6455-4961-91f1-0abd15559afe
v = rand(Normal(v̄, σin), 500)u"km/s"

# ╔═╡ 35f30d99-b77f-4495-9564-627be22609ce
md"""
Now we can plot a histogram of the velocity dataset.
"""

# ╔═╡ 14c6db7d-c635-4082-9bf9-e190f343daa5
md"""
!!! warning
    Makie doesn't do histograms with units, this part may need to be rewritten
"""

# ╔═╡ 9e92d1d8-fd03-4dc5-b41b-a24afb3f6bd0
hist(v .|> ustrip)

# ╔═╡ 54e2c96c-bc2f-4135-83df-2d4975842441
md"""
Now we can calculate the velocity dispersion of the galaxy. This demonstrates how you can perform basic operations like subtraction and division with `Quantity` objects, and also use them in standard Julia functions such as `mean()` and `size()`. They retain their units through these operations just as you would expect them to.
"""

# ╔═╡ 96270bd1-2072-4de2-9752-d2d00a5035ed
σ = √(sum((v .- mean(v)).^2) / length(v))

# ╔═╡ f7851441-d76f-4536-b057-8f5eba8be595
md"""
Now for the actual mass calculation. If a galaxy is pressure-supported (for example, an elliptical or dwarf spheroidal galaxy), its mass within the stellar extent can be estimated using a straightforward formula: ``M_{1/2} = 4σ^2 R_\text{eff}/G``. There are caveats to the use of this formula for science -- see Wolf et al. 2010 for details. For demonstrating `Quantity`, you can accept that this is often good enough. For the calculation, we can multiply the quantities together, and Unitful will keep track of the units.
"""

# ╔═╡ d6cdbf2a-3509-4fbb-b9ea-66116d613e70
M = 4 * σ^2 * Reff / G

# ╔═╡ 4b99f953-41de-45c6-b1a3-e136f1162615
M |> u"kg"

# ╔═╡ 07bb88c7-3fe2-42da-8bd7-20c2d8387a82
md"""
We can also easily express the mass in whatever form you like -- solar masses are common in astronomy, or maybe you want the default SI and CGS units.
"""

# ╔═╡ b9651667-9b63-4e43-84b6-6b7bf183cace
Print("""Galaxy mass
in solar units: $(M|>u"Msun")
SI units: $(M|>u"kg")
CGS units: $(M|>u"g")""")

# ╔═╡ 393cf7ca-df52-4b2c-a663-dfc6edcff89f
md"""
Or, if you want the log of the mass, you can just use the builtin `log10` as long as the logarithm's argument is dimensionless.
"""

# ╔═╡ 9860dfc9-ed59-4c87-b8db-2c05678acf57
log10(M/u"Msun")

# ╔═╡ 0d37b7f1-ff6f-434e-8a47-e0c65073e479
md"""
However, you can't take the log of something with units, as that is not mathematically sensible.
"""

# ╔═╡ 6b7a3aa8-3ac2-42c7-b56e-46ad8e0d4061
log10(M)

# ╔═╡ 754f4db1-516b-4063-b4c7-0a4bd268b935
md"""
## Exercises
"""

# ╔═╡ cbfa1f77-15b8-4d7f-bb08-987078745789
md"""
Use Kepler's law in the form given below to determine the (circular) orbital speed of the Earth around the Sun in km/s. No need to look up constants or conversion factors to do this calculation -- it's all in Unitful.jl and PhysicalConstants.jl.
```math
v = \sqrt{\frac{G M_⊙}{r}}
```
"""

# ╔═╡ f8ada759-3ed7-4e6a-9cb8-7e63d59d81a7


# ╔═╡ 5712ed02-49b7-4eb8-a1eb-21a6e25726e0
md"""
There's a much easier way to figure out the velocity of the Earth using just two units or quantities. Do that and then compare to the Kepler's law answer (the easiest way is probably to compute the percentage difference, if any).
"""

# ╔═╡ 59ec195b-00cd-4d0c-8e4f-c9bf11f92e1c


# ╔═╡ 043c8bdd-f486-41af-9e56-ffcb24872d9f
md"""
(Completely optional, but a good way to convince yourself of the value of `Quantity`:) Do the above calculations by hand -- you can use a calculator (or Julia just for its arithmatic) but look up all the appropriate conversion factors and use paper-and-pencil approaches for keeping track of them all. Which one took longer?
"""

# ╔═╡ 161ad10d-5c05-42aa-a64b-8343296b5617
md"""
## 2. Molecular cloud mass
"""

# ╔═╡ 981bada9-a648-4665-88d4-bb1fc9a448b5
md"""
In this second example, we will demonstrate how using Quantity objects can facilitate a full derivation of the total mass of a molecular cloud using radio observations of isotopes of Carbon Monoxide (CO).
"""

# ╔═╡ ed7d71f8-6b09-426c-9bb1-9a5243b77b32
md"""
### Setting up the data cube
"""

# ╔═╡ 8a71946c-9b3f-4b27-a8f3-5092f0e306bd
md"""
Let's assume that we've mapped the inner part of a molecular cloud in the J=1-0 rotational transition of C¹⁸O and are interested in measuring its total mass. The measurement produced a data cube with RA and Dec as spatial coordiates and velocity as the third axis. Each voxel in this data cube represents the brightness temperature of the emission at that position and velocity. Furthermore, we'll assume that we have an independent measurement of distance to the cloud _d_ = 250 pc and that the excitation temperature is known and constant throughout the cloud: _T_ₑₓ = 25 K.
"""

# ╔═╡ c28fc917-3cff-413d-a0bd-0970039b8bac
d = 250u"pc"

# ╔═╡ dedee06e-134b-4ec5-939d-5d16c87182c9
Tₑₓ = 25u"K"

# ╔═╡ 55d5a031-e21e-419e-ace7-0747f1e627d3
md"""
We'll generate a synthetic dataset, assuming the cloud follows a Gaussian distribution in each of RA, Dec and velocity. We start by creating a 100×100×300 array, such that the first coordinate is right ascension, the second is declination, and the third is velocity. We use an array comprehension to create data cubes for each of the three coordinates, and then use them in the formula for a Gaussian to generate an array with the synthetic data cube. In this cube, the cloud is positioned at the center of the cube, with _σ_ and the center in each dimension shown below. Note in particular that the _σ_ for RA and Dec have different units from the center, but astropy automatically does the relevant conversions before computing the exponential.
"""

# ╔═╡ f0c9baca-86e5-46a5-a1c7-f88335154f29
# Cloud's center
cen_ra, cen_dec, cen_v = 52.25°, 0.25°, 15u"km/s";

# ╔═╡ c297447b-e9ca-4fa6-bff6-f36050280b8e
# Cloud's size
sig_ra, sig_dec, sig_v = 3u"arcminute", 4u"arcminute", 3u"km/s";

# ╔═╡ 9c66da40-6184-4260-8c03-f6807efb4aed
begin
    #1D coordinate quantities
    ra_cube = range(52°, 52.5°, length = 100)
    dec_cube = range(0°, 0.5°, length = 100)
    v_cube = range(0u"km/s", 30u"km/s", length = 300)
end;

# ╔═╡ 237fade9-8273-493e-8cd0-abaed6f32985
data_gauss = [
    exp(
        -0.5 * (
            ((ra - cen_ra) / sig_ra)^2
            + ((dec - cen_dec) / sig_dec)^2
            + ((v - cen_v) / sig_v)^2
        )
    )
    for ra in ra_cube, dec in dec_cube, v in v_cube
]

# ╔═╡ b616ca13-d26b-4927-aab2-9f7f01e42155
md"""
The units of the exponential are dimensionless, so we multiply the data cube by K to get brightness temperature units. Radio astronomers use a rather odd set of units [K km/s] as of integrated intensity (that is, summing all the emission from a line over velocity). As an aside for experts, we're setting up our artificial cube on the main-beam temperature scale (``T_\text{MB}``) which is the closest we can normally get to the actual brightness temperature of our source.
"""

# ╔═╡ a519bf53-decd-416a-90f9-e25b64236e29
data = data_gauss*u"K";

# ╔═╡ 7f541f69-80e1-4289-bdfa-efd04401b769
md"""
We will also need to know the width of each velocity bin and the size of each pixel, so let's calculate that now.
"""

# ╔═╡ 5b7c780f-be18-4627-a254-3a7e890fdf9f
# Average pixel size in right ascension
# This is only right if dec ~ 0, because of the cos(dec)factor.
Δra = (maximum(ra_cube) - minimum(ra_cube)) / length(ra_cube) |> u"arcsecond"

# ╔═╡ e86e25e6-9a75-427c-ab03-565b37cee91e
# Average pixel size in declination
Δdec = (maximum(dec_cube) - minimum(dec_cube)) / length(dec_cube) |> u"arcsecond"

# ╔═╡ e13b9288-60d4-4be1-ab5c-9a521a14a9f3
# Average velocity bin width
Δv = (maximum(v_cube) - minimum(v_cube)) / length(v_cube)

# ╔═╡ 7e834830-f098-4629-9968-0025b9804577
md"""
We're interested in the integrated intensity over all of the velocity channels, so let's create a 2D quantity array by summing our data cube along the velocity axis (multiplying by the velocity width of a pixel).
"""

# ╔═╡ e5e8a468-2292-495b-a266-3439de08c0ed
intcloud = dropdims(sum(data .* Δv, dims = 3), dims = 3)

# ╔═╡ b7e084f1-9d4d-4455-88f8-ff3bb3cc6885
md"""
We can plot the 2D quantity using Makie's `heatmap` function, by passing the quantity's value.
"""
# Similarly, we can set the correct extent using the values of ``x_i`` and ``x_f``.
# Finally, we can set the colorbar label to have proper units.

# ╔═╡ 139cb59e-927b-45d1-892c-dd300eb71b12
md"""
!!! warning
    Lack of feature parity: the colorbar intensity does not seem to show units (issues: [Makie.jl#3893](https://github.com/MakieOrg/Makie.jl/issues/3893) and [Makie.jl#3897](https://github.com/MakieOrg/Makie.jl/issues/3897))
"""

# ╔═╡ 3f75b0cc-1b58-4172-ba07-70fd86767fe8
let
    fig, ax, hm = heatmap(
        ra_cube, dec_cube, intcloud,
        axis = (xlabel = "RA", ylabel = "Dec", xreversed = true),
    )
    Colorbar(fig[1,2], hm, label = "Intensity")
    fig
end

# ╔═╡ dbd270a4-dc73-4cdf-be5e-1ad6c2a98765
md"""
### Measuring The Column Density of CO
"""

# ╔═╡ bd33224b-96c4-4e0c-a93a-f454bfe1be0b
md"""
In order to calculate the mass of the molecular cloud, we need to measure its column density. A number of assumptions are required for the following calculation; the most important are that the emission is optically thin (typically true for C¹⁸O) and that conditions of local thermodynamic equilibrium hold along the line of sight. In the case where the temperature is large compared to the separation in energy levels for a molecule and the source fills the main beam of the telescope, the total column density for C¹³O is

```math
N = \frac{C}{1 - e^{-B}} \int T_B(V) \, dV
```

where the constants _C_ and _B_ are given by
```math
\begin{align*}
C &= 3.0 × {10}^{14} \; \mathrm{K^{-1} {cm}^{-2} {km}^{-1} s} \; {\left(\frac{ν}{ν_{13}}\right)}^2 \frac{A_{13}}{A} \\
B &= \frac{hν}{k_\text{B} T}
\end{align*}
```

(Rohlfs & Wilson [Tools for Radio Astronomy](https://www.springer.com/gp/book/9783662053942)).
"""

# ╔═╡ bafa7022-05c9-4907-ac52-24539fb54a24
md"""
Here we have given an expression for _C_ scaled to the values for C¹³O (_ν_₁₃ and _A_₁₃). In order to use this relation for C¹⁸O, we need to rescale the frequencies _ν_ and Einstein coefficients _A_. _C_ is in funny mixed units, but that's okay. We'll define it as a Quantities object and not have to worry about it.
"""

# ╔═╡ a82eea6b-f888-485b-99f1-6b365a6ecb94
md"""
First, we look up the wavelength for these emission lines and store them as quantities.
"""

# ╔═╡ c6d94073-e72a-4981-b8fa-535e44b18b91
λ₁₃ = 2.60076u"μm"; # variable name typed as \lambda<TAB> \_1<TAB> \_3<TAB>

# ╔═╡ e53f99bf-dff7-41c8-801a-13cce8c41178
λ₁₈ = 2.73079u"μm";

# ╔═╡ 2656341a-bfdb-4cd3-b9dd-504e27ba067f
md"""
Since the wavelength and frequency of light are related using the speed of light, we can convert between them. However, doing so just using the `uconvert()` function or its equivalent `|>` fails, as units of length and frequency are not convertible:
"""

# ╔═╡ e7896cb2-f47f-424d-a05d-e76cc8c3bc2c
λ₁₃ |> u"Hz"

# ╔═╡ 8f84afe1-bbf6-419f-b386-d2219e423d21
md"""
Fortunately, the Unitful ecosystem comes to the rescue by providing a feature called "unit equivalences." Equivalences provide a way to convert between two physically different units that are not normally equivalent, but in a certain context have a one-to-one mapping. For more on equivalencies, see [the documentation for UnitfulEquivalences.jl](https://sostock.github.io/UnitfulEquivalences.jl).
"""

# ╔═╡ 5b7fa766-f66e-44e7-89e1-dc9f60e0b78a
md"""
In this case, passing the `Spectral()` argument to `uconvert()` provides the equivalences necessary to handle conversions between wavelength and frequency:
"""

# ╔═╡ 34ace9db-a63e-45b3-b0c0-18e3fa638d2a
ν₁₃ = uconvert(u"Hz", λ₁₃, Spectral())

# ╔═╡ 75dc88c0-d1f1-43e4-881e-3c9c69810f7b
ν₁₈ = uconvert(u"Hz", λ₁₈, Spectral())

# ╔═╡ 6c4a7e27-acb8-4962-b723-ef40658b09d2
md"""
Next, we look up Einstein coefficients (in units of s⁻¹), and calculate the ratios in constant _C_. Note how the ratios of frequency and Einstein coefficient units are dimensionless, so the unit of _C_ is unchanged.
"""

# ╔═╡ dfc3ff19-c8cf-4d89-9cda-184cf9837cc1
A₁₃ = 7.4e-8/u"s";

# ╔═╡ a976247d-59e7-4cb5-8782-361dfb866397
A₁₈ =  8.8e-8/u"s";

# ╔═╡ 1272fd51-d070-4ff7-9dac-3113d5be275b
C = 3e14u"s/(K*cm^2*km)" * (ν₁₈/ν₁₃)^3 * (A₁₃/A₁₈)

# ╔═╡ a91cd537-2d8b-4e17-a7f3-edc297510a2d
md"""
Now we move on to calculate the constant _B_. This is given by the ratio of ``\dfrac{hν}{k_\text{B}T}``, where _h_ is Planck's constant, ``k_\text{B}`` is the Boltzmann's constant, _ν_ is the emission frequency, and _T_ is the excitation temperature. The constants were imported from PhysicalConstants.jl, and the other two values are already calculated, so here we just take the ratio.
"""

# ╔═╡ 00c4ec7f-d795-43b5-94e2-c5e50c6e0cd5
B = h * ν₁₈ / (k_B * Tₑₓ)

# ╔═╡ cfc2a2ae-17c9-4b9f-8a9e-256607e02dac
md"""
The units of _B_ are Hz sec, which can be converted to a raw float if you actually care about its value. Usually this is not necessary, though. Quantities are at their best if you use them without worrying about intermediate units, and only convert at the very end when you want a final answer.
"""

# ╔═╡ 1e63ee80-b9a0-4061-aafa-d3cfd1175c50
B |> NoUnits

# ╔═╡ 85431572-f106-4764-acce-2f8aa4a0bb6f
md"""
At this point we have all the ingredients to calculate the number density of CO molecules in this cloud. We already integrated (summed) over the velocity channels above to show the integrated intensity map, but we'll do it again here for clarity. This gives us the column density of CO for each spatial pixel in our map. We can then print out the peak column column density.
"""

# ╔═╡ 12ae26eb-afa5-44d9-b27d-e7e0d13b6932
NCO = let
    integrated = dropdims(sum(data * Δv, dims = 3), dims = 3)
    C * integrated / (1 - exp(-B))
end

# ╔═╡ d9e46330-597b-41fb-9b21-765145f864ec
Print("Peak CO Column density: $(maximum(NCO))")

# ╔═╡ 7deb8ee0-dca9-41a1-8531-cdbef591f83a
md"""
### CO to Total Mass
"""

# ╔═╡ 459751c9-37de-4a8e-8866-fb20f0b60ab3
md"""
We are using CO as a tracer for the much more numerous H₂, the quantity we are actually trying to infer. Since most of the mass is in H₂, we calculate its column density by multiplying the CO column density with the (known/assumed) H₂/CO ratio.
"""

# ╔═╡ 9833719a-f1d1-4c65-b92f-024c300ec01c
H₂_CO_ratio = 5.9e6;

# ╔═╡ d263fbba-f7cf-45af-90e0-ffcd0c4edecb
NH₂ = NCO * H₂_CO_ratio

# ╔═╡ fb547293-3860-455e-bdd8-5bd84d39521b
Print("Peak H₂ column density: $(maximum(NH₂))")

# ╔═╡ c9e0945b-dade-4f7e-b364-62e8ac180c28
md"""
That's a peak column density of roughly 50 magnitudes of visual extinction (assuming the conversion between ``N_\mathrm{H_2}`` and ``A_V`` from Bohlin et al. 1978), which seems reasonable for a molecular cloud.
"""

# ╔═╡ 2f81cb18-0199-47f5-9333-b89e72a8af6b
md"""
We obtain the mass column density by multiplying the number column density by the mass of an individual H₂ molecule.
"""

# ╔═╡ 140545d5-cad4-4729-8233-4c959bd3eaae
mH₂ = 2*1.008u"u"  # aka atomic mass unit/amu/Dalton

# ╔═╡ bc0d09b3-bd4c-47ba-9e99-10e2aaaaa39c
ρ = NH₂ * mH₂

# ╔═╡ 605d4fce-4f52-479f-bd02-273641215330
md"""
A final step in going from the column density to mass is summing up over the area. If we do this in the straightforward way of length × width of a pixel, this area is then in units of deg².
"""

# ╔═╡ c97465fc-c081-4756-bb72-4bdfbac48924
Δap = Δra * Δdec |> °^2

# ╔═╡ db43ab8e-263d-4104-bd3a-cc32ab995098
Δa = Δap * d^2 |> u"pc^2"

# ╔═╡ 67fd0710-4e7e-463f-87ff-9b69231111a7
Δa |> u"cm^2"

# ╔═╡ 698088de-a916-4b04-bb1d-d23e08df55ee
Δap |> u"sr"

# ╔═╡ ba72c6e2-48af-4a4c-a316-86d1cbf3a0e6
md"""
TODO: astropy has a problem here with treating angles dimensionally, but Unitful does it fine? investigate and write up proper instructions
"""

# ╔═╡ 14079708-26d1-490d-8038-4839e7dacc18
(Δap |> u"sr") * d^2 |> u"cm^2"

# ╔═╡ adc04cc3-c994-42b2-8ecd-fd50c047d811
cloud_mass = sum(ρ * Δa)

# ╔═╡ 406de0ba-fe94-4847-a198-3988e49a7e0a
cloud_mass |> u"Msun"

# ╔═╡ d8711a00-21a0-4662-a7be-8aacfe7a5087
md"""
## Exercises
"""

# ╔═╡ c43da11f-3d69-4661-b090-36ff5aefa541
md"""
The astro material was pretty heavy on that one, so let's focus on some associated statistics using Quantity's array capabililities. Compute the median and mean of the data with the np.mean and np.median functions. Why are their values so different?
"""

# ╔═╡ f1d9ca5c-a1e2-4270-b3dd-7386af419d72


# ╔═╡ 495be19b-2730-4ad7-b653-1aad5e87d072
md"""
Similarly, compute the standard deviation and variance (if you don't know the relevant functions, look it up in the numpy docs or just type np. and a code cell). Do they have the units you expect?
"""

# ╔═╡ 3ee76460-77b6-4703-8f67-65c708162073


# ╔═╡ a919ff70-8b4d-4cca-b695-647607e25857
md"""
# Notebook setup 🔧
"""

# ╔═╡ bbdd7d17-37cb-4d37-811a-5a6320d66f17
TableOfContents(; title = "On this page", depth = 4)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
PhysicalConstants = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
UnitfulEquivalences = "da9c4bc3-91c8-4f02-8a40-6b990d2a7e0c"

[compat]
CairoMakie = "~0.15.9"
Distributions = "~0.25.123"
PhysicalConstants = "~0.2.4"
PlutoUI = "~0.7.80"
Unitful = "~1.28.0"
UnitfulAstro = "~1.2.2"
UnitfulEquivalences = "~0.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "b90bc19d65f9fcf6fee1f26867b5bd649c8ff398"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "2eeb2c9bef11013efc6f8f97f32ee59b146b09fb"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.44"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "35ea197a51ce46fcd01c4a44befce0578a1aaeca"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.5.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "fa072933899aae6dc61dde934febed8254e66c6a"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.9"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d0efe2c6fdcdaa1c161d206aa8b933788397ec71"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.6+0"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSolve]]
git-tree-sha1 = "78ea4ddbcf9c241827e7035c3a03e2e456711470"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.6"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "3b4be73db165146d8a88e47924f464e55ab053cd"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.7"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e86f4a2805f7f19bec5129bc9150c38208e5dc23"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.4"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "c55f5a9fd67bdbc8e089b5a3111fe4292986a8e8"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.6"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "66381d7059b5f3f6162f28831854008040a4e905"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+1"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6522cfb3b8fe97bec632252263057996cbd3de20"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.18.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "1d86f920d19a2eb2b9528d62a834412f4822eaa1"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.5"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "67c6f1f085cb2671c93fe34244c9cccde30f7a26"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.5.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "4260cfc991b8885bf747801fb60dd4503250e478"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.11"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "68af66ec16af8b152309310251ecb4fbfe39869f"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.9"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "7eb8cdaa6f0e8081616367c10b31b9d9b34bb02a"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.7"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf"]
git-tree-sha1 = "cb47f69a1cab9dcec7ff4a5d6e163410d6905866"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.14.1"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsMakieExt = "Makie"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f2b3b9e52a5eb6a3434c8cca67ad2dde011194f4"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.30+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "135492b7e97fc86d9b132b96a54d2d3dd3e0c6a8"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.8+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.PhysicalConstants]]
deps = ["Measurements", "Roots", "Unitful"]
git-tree-sha1 = "9be04ab1a9e54508348a2aacb5b5a4f04c85d397"
uuid = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
version = "0.2.4"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "fbc875044d82c113a9dee6fc14e16cf01fd48872"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.80"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "b2f70f34eb9973572d55c332933c6a04c911f549"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.14"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2700b235561b0335d5bef7097a111dc513b8655e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.2"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "ad8002667372439f2e3611cfd14097e03fa4bccd"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.3"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "fa95b3b097bcef5845c142ea2e085f1b2591e92c"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.7.1"
weakdeps = ["Measurements", "StaticArraysCore", "Tables"]

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "08c10bc34f4e7743f530793d0985bf3c254e193d"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.8"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "57e1b2c9de4bd6f40ecb9de4ac1797b81970d008"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.28.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "79875b1f2e4bf918f0702a5980816955066d9ae2"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.7.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "fbe44a0ade62ae5ed0240ad314dfdd5482b90b40"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.2"

[[deps.UnitfulEquivalences]]
deps = ["Unitful"]
git-tree-sha1 = "76fc2f7fdc87531a1018eb7d647df7c29daf36b7"
uuid = "da9c4bc3-91c8-4f02-8a40-6b990d2a7e0c"
version = "0.2.0"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e2a7072fc0cdd7949528c1455a3e5da4122e1153"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.56+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"
"""

# ╔═╡ Cell order:
# ╟─053db95a-6175-11ef-2347-21d349323cd9
# ╟─c21af58f-45c8-4221-8832-ec9e9b2aa652
# ╟─c79f86dd-59d7-4a23-bf9e-534cdb68d4bf
# ╟─89a6629b-052f-4b35-8fe8-91e5cfc30115
# ╟─77233996-43ba-40c6-a582-f7e0d4681d5d
# ╟─18c6c7a5-a7d9-42ff-b7e0-330cbaf84b5b
# ╟─6a2f1337-6535-49bb-9da2-8cfca623085b
# ╟─76ed032d-6c8d-419d-8cf2-ecab5eab9652
# ╠═190d57d0-a588-48b5-9345-ada852ffb6b0
# ╟─3f0d31af-d6f1-44e3-bf93-9205f54e29ab
# ╠═f4b328d3-fa70-472f-9f38-fd740fe455f2
# ╟─022039cc-5b17-4a33-9e43-639bd79e1365
# ╠═cb2c8246-e576-4476-8640-6c7602245306
# ╟─639c4861-cb16-417f-8de3-9aa3c27d0db9
# ╟─7e302677-45cd-401d-945c-dcdea7349bb1
# ╟─4ece2ccf-5113-4eb0-bb12-b1574738d649
# ╟─8bb28f7a-a625-4fec-8677-0f0b41c1c432
# ╟─fc73c5b7-8fde-42f2-9ffd-d99f75dab037
# ╠═b291375e-c68a-4c92-a8ef-2528cf9be177
# ╟─30c66d47-c82c-487e-8e23-3786fbfe2f0e
# ╠═dc3bdf39-c7b4-402d-851a-d1971d79d7bc
# ╟─6eb5071d-327c-4ebb-8167-f682c6900faf
# ╠═91eb8613-7200-4176-9fb3-309450a9a052
# ╟─9bacca3a-ea0a-4de2-a078-39ea6767639a
# ╠═beb01856-9d91-4fd8-8253-f73b87ee47d1
# ╠═39fa4bb9-312a-491d-bc22-1c8a4488cb47
# ╟─38f72bb4-8fac-456d-825d-205a47f46042
# ╠═7173b8c2-0599-4170-a5f0-c2258dcee259
# ╠═e39d8e4c-6455-4961-91f1-0abd15559afe
# ╟─35f30d99-b77f-4495-9564-627be22609ce
# ╟─14c6db7d-c635-4082-9bf9-e190f343daa5
# ╠═9e92d1d8-fd03-4dc5-b41b-a24afb3f6bd0
# ╟─54e2c96c-bc2f-4135-83df-2d4975842441
# ╠═96270bd1-2072-4de2-9752-d2d00a5035ed
# ╟─f7851441-d76f-4536-b057-8f5eba8be595
# ╠═d6cdbf2a-3509-4fbb-b9ea-66116d613e70
# ╠═4b99f953-41de-45c6-b1a3-e136f1162615
# ╟─07bb88c7-3fe2-42da-8bd7-20c2d8387a82
# ╠═b9651667-9b63-4e43-84b6-6b7bf183cace
# ╟─393cf7ca-df52-4b2c-a663-dfc6edcff89f
# ╠═9860dfc9-ed59-4c87-b8db-2c05678acf57
# ╟─0d37b7f1-ff6f-434e-8a47-e0c65073e479
# ╠═6b7a3aa8-3ac2-42c7-b56e-46ad8e0d4061
# ╟─754f4db1-516b-4063-b4c7-0a4bd268b935
# ╟─cbfa1f77-15b8-4d7f-bb08-987078745789
# ╠═f8ada759-3ed7-4e6a-9cb8-7e63d59d81a7
# ╟─5712ed02-49b7-4eb8-a1eb-21a6e25726e0
# ╠═59ec195b-00cd-4d0c-8e4f-c9bf11f92e1c
# ╟─043c8bdd-f486-41af-9e56-ffcb24872d9f
# ╟─161ad10d-5c05-42aa-a64b-8343296b5617
# ╟─981bada9-a648-4665-88d4-bb1fc9a448b5
# ╟─ed7d71f8-6b09-426c-9bb1-9a5243b77b32
# ╟─8a71946c-9b3f-4b27-a8f3-5092f0e306bd
# ╠═c28fc917-3cff-413d-a0bd-0970039b8bac
# ╠═dedee06e-134b-4ec5-939d-5d16c87182c9
# ╟─55d5a031-e21e-419e-ace7-0747f1e627d3
# ╠═98b45be6-2e62-4aad-92e5-ef5a33f5af84
# ╠═f0c9baca-86e5-46a5-a1c7-f88335154f29
# ╠═c297447b-e9ca-4fa6-bff6-f36050280b8e
# ╠═9c66da40-6184-4260-8c03-f6807efb4aed
# ╠═237fade9-8273-493e-8cd0-abaed6f32985
# ╟─b616ca13-d26b-4927-aab2-9f7f01e42155
# ╠═a519bf53-decd-416a-90f9-e25b64236e29
# ╟─7f541f69-80e1-4289-bdfa-efd04401b769
# ╠═5b7c780f-be18-4627-a254-3a7e890fdf9f
# ╠═e86e25e6-9a75-427c-ab03-565b37cee91e
# ╠═e13b9288-60d4-4be1-ab5c-9a521a14a9f3
# ╟─7e834830-f098-4629-9968-0025b9804577
# ╠═e5e8a468-2292-495b-a266-3439de08c0ed
# ╟─b7e084f1-9d4d-4455-88f8-ff3bb3cc6885
# ╟─139cb59e-927b-45d1-892c-dd300eb71b12
# ╠═3f75b0cc-1b58-4172-ba07-70fd86767fe8
# ╟─dbd270a4-dc73-4cdf-be5e-1ad6c2a98765
# ╟─bd33224b-96c4-4e0c-a93a-f454bfe1be0b
# ╟─bafa7022-05c9-4907-ac52-24539fb54a24
# ╟─a82eea6b-f888-485b-99f1-6b365a6ecb94
# ╠═c6d94073-e72a-4981-b8fa-535e44b18b91
# ╠═e53f99bf-dff7-41c8-801a-13cce8c41178
# ╟─2656341a-bfdb-4cd3-b9dd-504e27ba067f
# ╠═e7896cb2-f47f-424d-a05d-e76cc8c3bc2c
# ╟─8f84afe1-bbf6-419f-b386-d2219e423d21
# ╟─5b7fa766-f66e-44e7-89e1-dc9f60e0b78a
# ╠═34ace9db-a63e-45b3-b0c0-18e3fa638d2a
# ╠═75dc88c0-d1f1-43e4-881e-3c9c69810f7b
# ╟─6c4a7e27-acb8-4962-b723-ef40658b09d2
# ╠═dfc3ff19-c8cf-4d89-9cda-184cf9837cc1
# ╠═a976247d-59e7-4cb5-8782-361dfb866397
# ╠═1272fd51-d070-4ff7-9dac-3113d5be275b
# ╟─a91cd537-2d8b-4e17-a7f3-edc297510a2d
# ╠═00c4ec7f-d795-43b5-94e2-c5e50c6e0cd5
# ╟─cfc2a2ae-17c9-4b9f-8a9e-256607e02dac
# ╠═1e63ee80-b9a0-4061-aafa-d3cfd1175c50
# ╟─85431572-f106-4764-acce-2f8aa4a0bb6f
# ╠═12ae26eb-afa5-44d9-b27d-e7e0d13b6932
# ╠═d9e46330-597b-41fb-9b21-765145f864ec
# ╟─7deb8ee0-dca9-41a1-8531-cdbef591f83a
# ╟─459751c9-37de-4a8e-8866-fb20f0b60ab3
# ╠═9833719a-f1d1-4c65-b92f-024c300ec01c
# ╠═d263fbba-f7cf-45af-90e0-ffcd0c4edecb
# ╠═fb547293-3860-455e-bdd8-5bd84d39521b
# ╟─c9e0945b-dade-4f7e-b364-62e8ac180c28
# ╟─2f81cb18-0199-47f5-9333-b89e72a8af6b
# ╠═140545d5-cad4-4729-8233-4c959bd3eaae
# ╠═bc0d09b3-bd4c-47ba-9e99-10e2aaaaa39c
# ╟─605d4fce-4f52-479f-bd02-273641215330
# ╠═c97465fc-c081-4756-bb72-4bdfbac48924
# ╠═db43ab8e-263d-4104-bd3a-cc32ab995098
# ╠═67fd0710-4e7e-463f-87ff-9b69231111a7
# ╠═698088de-a916-4b04-bb1d-d23e08df55ee
# ╟─ba72c6e2-48af-4a4c-a316-86d1cbf3a0e6
# ╠═14079708-26d1-490d-8038-4839e7dacc18
# ╠═adc04cc3-c994-42b2-8ecd-fd50c047d811
# ╠═406de0ba-fe94-4847-a198-3988e49a7e0a
# ╟─d8711a00-21a0-4662-a7be-8aacfe7a5087
# ╟─c43da11f-3d69-4661-b090-36ff5aefa541
# ╠═f1d9ca5c-a1e2-4270-b3dd-7386af419d72
# ╟─495be19b-2730-4ad7-b653-1aad5e87d072
# ╠═3ee76460-77b6-4703-8f67-65c708162073
# ╟─a919ff70-8b4d-4cca-b695-647607e25857
# ╠═bbdd7d17-37cb-4d37-811a-5a6320d66f17
# ╠═713d2246-e45d-478f-9e4b-eceb25f31467
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

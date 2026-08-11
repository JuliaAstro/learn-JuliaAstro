### A Pluto.jl notebook ###
# v1.0.2

#> [frontmatter]
#> title = "Units and integration"
#> layout = "layout.jlhtml"
#> date = "2025-12-31"
#> description = "TODO"
#> tags = ["modeling", "units", "synphot", "LaTeX", "astrostatistics", "Makie", "units", "physics"]

using Markdown
using InteractiveUtils

# ╔═╡ 7fb0e7b3-e4e3-44d0-bed8-91dec44105f9
begin
    using Pkg: Pkg
    Pkg.activate(; temp = true)

    # TODO: Can remove when https://github.com/MakieOrg/Makie.jl/pull/5623
    # is upstreamed

    # Data viz
    Pkg.add(
        [
            Pkg.PackageSpec(;
                url = "https://github.com/icweaver/Makie.jl",
                subdir = "Makie",
                rev = "units-matrix",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "CairoMakie",
                rev = "ff/breaking-0.25",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "ComputePipeline",
                rev = "ff/breaking-0.25",
            ),
        ]
    )

    Pkg.add(["PlutoUI", "QuadGK", "NumericalIntegration", "InitialMassFunctions", "Unitful", "UnitfulAstro", "UnitfulEquivalences", "PhysicalConstants"])

    using CairoMakie
    using NumericalIntegration
    using QuadGK
    using InitialMassFunctions: Salpeter1955, pdf
    using Unitful
    using Unitful: AbstractQuantity, 𝐋, 𝐓
    using PhysicalConstants.CODATA2018: c_0, h, k_B
    using UnitfulEquivalences: Spectral
    using QuadGK: quadgk
end

# ╔═╡ 3c6d6e2f-d36d-4b19-8214-6c3aa689dd05
begin
    using Pluto: frontmatter
    using PlutoUI: TableOfContents
    using Test: @test
end

# ╔═╡ 1c8b21e7-ca67-44b0-ad7f-b5614bae018b
md"""
## Summary
In this tutorial, we will use the examples of the Planck function and the stellar initial mass function (IMF) to illustrate how to integrate numerically, using the trapezoidal approximation and Gaussian quadrature, respectivelly. We will also explore making a custom Julia object (struct), an instance of which is callable in the same way as a function. In addition, we will encounter unit handling with Unitful.jl, and get a first taste of how to convert between them. Finally, we will use ``\LaTeX`` to make our figure axis labels easy to read.
"""

# ╔═╡ 8caa3ce0-f079-4274-9b6d-deb19a009cd1
md"""
## Packages 📦
"""

# ╔═╡ 035d1d8f-f9f6-426a-9c47-43e0d6476251
md"""
## The Planck function

The Planck function describes how a blackbody radiates energy. We will explore how to find bolometric luminosity using the Planck function in both frequency and wavelength space.

Let's say we have a black-body at 5000 Kelvin. We can find out the total intensity (bolometric) from this object by integrating the Planck function. The simplest way to do this is by approximating the integral using the trapezoid rule. Let's do this first using the frequency definition of the Planck function.

We'll define a photon frequency grid, and evaluate the Planck function at those frequencies. Those will be used to numerically integrate using the trapezoidal rule. By multiplying a collection by a unit (e.g., `[1, 2, 3] * u"m"`), we get a collection of quantities:

!!! note
	There are a wide range of other "quantity" types in Unitful.jl. See the manual for more.
"""

# ╔═╡ 219d613e-1b22-428a-ba2b-2a21ebc9dee4
νs = range(1.0, 3_000.0, length = 1000) * u"THz"

# ╔═╡ d20c3501-a6ff-4ef9-9bcf-78aaadb518f6
md"""
We next compute our `blackbody` function in this frequency space:
"""

# ╔═╡ 3cd6ca5b-a482-4a55-8550-cbd458e378fd
function blackbody(ν::AbstractQuantity{V, inv(𝐓)}, T) where {V}
    return 2 * h * ν^3 / c_0^2 / expm1(h * ν / (k_B * T))
end

# ╔═╡ 81a0a6d9-6e3c-41f7-aa1c-8c48f610474c
md"""
!!! note
	Note the type signature of the first argument in `blackbody`. This specifies that only quantities specifying a frequency will be accepted.
"""

# ╔═╡ 65b93893-3b89-4a6e-8fc9-eb2ed7b7ab20
md"""
And plot the results with Makie.jl for [automatic unit support](https://docs.makie.org/stable/explanations/dim-converts):
"""

# ╔═╡ 91b78fae-df98-419f-8469-ade77cd6df53
md"""
!!! warning "TODO"
	Latex strings and rich strings issue being [tracked upstream here](https://github.com/MakieOrg/Makie.jl/pull/5484#issuecomment-4445754898).
"""

# ╔═╡ 01b1c1b8-e390-42d8-9b7d-325eda02e482
md"""
!!! note "Using LaTeX for axis labels"
	Here, we’ve used ``\LaTeX`` markup to add nice-looking axis labels. To do that, we enclose ``\LaTeX`` markup in the [`L""` string macro](https://github.com/JuliaStrings/LaTeXStrings.jl#latexstrings) that is automatically exported by Makie. Alternatively, we could have also used Makie's [`rich` text formatting](https://docs.makie.org/stable/reference/plots/text#Rich-text) or Julia's built-in [unicode support](https://docs.julialang.org/en/v1/manual/unicode-input/), which is useful for simpler typsetting.
"""

# ╔═╡ 892358f7-2176-4a7f-97b8-870f3a480e1d
md"""
### Integrating using the trapezoid method

Now we numerically integrate using the default trapezoid rule provided by [NumericalIntegration.jl](https://github.com/JuliaMath/NumericalIntegration.jl):
"""

# ╔═╡ 5baa9b54-81d6-406e-a612-1745e52dcc6b
md"""
Now we can do something similar, but for a wavelength grid. We want to integrate over an equivalent wavelength range to the frequency range we did earlier. We can transform the maximum frequency into the corresponding (minimum) wavelength by using the `uconvert()` method, with the addition of an equivalency.
"""

# ╔═╡ 6b4ff602-7806-4a7f-a230-74827206803d
λs = range(
    uconvert(u"Å", νs[end], Spectral()),
    uconvert(u"Å", νs[begin], Spectral());
    length = length(νs)
)

# ╔═╡ 9ccbe520-ab32-48bc-b0a4-4071247161a7
function blackbody(λ::AbstractQuantity{V, 𝐋}, T) where {V}
    return 2 * h * c_0^2 / λ^5 / expm1(h * c_0 / (λ * k_B * T))
end

# ╔═╡ 6ead29e6-301b-46dc-831a-8f8e42bf57ea
bb5000K_ν = blackbody.(νs, 5000u"K") .|> u"erg/cm^2"

# ╔═╡ 5c2f92bd-05b9-4ef7-ac0c-7f3a33c5985e
lines(
    νs, bb5000K_ν;
    axis = (;
        xlabel = L"\nu",
        ylabel = L"I_\nu",
        title = "Planck function in frequency",
    ),
)

# ╔═╡ aa288311-794c-4f03-a2b1-23be8a729702
integrate(νs, bb5000K_ν) |> u"erg/(s * cm^2 * sr)"

# ╔═╡ 69555b33-11ba-4663-b419-5e9f1e1ef558
md"""
!!! note
	Thanks to Julia's [multiple dispatch paradigm](https://docs.julialang.org/en/v1/manual/performance-tips/#Break-functions-into-multiple-definitions), we can define another method for the `blackbody` function which only accepts wavelength quantities. In other words we can use the same spelling of `blackbody` and dispatch to the correct method based on the type of the inputs it receives instead of needing to, e.g., define separate `blackbody_ν` and `blackbody_λ` named functions.
"""

# ╔═╡ ff0c4954-0a78-41aa-8694-9b1b7c28b70b
bb5000K_λ = blackbody.(λs, 5000u"K") .|> u"erg/s / (Å * sr * cm^2)"

# ╔═╡ 88e247da-d821-446d-a371-175290f511bf
let
    fig, ax, p = lines(
        λs, bb5000K_λ;
        axis = (;
            xlabel = "λ",
            ylabel = "I_λ",
            title = "Planck function in wavelength",
            xtickformat = "{:.0f}", # Disable scientific notation on x-axis
        ),
    )

    xlims!(ax, 1.0e3u"Å", 5.5e4u"Å")

    fig
end

# ╔═╡ 05f5b4be-8e4d-4412-a08b-f6f185d5f203
md"""
Integrating the above function then gives:
"""

# ╔═╡ 82786592-f216-431f-895a-31c64dc76332
integrate(λs, bb5000K_λ) |> u"erg/(s * cm^2 * sr)"

# ╔═╡ b987970b-815c-44e5-b481-ebfc36874dc3
md"""
which is within a couple percent of the answer we got in frequency space, despite our bad sampling at small wavelengths!
"""

# ╔═╡ aef35e93-1941-4b84-a179-307b7a8d6a31
md"""
### How to simulate actual observations

Thanks to the composablity of Julia's ecosystem, we can reach for the package we need for the specific part of the simulated observation pipeline that we need. Below is an inexhaustive list of useful packages in this space:

- [Korg.jl](https://ajwheeler.github.io/Korg.jl/stable/): generate intrinsic stellar spectra.
- [DustExctinction.jl](https://juliaastro.org/DustExtinction/stable/): redden/de-redden each spectrum.
- [PhotometricFilters.jl](https://juliaastro.org/PhotometricFilters/stable/): observe each spectrum through their relevant bandpass.
"""

# ╔═╡ fbac0e17-35ae-40c8-896f-b2dcdcfd16d8
md"""
## The stellar initial mass function (IMF)

The stellar initial mass function tells us how many of each mass of stars are formed. In particular, low-mass stars are much more abundant than high-mass stars are. Let’s explore more of the functionality of Julia using this concept.

People generally think of the IMF as a power-law probability density function. In other words, if you count the stars that have been born recently from a cloud of gas, their distribution of masses will follow the IMF. Let’s write a little struct to help us keep track of that:
"""

# ╔═╡ 93fb1f1c-6918-484e-aca9-e656eb5d01ed
struct PowerLawPDF_0
    γ
    B
end

# ╔═╡ 6ffb18ae-b965-4f2f-a5c3-176a7783cca2
md"""
We can now create a `PowerLawPDF` object:
"""

# ╔═╡ bd6ee2a6-7aee-4450-be9e-8612cd23d1e1
PowerLawPDF_0(1, 2)

# ╔═╡ 953efec5-b0f8-4bc7-b920-022dc930cc3c
md"""
Let's make this a bit more convenient to work with by making the following additions:

1. Use the built-in [`Base.@kwdefs`](https://docs.julialang.org/en/v1/base/base/#Base.@kwdef) macro to define some default values, and also allow assignment via keyword:

   ```julia
   Base.@kwdef struct PowerLawPDF
       γ = 1.0
       B = 1.0
   end

   PowerLawPDF() # PowerLawPDF(1.0, 1.0)
   PowerLawPDF(; γ = 2.0) # PowerLawPDF(2.0, 1.0)
   ```

1. [Make the struct parametric](https://docs.julialang.org/en/v1/manual/types/#Parametric-Types) to leverage [type stability](https://docs.julialang.org/en/v1/manual/faq/#man-type-stability) for increased performance:

   ```julia
   Base.@kwdef struct PowerLawPDF{T}
       γ::T = 1.0
       B::T = 1.0
   end

   PowerLawPDF() # PowerLawPDF{Float64}(1.0, 1.0)
   ```

1. Make our struct [callable](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects) so that we can use it like any other function:

   ```julia
   Base.@kwdef struct PowerLawPDF{T}
       γ::T = 1.0
       B::T = 1.0
   end

   (p::PowerLawPDF)(x) = x^p.γ / p.B	

   my_power_law = PowerLawPDF(γ = 2.0) # PowerLawPDF(2.0, 1.0)
   my_power_law(3) # Same as 3^2.0 / 3^1.0 = 9.0
   ```
"""

# ╔═╡ e2238279-1628-4e86-b3e8-9a7ddb9f7b24
begin
    Base.@kwdef struct PowerLawPDF{T}
        γ::T = 1.0
        B::T = 1.0
    end
    (p::PowerLawPDF)(x) = x^p.γ / p.B
end

# ╔═╡ ac7c4b7e-646d-4658-9bf4-0e3d5647fd44
md"""
Let's use this struct to demonstrate the second method of numerical integration that we will explore in this tutorial: Gaussian quadrature.
"""

# ╔═╡ 09c54910-49bd-4f68-89b7-e245dd8243c9
md"""
### Integrating using Gaussian quadrature

In this section, we’ll explore a method of numerical integration that does not require having your sampling grid set-up already. [`QuadGK.quadgk`](https://juliamath.github.io/QuadGK.jl/stable/api/#quadgk) takes a function and both a lower and upper bound, and our `PowerLawPDF` class takes care of this just fine:
"""

# ╔═╡ 05644f4d-5e3c-4e2d-98f4-ba7fec037af9
salpeter_unnormalized = PowerLawPDF(; γ = -2.35)

# ╔═╡ 5140c235-6068-4362-8f8c-61cf2fc1c4bc
md"""
Now we can use our new struct to normalize our IMF given the mass bounds. This amounts to normalizing a probability density function. We’ll use Gaussian quadrature (`quadgk`) to find the integral. `quadgk` returns the numerical value of the integral and its uncertainty. We only care about the numerical value, so we’ll pack the uncertainty into `_` (a placeholder variable). We immediately throw the integral into our IMF object and use it for normalizing!
"""

# ╔═╡ ce09712b-46ab-4366-9088-13ed9f1bac59
B, _ = quadgk(salpeter_unnormalized, 0.01, 100.0)

# ╔═╡ 417552c2-5e66-43c6-8fa8-3a5b60f273ae
salpeter = PowerLawPDF(salpeter_unnormalized.γ, B)

# ╔═╡ ed74a7e7-4754-4f18-a67d-83edfe2a8440
md"""
Next, we plot the resulting normalized curve evaluated over our grid of stellar masses (in log space):
"""

# ╔═╡ cb5af968-e136-457e-826f-eca91825627b
m_grid = logrange(10^-2, 10^2, length = 100)

# ╔═╡ 4b93a32d-61cd-44a9-84cb-c852c4da4b6b
lines(
    m_grid, salpeter.(m_grid);
    axis = (;
        xscale = log10,
        yscale = log10,
        xlabel = "Stellar mass",
        ylabel = "Probability density",
    )
)

# ╔═╡ dbbba33c-fed7-4c82-abd4-9bea6b1a8e51
md"""
### How many more M stars are there than O stars?

Let’s compare the number of M dwarf stars (mass less than 60% solar) created by the IMF, to the number of O stars (mass more than 15 times solar).
"""

# ╔═╡ e7cfec8a-b940-43bb-badb-17e2d13ae0ec
n_M, _ = quadgk(salpeter, 0.01, 0.6)

# ╔═╡ 2e076a6d-0a68-4ea9-92c2-aec19623c92b
n_O, _ = quadgk(salpeter, 15.0, 100.0)

# ╔═╡ 13d91bf7-5f31-4cd3-bef6-4283c37af3f1
r_n = n_M / n_O

# ╔═╡ cef0cdcf-2aa1-43d9-9bab-d35babd17d8d
md"""
There are almost $(round(Int, ceil(r_n; digits=-3))) as many low-mass stars born as there are high-mass stars!
"""

# ╔═╡ 91a68f76-4cab-4393-a5d6-0a7bf9e04828
md"""
### Where is all the mass?

Now let’s compute the relative total masses for all O stars and all M stars born. To do this, weight the Salpeter IMF by mass (i.e., add an extra factor of mass to the integral).

Mathematically, the integral for the M stars is

```math
m^M = ∫_{0.01 M_⊙}^{0.6M_⊙} m \, \text{IMF}(m) \, \text{d}m
```

and it amounts to weighting the probability density function (the IMF) by mass. More generally, you find the value of some property ``\rho`` that depends on ``m`` by calculating

```math
ρ(m)^M = ∫_{0.01M_⊙}^{0.6M_⊙} ρ(m) \, \text{IMF}(m) \, \text{d}m
```
"""

# ╔═╡ e42d809a-60f5-486f-bbd6-cd331556300e
m_M, _ = quadgk(m -> m * salpeter(m), 0.01, 0.6)

# ╔═╡ 48a290c8-bc44-4e48-b41d-9da8ce5f5396
md"""
Similarly for the O stars:
"""

# ╔═╡ 26688cc1-1f14-4322-b186-6614e32de910
m_O, _ = quadgk(m -> m * salpeter(m), 15, 100)

# ╔═╡ 0e9099fd-1dee-4a51-b5f7-bd6475dec8ef
r_m = m_M / m_O

# ╔═╡ 4c62d130-da10-4d14-994d-ae22d8394253
md"""
So about $(round(Int, r_m)) times as much mass is tied up in M stars as in O stars.
"""

# ╔═╡ b5846dbd-59a5-40be-acf1-2418ac80b090
md"""
### InitialMassFunctions.jl

For creating Salpeter and other published IMF functions, the [InitialMassFunctions.jl](https://cgarling.github.io/InitialMassFunctions.jl/stable/) package is available. Here is the above example again, which integrates to within machine precision (`eps()` = $(eps())) of our previous results:
"""

# ╔═╡ d5e03149-6992-43ee-974f-cc947ac5bda5
# Normalized between 0.01 and 100.0 as above
salpeter_imf = Salpeter1955(0.01, 100.0)

# ╔═╡ 305b6672-3281-4759-b39a-80bc3985615d
@test quadgk(salpeter, 0.01, 0.6)[1] ≈ quadgk(x -> pdf(salpeter_imf, x), 0.01, 0.6)[1] atol = eps()

# ╔═╡ 97a89580-f254-4c44-b55f-1d294535f1c9
@test quadgk(salpeter, 15.0, 100.0)[1] ≈ quadgk(x -> pdf(salpeter_imf, x), 15.0, 100.0)[1] atol = eps()

# ╔═╡ 085c3464-b272-479b-bf00-c8f28085176e
md"""
### Extras

- Now compare the total luminosity from all O stars to total luminosity from all M stars. This requires a mass-luminosity relation, like this one which you will use as ``ρ(m)``:
  ```math
  \tilde{L}(\tilde{M}) = \begin{cases}
      0.23 \tilde{M}^{2.3}, & 0.1 < \tilde{M} < 0.43, \\
      \tilde{M}^4, & 0.43 < \tilde{M} < 2, \\
      1.5 \tilde{M}^{3.5}, & 2 < \tilde{M} < 20, \\
      3200 \tilde{M}, & 20 < \tilde{M} < 100.
  \end{cases}
  ```
  where ``\tilde{L} = L/L_⊙`` and ``\tilde{M} = M/M_⊙``.
- Think about which stars are producing most of the light, and which stars have most of the mass. How might this result in difficulty inferring stellar masses from the light they produce? If you’re interested in learning more, see [this review article](https://ned.ipac.caltech.edu/level5/Sept14/Courteau/Courteau_contents.html).
"""

# ╔═╡ b76f44ba-62b3-4e74-9843-43f1765b508f
md"""
## Challenge problems

- Right now, we aren’t worried about the bounds of the power law, but the IMF should drop off to zero probability at masses below .01 solar masses and above 100 solar masses. Modify `PowerLawPDF` in a way that allows both float and array inputs.

- Modify the `PowerLawPDF` class to explicitly use units.

- Derive a relationship between recent star-formation rate and ``H\alpha`` luminosity. In other words, find a value of ``C`` for the function:

  ```math
  \text{SFR}\left[\frac{M_⊙}{\text{yr}}\right] =
  C L_\mathrm{H\alpha} \left[\mathrm{\frac{erg}{s}}\right] \quad .
  ```

- How does this depend on the slope and endpoints of the IMF?

- Take a look at Appendix B of [Hunter & Elmegreen 2004, AJ, 128, 2170](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:astro-ph/0408229)

- What effect does changing the power-law index or upper mass limit of the IMF have on the value of ``C``?

- Predict the effect on the value of ``C`` using a different form of the IMF, like Kroupa or Chabrier (both are lighter on the low-mass end).

- If you’re not tired of IMFs yet, try defining a new class that implements a broken-power-law (Kroupa) or log-parabola (Chabrier) IMF. Perform the same calculations as above.
"""

# ╔═╡ c11661a9-b034-41dc-afd2-4bae2268a037
md"""
# Notebook setup 🔧
"""

# ╔═╡ 035a2be3-1056-4fc5-952e-93fbe07039b3
TableOfContents()

# ╔═╡ 51032b71-7a7f-4918-bef6-be623d6058bd
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ a4c68d00-d91e-4c27-869b-62cdacde8f9c
md"""
# Integration with units and blackbody curves -- Unitful

This tutorial is modified from <https://learn.astropy.org/tutorials/units-and-integration.html>.

_Original authors: Zach Pace, Lia Corrales, Stephanie T. Douglas_

!!! tip "Learning goals"

	- perform numerical integration
	- make trapezoidal approximations
	- use gaussian quadrature
	- define blackbody curves
	- understand how units interact with one another
	- define a Julia struct type
	- explore how callable objects works
	- add ``\LaTeX`` labels to plot figures using the `L"..."` string macro

$(keywords())

!!! warning "Companion content"
	- [learn.JuliaAstro > Unit handling](https://learn.juliaastro.org/tutorials/units/)
"""

# ╔═╡ Cell order:
# ╟─a4c68d00-d91e-4c27-869b-62cdacde8f9c
# ╟─1c8b21e7-ca67-44b0-ad7f-b5614bae018b
# ╟─8caa3ce0-f079-4274-9b6d-deb19a009cd1
# ╠═7fb0e7b3-e4e3-44d0-bed8-91dec44105f9
# ╟─035d1d8f-f9f6-426a-9c47-43e0d6476251
# ╠═219d613e-1b22-428a-ba2b-2a21ebc9dee4
# ╟─d20c3501-a6ff-4ef9-9bcf-78aaadb518f6
# ╠═3cd6ca5b-a482-4a55-8550-cbd458e378fd
# ╠═6ead29e6-301b-46dc-831a-8f8e42bf57ea
# ╟─81a0a6d9-6e3c-41f7-aa1c-8c48f610474c
# ╟─65b93893-3b89-4a6e-8fc9-eb2ed7b7ab20
# ╠═5c2f92bd-05b9-4ef7-ac0c-7f3a33c5985e
# ╟─91b78fae-df98-419f-8469-ade77cd6df53
# ╟─01b1c1b8-e390-42d8-9b7d-325eda02e482
# ╟─892358f7-2176-4a7f-97b8-870f3a480e1d
# ╠═aa288311-794c-4f03-a2b1-23be8a729702
# ╟─5baa9b54-81d6-406e-a612-1745e52dcc6b
# ╠═6b4ff602-7806-4a7f-a230-74827206803d
# ╠═9ccbe520-ab32-48bc-b0a4-4071247161a7
# ╟─69555b33-11ba-4663-b419-5e9f1e1ef558
# ╠═ff0c4954-0a78-41aa-8694-9b1b7c28b70b
# ╟─88e247da-d821-446d-a371-175290f511bf
# ╟─05f5b4be-8e4d-4412-a08b-f6f185d5f203
# ╠═82786592-f216-431f-895a-31c64dc76332
# ╟─b987970b-815c-44e5-b481-ebfc36874dc3
# ╟─aef35e93-1941-4b84-a179-307b7a8d6a31
# ╟─fbac0e17-35ae-40c8-896f-b2dcdcfd16d8
# ╠═93fb1f1c-6918-484e-aca9-e656eb5d01ed
# ╟─6ffb18ae-b965-4f2f-a5c3-176a7783cca2
# ╠═bd6ee2a6-7aee-4450-be9e-8612cd23d1e1
# ╟─953efec5-b0f8-4bc7-b920-022dc930cc3c
# ╠═e2238279-1628-4e86-b3e8-9a7ddb9f7b24
# ╟─ac7c4b7e-646d-4658-9bf4-0e3d5647fd44
# ╟─09c54910-49bd-4f68-89b7-e245dd8243c9
# ╠═05644f4d-5e3c-4e2d-98f4-ba7fec037af9
# ╟─5140c235-6068-4362-8f8c-61cf2fc1c4bc
# ╠═ce09712b-46ab-4366-9088-13ed9f1bac59
# ╠═417552c2-5e66-43c6-8fa8-3a5b60f273ae
# ╟─ed74a7e7-4754-4f18-a67d-83edfe2a8440
# ╠═cb5af968-e136-457e-826f-eca91825627b
# ╠═4b93a32d-61cd-44a9-84cb-c852c4da4b6b
# ╟─dbbba33c-fed7-4c82-abd4-9bea6b1a8e51
# ╠═e7cfec8a-b940-43bb-badb-17e2d13ae0ec
# ╠═2e076a6d-0a68-4ea9-92c2-aec19623c92b
# ╠═13d91bf7-5f31-4cd3-bef6-4283c37af3f1
# ╟─cef0cdcf-2aa1-43d9-9bab-d35babd17d8d
# ╟─91a68f76-4cab-4393-a5d6-0a7bf9e04828
# ╠═e42d809a-60f5-486f-bbd6-cd331556300e
# ╟─48a290c8-bc44-4e48-b41d-9da8ce5f5396
# ╠═26688cc1-1f14-4322-b186-6614e32de910
# ╠═0e9099fd-1dee-4a51-b5f7-bd6475dec8ef
# ╟─4c62d130-da10-4d14-994d-ae22d8394253
# ╟─b5846dbd-59a5-40be-acf1-2418ac80b090
# ╠═d5e03149-6992-43ee-974f-cc947ac5bda5
# ╠═305b6672-3281-4759-b39a-80bc3985615d
# ╠═97a89580-f254-4c44-b55f-1d294535f1c9
# ╟─085c3464-b272-479b-bf00-c8f28085176e
# ╟─b76f44ba-62b3-4e74-9843-43f1765b508f
# ╟─c11661a9-b034-41dc-afd2-4bae2268a037
# ╠═035a2be3-1056-4fc5-952e-93fbe07039b3
# ╟─51032b71-7a7f-4918-bef6-be623d6058bd
# ╠═3c6d6e2f-d36d-4b19-8214-6c3aa689dd05

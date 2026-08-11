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

# ╔═╡ 8d51b4dd-fca2-4ac0-a706-90261ce1035c
begin
    using Pkg: Pkg
    Pkg.activate(; temp = true)

    # Data viz
    Pkg.add(
        [
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "Makie",
                rev = "ff/breaking-0.25",
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

    Pkg.add(["DynamicQuantities", "PlutoUI", "QuadGK", "NumericalIntegration", "InitialMassFunctions"])

    using DynamicQuantities
    using CairoMakie
    using DynamicQuantities.Constants: c, h, k_B
    using NumericalIntegration
    using QuadGK
    using InitialMassFunctions: Salpeter1955, pdf
end

# ╔═╡ 00b0c01d-3d3c-4822-9738-e33efe95cef8
begin
    using Pluto: frontmatter
    using PlutoUI: TableOfContents
    using Test: @test
end

# ╔═╡ 4a95f7da-6c3e-4290-9985-9f5209c9e5e0
md"""
## Summary
In this tutorial, we will use the examples of the Planck function and the stellar initial mass function (IMF) to illustrate how to integrate numerically, using the trapezoidal approximation and Gaussian quadrature, respectivelly. We will also explore making a custom Julia object (struct), an instance of which is callable in the same way as a function. In addition, we will encounter unit handling with DynamicQuantities.jl, and get a first taste of how to convert between them. Finally, we will use ``\LaTeX`` to make our figure axis labels easy to read.
"""

# ╔═╡ cd733dc4-a6b7-415f-b619-ebced71f8e27
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
	There are a wide range of other "quantity" types in DynamicQuantities.jl. See the manual for more.
"""

# ╔═╡ 219d613e-1b22-428a-ba2b-2a21ebc9dee4
νs = range(1.0, 3_000.0; length = 1000) * u"THz"

# ╔═╡ c7ba06c3-121a-431b-9d73-39c86d4d154c
md"""
We next compute our `blackbody` function in this frequency space:
"""

# ╔═╡ a27c20a0-aa84-4d2a-95e6-1c830be322f8
blackbody_ν(ν, T) = 2 * h * ν^3 / c^2 / expm1(h * ν / (k_B * T))

# ╔═╡ 6ead29e6-301b-46dc-831a-8f8e42bf57ea
bb5000K_ν = blackbody_ν.(νs, 5_000u"K")

# ╔═╡ 7b87bf56-9465-4ba6-a461-42af2208a99a
md"""
And plot the results with Makie.jl for [automatic unit support](https://docs.makie.org/stable/explanations/dim-converts):
"""

# ╔═╡ 5cd159cc-3939-4827-883f-d61f9f4931c5
lines(
    νs, bb5000K_ν;
    axis = (;
        xlabel = L"\nu",
        ylabel = L"I_\nu",
        title = "Planck function in frequency",
        dim1_conversion = Makie.DQConversion(us"THz"),
        dim2_conversion = Makie.DQConversion(us"erg/Hz/s/sr/cm^2"),
    ),
)

# ╔═╡ 25efd32c-5beb-47d2-9d75-32ce1cf9960b
md"""
!!! warning "TODO"
	Latex strings and rich strings issue being [tracked upstream here](https://github.com/MakieOrg/Makie.jl/pull/5484#issuecomment-4445754898).
"""

# ╔═╡ 93f0e1a7-c749-4a6c-b40b-30e93a099390
md"""
!!! note "Using LaTeX for axis labels"
	Here, we’ve used ``\LaTeX`` markup to add nice-looking axis labels. To do that, we enclose ``\LaTeX`` markup in the [`L""` string macro](https://github.com/JuliaStrings/LaTeXStrings.jl#latexstrings) that is automatically exported by Makie. Alternatively, we could have also used Makie's [`rich` text formatting](https://docs.makie.org/stable/reference/plots/text#Rich-text) or Julia's built-in [unicode support](https://docs.julialang.org/en/v1/manual/unicode-input/), which is useful for simpler typsetting.
"""

# ╔═╡ 15dfa732-d707-4556-95fa-fdc116e178a3
md"""
### Integrating using the trapezoid method

Now we numerically integrate using the default trapezoid rule provided by [NumericalIntegration.jl](https://github.com/JuliaMath/NumericalIntegration.jl):
"""

# ╔═╡ aa288311-794c-4f03-a2b1-23be8a729702
integrate(νs, bb5000K_ν) |> us"erg/(s * cm^2 * sr)"

# ╔═╡ 5baa9b54-81d6-406e-a612-1745e52dcc6b
md"""
Now we can do something similar, but for a wavelength grid. We want to integrate over an equivalent wavelength range to the frequency range we did earlier. We can transform the maximum frequency into the corresponding (minimum) wavelength, and vice-versa, using the relationship ``c = λν`` to create our uniform grid in wavelength space:
"""

# ╔═╡ ff27f02a-2efe-4b99-a9ab-8c851c2cc55b
λs = range(c / νs[end], c / νs[begin]; length = length(νs))

# ╔═╡ 680ec021-58e9-4518-b992-a31d627a4d57
md"""
We next compute and plot the corresponding blackbody curve:
"""

# ╔═╡ cc3fd5c1-12b1-4e94-afe2-7cc440dd7bf0
blackbody_λ(λ, T) = 2 * h * c^2 / λ^5 / expm1(h * c / (λ * k_B * T))

# ╔═╡ ff0c4954-0a78-41aa-8694-9b1b7c28b70b
bb5000K_λ = blackbody_λ.(λs, 5000u"K")

# ╔═╡ 07a699ba-9a78-4baa-97d6-77388c3578ae
md"""
!!! warning "TODO"
	[See here](https://github.com/JuliaPhysics/DynamicQuantities.jl/issues/197#issuecomment-4445192280) for a discussion on potentially making this trait-based so that a single `blackbody` function can be called instead.
"""

# ╔═╡ 88e247da-d821-446d-a371-175290f511bf
let
    fig, ax, p = lines(
        λs, bb5000K_λ;
        axis = (;
            xlabel = "λ",
            ylabel = "I_λ",
            title = "Planck function in wavelength",
            xtickformat = "{:.0f}", # Disable scientific notation on x-axis
            dim1_conversion = Makie.DQConversion(us"Å"),
            dim2_conversion = Makie.DQConversion(us"erg/(s * Å * sr * cm^2)"),
        ),
    )

    xlims!(ax, 1.0e3u"Å", 5.5e4u"Å")

    fig
end

# ╔═╡ 23f12bab-9cc9-48ac-854c-ae06b4ad7b85
md"""
Integrating the above function then gives:
"""

# ╔═╡ 82786592-f216-431f-895a-31c64dc76332
integrate(λs, bb5000K_λ) |> us"erg/(s * cm^2 * sr)"

# ╔═╡ 5a4a0548-f2b9-4d68-8261-1b76f3755ede
md"""
which is within a couple percent of the answer we got in frequency space, despite our bad sampling at small wavelengths!
"""

# ╔═╡ 0e9d0c6b-f4e4-431b-bd31-41490a63aa9c
md"""
### How to simulate actual observations

Thanks to the composablity of Julia's ecosystem, we can reach for the package we need for the specific part of the simulated observation pipeline that we need. Below is an inexhaustive list of useful packages in this space:

- [Korg.jl](https://ajwheeler.github.io/Korg.jl/stable/): generate intrinsic stellar spectra.
- [DustExctinction.jl](https://juliaastro.org/DustExtinction/stable/): redden/de-redden each spectrum.
- [PhotometricFilters.jl](https://juliaastro.org/PhotometricFilters/stable/): observe each spectrum through their relevant bandpass.
"""

# ╔═╡ 11773ea5-a065-4893-a009-26524462c25c
md"""
## The stellar initial mass function (IMF)

The stellar initial mass function tells us how many of each mass of stars are formed. In particular, low-mass stars are much more abundant than high-mass stars are. Let’s explore more of the functionality of Julia using this concept.

People generally think of the IMF as a power-law probability density function. In other words, if you count the stars that have been born recently from a cloud of gas, their distribution of masses will follow the IMF. Let’s write a little struct to help us keep track of that:
"""

# ╔═╡ de615b56-ee1e-4fd7-96fb-a0a6c0cb51d3
struct PowerLawPDF_0
    γ
    B
end

# ╔═╡ f42bcbe2-6b97-471d-abc1-8bb4ef583b63
md"""
We can now create a `PowerLawPDF` object:
"""

# ╔═╡ 7aacf5ad-d325-4903-be27-5c50955c90dd
PowerLawPDF_0(1, 2)

# ╔═╡ 7ffacce1-dd54-42aa-b838-daad690e7669
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

# ╔═╡ 3fb419ea-dfc4-4a4b-96c7-cdb6bb533853
begin
    Base.@kwdef struct PowerLawPDF{T}
        γ::T = 1.0
        B::T = 1.0
    end
    (p::PowerLawPDF)(x) = x^p.γ / p.B
end

# ╔═╡ 5c60e813-f63a-4046-8de5-3d8654d7ef0c
md"""
Let's use this struct to demonstrate the second method of numerical integration that we will explore in this tutorial: Gaussian quadrature.
"""

# ╔═╡ 06b8e0db-7683-47d0-a0cf-ab2b4c262a31
md"""
### Integrating using Gaussian quadrature

In this section, we’ll explore a method of numerical integration that does not require having your sampling grid set-up already. [`QuadGK.quadgk`](https://juliamath.github.io/QuadGK.jl/stable/api/#quadgk) takes a function and both a lower and upper bound, and our `PowerLawPDF` class takes care of this just fine:
"""

# ╔═╡ c9b2f21e-5d92-4058-b96b-43395facd11e
salpeter_unnormalized = PowerLawPDF(; γ = -2.35)

# ╔═╡ 0c031a31-4a63-4c98-91cc-b0173efaed7e
md"""
Now we can use our new struct to normalize our IMF given the mass bounds. This amounts to normalizing a probability density function. We’ll use Gaussian quadrature (`quadgk`) to find the integral. `quadgk` returns the numerical value of the integral and its uncertainty. We only care about the numerical value, so we’ll pack the uncertainty into `_` (a placeholder variable). We immediately throw the integral into our IMF object and use it for normalizing!
"""

# ╔═╡ 0fafa682-38cf-4caa-9535-afe9eaf6f1f7
B, _ = quadgk(salpeter_unnormalized, 0.01, 100.0)

# ╔═╡ e768ae7e-125b-4bc0-8398-aaa1207ff557
salpeter = PowerLawPDF(salpeter_unnormalized.γ, B)

# ╔═╡ 3e9152d6-ec72-4a82-a9e1-573b76011998
md"""
Next, we plot the resulting normalized curve evaluated over our grid of stellar masses (in log space):
"""

# ╔═╡ 5d5cdec6-51db-4747-890b-959ce7949550
m_grid = logrange(10^-2, 10^2, length = 100)

# ╔═╡ ddbce299-cd22-4d87-ae5f-c44b50f46f1a
lines(
    m_grid, salpeter.(m_grid);
    axis = (;
        xscale = log10,
        yscale = log10,
        xlabel = "Stellar mass",
        ylabel = "Probability density",
    )
)

# ╔═╡ f4648f97-de9f-44e7-857f-aff8d0420b76
md"""
### How many more M stars are there than O stars?

Let’s compare the number of M dwarf stars (mass less than 60% solar) created by the IMF, to the number of O stars (mass more than 15 times solar).
"""

# ╔═╡ 01c05dc0-a1b6-4782-8a9a-da2b9ad5cfa8
n_M, _ = quadgk(salpeter, 0.01, 0.6)

# ╔═╡ d4ca3c29-fb21-4118-a312-df956454c90a
n_O, _ = quadgk(salpeter, 15.0, 100.0)

# ╔═╡ a6d251e8-6c12-4847-8336-f5fa1040c217
r_n = n_M / n_O

# ╔═╡ 43a8fe9a-4bb8-479f-ac10-08ee28d1d4ec
md"""
There are almost $(round(Int, ceil(r_n; digits=-3))) as many low-mass stars born as there are high-mass stars!
"""

# ╔═╡ 57cb5e9f-ba5e-4eaf-94cd-7538a55fce59
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

# ╔═╡ b90e4072-243c-450f-995b-5ffed001d634
m_M, _ = quadgk(m -> m * salpeter(m), 0.01, 0.6)

# ╔═╡ 5a97b7d1-f331-49aa-b890-103bfb825d67
md"""
Similarly for the O stars:
"""

# ╔═╡ f6ebbfaa-08ed-445b-ac81-e55c58cb3081
m_O, _ = quadgk(m -> m * salpeter(m), 15, 100)

# ╔═╡ e5340a7a-9654-414a-be62-1cc416ca1fdc
r_m = m_M / m_O

# ╔═╡ 8ce0941a-3405-46aa-a545-8636c4b8b6f6
md"""
So about $(round(Int, r_m)) times as much mass is tied up in M stars as in O stars.
"""

# ╔═╡ a1423be7-8649-446a-bc12-2988a3b1734a
md"""
### InitialMassFunctions.jl

For creating Salpeter and other published IMF functions, the [InitialMassFunctions.jl](https://cgarling.github.io/InitialMassFunctions.jl/stable/) package is available. Here is the above example again, which integrates to within machine precision (`eps()` = $(eps())) of our previous results:
"""

# ╔═╡ ddcfecd0-77af-4892-910c-ea917c9a8377
# Normalized between 0.01 and 100.0 as above
salpeter_imf = Salpeter1955(0.01, 100.0)

# ╔═╡ 1e72e3f5-343c-40bd-9dbb-3bc4d9194103
@test quadgk(salpeter, 0.01, 0.6)[1] ≈ quadgk(x -> pdf(salpeter_imf, x), 0.01, 0.6)[1] atol = eps()

# ╔═╡ 77e0f5da-3c0d-4276-ba09-8b465862509e
@test quadgk(salpeter, 15.0, 100.0)[1] ≈ quadgk(x -> pdf(salpeter_imf, x), 15.0, 100.0)[1] atol = eps()

# ╔═╡ 8f4baa15-9471-4f6a-be8b-ba47112128dc
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

# ╔═╡ 40ae6e9e-7d72-4363-a94c-ff1dec656ff6
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

# ╔═╡ 096ede9d-7900-4cdc-b66f-ef3703ac5898
md"""
# Notebook setup 🔧
"""

# ╔═╡ 508c13fe-39cc-47e9-a792-92ca463de0dd
TableOfContents()

# ╔═╡ be7aa3a4-7adc-4138-84de-e52ea23be4b9
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ 3080205e-60c5-11ef-3df9-cd4326814e34
md"""
# Integration with units and blackbody curves -- DQ

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
# ╟─3080205e-60c5-11ef-3df9-cd4326814e34
# ╟─4a95f7da-6c3e-4290-9985-9f5209c9e5e0
# ╟─cd733dc4-a6b7-415f-b619-ebced71f8e27
# ╠═8d51b4dd-fca2-4ac0-a706-90261ce1035c
# ╟─035d1d8f-f9f6-426a-9c47-43e0d6476251
# ╠═219d613e-1b22-428a-ba2b-2a21ebc9dee4
# ╟─c7ba06c3-121a-431b-9d73-39c86d4d154c
# ╠═a27c20a0-aa84-4d2a-95e6-1c830be322f8
# ╠═6ead29e6-301b-46dc-831a-8f8e42bf57ea
# ╟─7b87bf56-9465-4ba6-a461-42af2208a99a
# ╠═5cd159cc-3939-4827-883f-d61f9f4931c5
# ╟─25efd32c-5beb-47d2-9d75-32ce1cf9960b
# ╟─93f0e1a7-c749-4a6c-b40b-30e93a099390
# ╟─15dfa732-d707-4556-95fa-fdc116e178a3
# ╠═aa288311-794c-4f03-a2b1-23be8a729702
# ╟─5baa9b54-81d6-406e-a612-1745e52dcc6b
# ╠═ff27f02a-2efe-4b99-a9ab-8c851c2cc55b
# ╟─680ec021-58e9-4518-b992-a31d627a4d57
# ╠═cc3fd5c1-12b1-4e94-afe2-7cc440dd7bf0
# ╠═ff0c4954-0a78-41aa-8694-9b1b7c28b70b
# ╟─07a699ba-9a78-4baa-97d6-77388c3578ae
# ╠═88e247da-d821-446d-a371-175290f511bf
# ╟─23f12bab-9cc9-48ac-854c-ae06b4ad7b85
# ╠═82786592-f216-431f-895a-31c64dc76332
# ╟─5a4a0548-f2b9-4d68-8261-1b76f3755ede
# ╟─0e9d0c6b-f4e4-431b-bd31-41490a63aa9c
# ╟─11773ea5-a065-4893-a009-26524462c25c
# ╠═de615b56-ee1e-4fd7-96fb-a0a6c0cb51d3
# ╟─f42bcbe2-6b97-471d-abc1-8bb4ef583b63
# ╠═7aacf5ad-d325-4903-be27-5c50955c90dd
# ╟─7ffacce1-dd54-42aa-b838-daad690e7669
# ╠═3fb419ea-dfc4-4a4b-96c7-cdb6bb533853
# ╟─5c60e813-f63a-4046-8de5-3d8654d7ef0c
# ╟─06b8e0db-7683-47d0-a0cf-ab2b4c262a31
# ╠═c9b2f21e-5d92-4058-b96b-43395facd11e
# ╟─0c031a31-4a63-4c98-91cc-b0173efaed7e
# ╠═0fafa682-38cf-4caa-9535-afe9eaf6f1f7
# ╠═e768ae7e-125b-4bc0-8398-aaa1207ff557
# ╟─3e9152d6-ec72-4a82-a9e1-573b76011998
# ╠═5d5cdec6-51db-4747-890b-959ce7949550
# ╠═ddbce299-cd22-4d87-ae5f-c44b50f46f1a
# ╟─f4648f97-de9f-44e7-857f-aff8d0420b76
# ╠═01c05dc0-a1b6-4782-8a9a-da2b9ad5cfa8
# ╠═d4ca3c29-fb21-4118-a312-df956454c90a
# ╠═a6d251e8-6c12-4847-8336-f5fa1040c217
# ╟─43a8fe9a-4bb8-479f-ac10-08ee28d1d4ec
# ╟─57cb5e9f-ba5e-4eaf-94cd-7538a55fce59
# ╠═b90e4072-243c-450f-995b-5ffed001d634
# ╟─5a97b7d1-f331-49aa-b890-103bfb825d67
# ╠═f6ebbfaa-08ed-445b-ac81-e55c58cb3081
# ╠═e5340a7a-9654-414a-be62-1cc416ca1fdc
# ╟─8ce0941a-3405-46aa-a545-8636c4b8b6f6
# ╟─a1423be7-8649-446a-bc12-2988a3b1734a
# ╠═ddcfecd0-77af-4892-910c-ea917c9a8377
# ╠═1e72e3f5-343c-40bd-9dbb-3bc4d9194103
# ╠═77e0f5da-3c0d-4276-ba09-8b465862509e
# ╟─8f4baa15-9471-4f6a-be8b-ba47112128dc
# ╟─40ae6e9e-7d72-4363-a94c-ff1dec656ff6
# ╟─096ede9d-7900-4cdc-b66f-ef3703ac5898
# ╠═508c13fe-39cc-47e9-a792-92ca463de0dd
# ╟─be7aa3a4-7adc-4138-84de-e52ea23be4b9
# ╠═00b0c01d-3d3c-4822-9738-e33efe95cef8

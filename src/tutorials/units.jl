### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> tags = ["units", "histograms", "radio astronomy", "data cubes"]
#> title = "Unit handling"
#> description = "Work with units in astrophysical calculations."
#> date = "2025-11-19"
#> layout = "layout.jlhtml"
#> image = "/assets/units.png"
#> order = 1

using Markdown
using InteractiveUtils

# ╔═╡ fd88a6c1-0abe-4a5a-9414-bb15730c9d18
begin
    # Some packages are pinned to development branches until Makie v0.25 is merged:
    # https://github.com/MakieOrg/Makie.jl/pull/5673

    # Not used directly, but imported so that Pluto keeps them pinned:
    import ComputePipeline

    # Analysis
    using StatsBase: mean
    using Distributions: Normal
    using Random: Xoshiro

    # Units -- DynamicQuantities
    using DynamicQuantities: DynamicQuantities as DQ
    using DynamicQuantities: SymbolicConstants as C_DQ

    # Units -- Unitful
    using Unitful: Unitful as U
    using PhysicalConstants: CODATA2018 as C_U
    using UnitfulEquivalences: UnitfulEquivalences as UE
    import UnitfulAstro
    import DimensionfulAngles

    # Data viz
    using CairoMakie: Colorbar, stephist, heatmap
    using Makie: Makie as M
end

# ╔═╡ 90fd3ea9-e115-41a5-b020-b8ade6cc6398
begin
    using PlutoUI: TableOfContents, details
    using HypertextLiteral: @htl
end

# ╔═╡ c6ad0267-65d1-4372-a538-22acd9b5d02b
md"""
# Using units in astrophysical calculations

This notebook is modified from <https://learn.astropy.org/tutorials/quantities.html>.

_Original authors: Ana Bonaca, Erik Tollerud, Jonathan Foster, Lia Corrales, Kris Stern, Stephanie T. Douglas_

!!! tip "Learning goals"
    - Estimate a hypothetical galaxy's mass with units
    - Take advantage of constants in the various units packages
    - Print formatted unit strings
    - Plot objects with unit labels, using Makie.jl
    - Do math with units
    - Convert quantities
    - Convert between wavelength and energy
    - Write functions that take objects with units instead of plain arrays
    - Make synthetic radio observations
    - Use objects with units such as data cubes to facilitate a full derivation of the total mass of a molecular cloud

!!! warning "Companion content"
    Content here.
"""

# ╔═╡ 2670973e-750c-47f0-9c6a-a58aadf87682
md"""
## Summary
In this tutorial we present some examples showing how objects with units can make astrophysics calculations easier. The examples include calculating the mass of a galaxy from its velocity dispersion and determining masses of molecular clouds from ``\mathrm{CO}`` intensity maps. We end with an example of good practices for using quantities in functions you might distribute to other people.

This notebook will cover the main usage for the two most common packages for units in Julia: [DynamicQuantities.jl](https://juliaphysics.github.io/DynamicQuantities.jl/stable/) and [Unitful.jl](https://juliaphysics.github.io/Unitful.jl/stable/). Each package has their own trade-offs, which can be read [about here](https://discourse.julialang.org/t/ann-dynamicquantities-jl-type-stable-physical-quantities/99963). For convenience, the usage for each package will be shown side-by-side in this notebook (if your screen is wide enough), with variables and functions using each appended with `_DQ` and `_U`, respectively. Similarly, functions from each package will be qualified with `DQ` or `U`, respectively.
"""

# ╔═╡ 05b485e7-115a-4dbb-aa73-0ca6ace2f5c0
md"""
### Packages 📦
"""

# ╔═╡ 1d2293bd-a236-4b41-a9d7-9c27463b5062
Reff_DQ = 29 * DQ.us"Constants.pc"

# ╔═╡ 426c37ea-83dd-4ef4-8a19-f0e16536c034
Reff_U = 29 * U.u"pc"

# ╔═╡ 1b239188-bba1-44d1-bc9c-10b60f762e0d
DQ.ustrip(DQ.u"Constants.pc", Reff_DQ), DQ.dimension(Reff_DQ)

# ╔═╡ f2c441b3-374b-439d-ba8c-4b28c0765d45
U.ustrip(U.u"pc", Reff_U), U.unit(Reff_U)

# ╔═╡ 7ba07e6b-34e0-4206-9cf7-5006c43635a6
DQ.uconvert(DQ.us"km", Reff_DQ)

# ╔═╡ d6e4fd70-8750-457f-a17f-9a044c25c482
U.uconvert(U.u"m", Reff_U)

# ╔═╡ 2c7a71f6-e618-40db-80b7-1ad5f09e59d5
Reff_DQ |> DQ.us"km"

# ╔═╡ 2a11aed4-6516-40dc-94fd-30aa94e76744
Reff_U |> U.u"km"

# ╔═╡ 984fcb2e-049b-4b35-a4b6-5c8c686921a3
md"""
!!! note "Symbolic units"
    Note the use of `us""` for DynamicQuantities.jl, which is required for working with [symbolic units](https://juliaphysics.github.io/DynamicQuantities.jl/stable/symbolic_units/). Will will leverage this later in the notebook to make working with angles more convenient.
"""

# ╔═╡ c5b4f340-c774-4f09-af4c-f326afce5de3
md"""
Next, we'll first create a synthetic dataset of radial velocity measurements, assuming a normal distribution with a mean velocity of 206 km/s and a velocity dispersion of 4.3 km/s:
"""

# ╔═╡ 2da8169e-a39d-4314-87c0-60b57f80ae96
v̄, σ_in = 206, 4.3

# ╔═╡ 892f19a6-af21-4353-9987-9de795ba7ad7
v_DQ = rand(Xoshiro(0), Normal(v̄, σ_in), 500) * DQ.us"km/s"

# ╔═╡ 64774c2c-0991-48b2-af4a-ced8fedd115b
v_U = rand(Xoshiro(0), Normal(v̄, σ_in), 500) * U.u"km/s"

# ╔═╡ 9f14f051-d1e9-4638-8008-840a886b50ab
stephist(v_DQ)

# ╔═╡ 0cd8f84d-a6aa-46c7-9e13-535e8e0cea1f
stephist(v_U)

# ╔═╡ 2d236d32-faf4-41ef-84de-1d40c02fb238
σ_DQ = √(sum((v_DQ .- mean(v_DQ)) .^ 2) / length(v_DQ))

# ╔═╡ b388f8be-3074-443c-97a9-72eb7098b5e3
σ_U = √(sum((v_U .- mean(v_U)) .^ 2) / length(v_U))

# ╔═╡ a3baf049-5419-4773-8b7f-0ba07a9f1728
M_DQ = 4 * σ_DQ^2 * Reff_DQ / C_DQ.G

# ╔═╡ efa15492-4811-4896-8947-b90191651952
M_U = 4 * σ_U^2 * Reff_U / C_U.G

# ╔═╡ bb2622f9-98a5-48ca-92cb-bbce84d797f1
M_DQ .|> (DQ.us"Constants.M_sun", DQ.us"g")

# ╔═╡ 14909ad6-eb42-48ab-afa7-8b40b700fe3b
M_U .|> (U.u"Msun", U.u"g")

# ╔═╡ 9e13bde0-eccf-4d21-b603-fd186e87b7d0
(log10 ∘ DQ.ustrip)(DQ.u"Constants.M_sun", M_DQ)

# ╔═╡ 77a81c41-0034-4a18-bb70-d2a1f9379593
(log10 ∘ U.ustrip)(U.u"Msun", M_U)

# ╔═╡ a13de611-145a-40da-9414-8f3f8f85ad98
(log10 ∘ DQ.ustrip)(M_DQ)

# ╔═╡ 1e24cdfa-f4a9-4f9e-b95d-af83e724cd90
(log10 ∘ U.ustrip)(M_U)

# ╔═╡ 095e7efe-f4c6-4a8a-a029-03ce97cd15bd
log10(M_DQ)

# ╔═╡ e4a98c65-70de-42fc-b490-f019dd93d3fb
log10(M_U)

# ╔═╡ 2648b454-2762-4941-b13c-2303ffcd6521
let
    sol = details(
        "Example solution",
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
        """
    )

    md"""
    !!! warn "Exercise"

        Use Kepler's law in the form given below to determine the (circular) orbital speed of the Earth around the Sun in km/s:

        ```math
        v = \sqrt{\frac{G M_⊙}{r}}
        ```

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

# ╔═╡ 2c76c406-d15c-4271-baa4-ebebf5299429
d_DQ, T_ex_DQ = 250 * DQ.us"Constants.pc", 25 * DQ.u"K"

# ╔═╡ be30f3c4-0c20-4949-bcc8-f3812b39befc
d_U, T_ex_U = 250 * U.u"pc", 25 * U.u"K"

# ╔═╡ aab4113d-d7ad-4521-a208-e192ad218cf0
begin
    # 1D coordinate quantities
    ras_DQ = range(52, 52.5; length = 100) * DQ.us"deg"
    decs_DQ = range(0, 0.5; length = 100) * DQ.us"deg"
    vs_DQ = range(0, 30; length = 300) * DQ.us"km/s"
end;

# ╔═╡ af22e5da-ff9b-4278-95d9-491f84a55add
begin
    # 1D coordinate quantities
    ras_U = range(52, 52.5; length = 100) * U.u"deg"
    decs_U = range(0, 0.5; length = 100) * U.u"deg"
    vs_U = range(0, 30; length = 300) * U.u"km/s"
end;

# ╔═╡ 539290af-787d-4deb-928c-50e4e9f28173
data_DQ = let
    # Cloud's center
    ra_0 = 52.25 * DQ.u"deg"
    d_0 = 0.25 * DQ.u"deg"
    v_0 = 15 * DQ.u"km/s"

    # Cloud's size
    ra_σ = 3 * DQ.u"arcmin"
    d_σ = 4 * DQ.u"arcmin"
    v_σ = 3 * DQ.us"km/s"

    A = [
        exp(
            -0.5 * (
                ((ra - ra_0) / ra_σ)^2
                    + ((d - d_0) / d_σ)^2
                    + ((v - v_0) / v_σ)^2
            )
        )
            for ra in ras_DQ, d in decs_DQ, v in vs_DQ
    ] * DQ.us"K"
end

# ╔═╡ 86d608ec-aca2-4738-9c3e-656ccd562da6
data_U = let
    # Cloud's center
    ra_0 = 52.25 * U.u"deg"
    d_0 = 0.25 * U.u"deg"
    v_0 = 15 * U.u"km/s"

    # Cloud's size
    ra_σ = 3 * U.u"arcminute"
    d_σ = 4 * U.u"arcminute"
    v_σ = 3 * U.u"km/s"

    A = [
        exp(
            -0.5 * (
                ((ra - ra_0) / ra_σ)^2
                    + ((d - d_0) / d_σ)^2
                    + ((v - v_0) / v_σ)^2
            )
        )
            for ra in ras_U, d in decs_U, v in vs_U
    ] * U.u"K"
end

# ╔═╡ e86d6cf7-2274-403a-b1db-017973f33fb7
md"""
!!! note
    The units of the exponential are dimensionless, so we multiplied the data cube by ``\mathrm{K}`` to get brightness temperature units. As an aside for experts, we're setting up our artificial cube on the main-beam temperature scale ``\left(T_\text{MB}\right)``, which is the closest we can normally get to the actual brightness temperature of our source.
"""

# ╔═╡ 1f6a0d6c-832b-474d-8fbf-a9de6e821d80
# Average pixel size
# This is only right if dec ~ 0, because of the cos(dec) factor.
Δra_DQ, Δdec_DQ = (
    (maximum(ras_DQ) - minimum(ras_DQ)) / length(ras_DQ), # Typed \Delta<TAB>
    (maximum(decs_DQ) - minimum(decs_DQ)) / length(decs_DQ),
)

# ╔═╡ efa21003-ee7a-4834-b6f7-7625a7820f70
# Average pixel size
# This is only right if dec ~ 0, because of the cos(dec) factor.
Δra_U, Δdec_U = (
    (maximum(ras_U) - minimum(ras_U)) / length(ras_U), # Typed \Delta<TAB>
    (maximum(decs_U) - minimum(decs_U)) / length(decs_U),
)

# ╔═╡ 6fa07efa-fd25-4bca-bd18-8ba5c557d146
# Average velocity bin width
Δv_DQ = (maximum(vs_DQ) - minimum(vs_DQ)) / length(vs_DQ)

# ╔═╡ 95bf572e-23a0-4554-8ee3-b7932721e850
# Average velocity bin width
Δv_U = (maximum(vs_U) - minimum(vs_U)) / length(vs_U)

# ╔═╡ 6ddeb42d-37a5-48a3-ac82-47e4ec2ca541
intcloud_DQ = let
    A = data_DQ * Δv_DQ
    sum(eachslice(A; dims = 3))
    # sum(A, dims = 3; init = zero(first(A)))[:, :, begin]
end

# ╔═╡ 9537c91e-03a7-44d7-95ae-e60ea8d5496a
intcloud_U = let
    A = data_U * Δv_U
    sum(eachslice(A; dims = 3))
end

# ╔═╡ 4d4794fb-804d-4b4a-8c32-a32d88e43e30
md"""
!!! todo
    Improve `Base.sum` support for DQ. Discussion here: <https://github.com/JuliaPhysics/DynamicQuantities.jl/issues/76#issuecomment-3614719247>

    This works too, but is more verbose:

    ```julia
    sum(A, dims = 3; init = zero(first(A)))[:, :, begin]
    ```
"""

# ╔═╡ 1fa8010d-9a5d-4297-9665-9fa8795ef5f7
md"""
!!! note
    Radio astronomers use a rather odd set of units ``[\mathrm{K\, km/s}]`` for integrated intensity (that is, summing all the emission from a line over velocity).
"""

# ╔═╡ f2b222ec-0783-487e-9c52-835976a555b6
# let
#     fig, ax, p = heatmap(
#         ras_DQ, decs_DQ, intcloud_DQ,
#         axis = (;
#             xreversed = true,
#             xlabel = "RA",
#             ylabel = "Dec",
#         ),
#         colormap = :cividis,
#     )

#     Colorbar(fig[1, 2], p; label = "Intensity")

#     fig
# end

# ╔═╡ f2c3768f-1c81-456d-ba81-6a91fc09e81b
# let
#     fig, ax, p = heatmap(
#         ras_U, decs_U, intcloud_U,
#         axis = (;
#             xreversed = true,
#             xlabel = "RA",
#             ylabel = "Dec",
#         ),
#         colormap = :cividis,
#     )

#     Colorbar(fig[1, 2], p; label = "Intensity")

#     fig
# end

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

# ╔═╡ bf618161-ef96-4445-8fc8-25dc5f662242
λ_13_DQ, λ_18_DQ = (2.60076, 2.73079) .* DQ.us"mm"

# ╔═╡ 182861a3-2902-4e56-8fff-3bfc182268c7
λ_13_U, λ_18_U = (2.60076, 2.73079) .* U.u"mm"

# ╔═╡ 321fa5f5-be6a-4e5f-ba7f-4c117a240253
λ_13_DQ |> DQ.us"Hz"

# ╔═╡ 2962b22a-5e4a-4dc1-9be6-8e57979ecba5
λ_13_U |> U.u"Hz"

# ╔═╡ ee46b865-dbd4-4939-b241-514941dd138d
ν_13_DQ, ν_18_DQ = C_DQ.c ./ (λ_13_DQ, λ_18_DQ) .|> DQ.us"Hz"

# ╔═╡ 554321c7-cdcb-415c-8ac2-c63b42b6889a
ν_13_U, ν_18_U = U.uconvert.(U.u"Hz", (λ_13_U, λ_18_U), UE.Spectral())

# ╔═╡ 0af912a5-18ad-4974-a0b5-c1f66c8fa37a
md"""
!!! note
    DynamicQuantities.jl does not have similar unit eqivalence functionality at this time because it is cheap to just roll our own conversions since the physical constants needed are already included in the package. For this reason, we just compute the ``c / λ`` relation directly.
"""

# ╔═╡ 3a24f4aa-b074-4beb-b77d-6778e2fe580a
CC_DQ = let
    A_13, A_18 = (7.4e-8, 8.8e-8) ./ DQ.us"s"
    3.0e14 * DQ.us"s/(K*cm^2*km)" * (ν_18_DQ / ν_13_DQ)^3 * (A_13 / A_18)
end

# ╔═╡ 98869c25-2644-47d1-b8cc-05699292f2a8
CC_U = let
    A_13, A_18 = (7.4e-8, 8.8e-8) ./ U.u"s"
    3.0e14 * U.u"s/(K*cm^2*km)" * (ν_18_U / ν_13_U)^3 * (A_13 / A_18)
end

# ╔═╡ 4092a893-818e-49f4-93d0-be7bd652dddc
B_DQ = C_DQ.h * ν_18_DQ / (C_DQ.k_B * T_ex_DQ)

# ╔═╡ a0967b01-4eb8-4160-af68-a908734dc25e
B_U = C_U.h * ν_18_U / (C_U.k_B * T_ex_U)

# ╔═╡ 82f11626-722b-46fa-9173-8b5d7a80a190
B_DQ

# ╔═╡ 392d3055-2ce0-4db7-8702-f79395f24b3b
B_U |> U.NoUnits

# ╔═╡ ef498536-8fb8-46e3-9d8d-f7eb3994a3f5
NCO_DQ = CC_DQ * intcloud_DQ / (1 - exp(-B_DQ))

# ╔═╡ 84f948ad-f903-46a6-a520-2499864f348f
NCO_U = CC_U * intcloud_U / (1 - exp(-B_U))

# ╔═╡ 25598fcd-761d-40b8-95ac-8b68a00026da
md"""
!!! note ""
    **Peak CO Column density (DQ):** $(maximum(NCO_DQ))
"""

# ╔═╡ 336f8626-7b0c-459d-bedb-281d010fcbd2
md"""
!!! note ""
    **Peak CO Column density (U):** $(maximum(NCO_U))
"""

# ╔═╡ eba2f06f-ebe9-492d-81d2-1cc4fccd5b0a
md"""
### ``\mathrm{CO}`` to Total Mass

We are using ``\mathrm{CO}`` as a tracer for the much more numerous ``\mathrm{H}_2``, the quantity we are actually trying to infer. Using the (known/assumed) ``\mathrm{H}_2 / \mathrm{CO}`` ratio:
"""

# ╔═╡ de22c778-1529-4177-85e6-a0178f437a8c
H₂_CO_ratio = 5.9e6

# ╔═╡ c262aa56-b9b6-4da4-8810-d6120f0724c6
NH₂_DQ = NCO_DQ * H₂_CO_ratio

# ╔═╡ 9a98c6bb-d0d5-426e-9137-9fc94fc944be
NH₂_U = NCO_U * H₂_CO_ratio

# ╔═╡ 1e8ec5c5-3d41-44a7-8dfa-81d981750d9e
md"""
!!! note ""
    **Peak ``\mathrm{H}_2`` column density (DQ):** $(maximum(NH₂_DQ))
"""

# ╔═╡ dd22bc93-28bf-412e-a52b-d4332475daa2
md"""
!!! note ""
    **Peak ``\mathrm{H}_2`` column density (U):** $(maximum(NH₂_U))
"""

# ╔═╡ 522d88d6-5f60-401b-8786-0236c0859eda
ρ_DQ = let
    mH₂ = 2 * 1.008 * C_DQ.u
    NH₂_DQ * mH₂
end

# ╔═╡ 014b94b3-a6fa-4627-a24d-fe8c5e275582
ρ_U = let
    mH₂ = 2 * 1.008 * U.u"u"
    NH₂_U * mH₂
end

# ╔═╡ 0eefe540-1669-4510-bef7-8fb7f8be68f1
Δap_DQ = Δra_DQ * Δdec_DQ

# ╔═╡ 085c04f2-0762-4213-b874-14f3ce1f0b49
Δap_U = Δra_U * Δdec_U

# ╔═╡ b057812d-1a8e-479c-8192-1c1acdb97d23
Δa_DQ = Δap_DQ * d_DQ^2

# ╔═╡ 98f0753f-635e-46e3-8045-ff40ca2e827b
Δa_U = Δap_U * d_U^2

# ╔═╡ 479b900b-b1bf-4d89-a53e-b1b20da7f5ad
Δa_DQ |> C_DQ.pc^2

# ╔═╡ 28237d12-6065-4fa3-bd40-6397dae7ecb6
Δa_U |> U.u"pc^2"

# ╔═╡ 9c0862d0-8028-4ef4-a867-18c7932dad91
md"""
!!! tip
    To treat angles as physical units instead of as dimensionless, see [this relevant section](https://juliaphysics.github.io/DynamicQuantities.jl/stable/examples/#3.-Using-dimensional-angles) in the DynamicQuantities.jl documentation, and this extension package for Unitful.jl: [DimensionfulAngles.jl](https://github.com/cmichelenstrofer/DimensionfulAngles.jl).
"""

# ╔═╡ 88966b7e-b261-45a6-9f12-6b74f24c03e4
M_cloud_DQ = sum(ρ_DQ * Δa_DQ) |> C_DQ.M_sun

# ╔═╡ 1aaaebb6-839a-4598-82b9-dc17efdf7623
M_cloud_U = sum(ρ_U * Δa_U) |> U.u"Msun"

# ╔═╡ 802f5cad-a0c0-4426-94f2-426f89dea7e1
md"""
!!! tip "Exercises"
    The astro material was pretty heavy on that one, so let's focus on some associated statistics using DynamicQuantities.jl's array capabililities. Compute the median and mean of the data with the `mean` and `median` functions. Why are their values so different?

    Similarly, compute the standard deviation and variance. Do they have the units you expect?
"""

# ╔═╡ 2eed0bc8-2306-42bb-9803-33ae131021e6
md"""
## 3. Using units with functions

Units are also a useful tool if you plan to share some of your code, either with collaborators or the wider community. By writing functions that take numbers with units objects instead of raw numbers or arrays, you can write code that is agnostic to the input unit. In this way, you may even be able to prevent [the destruction of Mars orbiters](http://en.wikipedia.org/wiki/Mars_Climate_Orbiter#Cause_of_failure). Below, we provide a simple example.

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

# ╔═╡ 1c3faf9f-adcf-4232-acf4-655d828623e3
function response_func_good_DQ(x, y)
    xscale = 0.9 * DQ.us"arcsec" # We use symbolic dimensions here
    yscale = 0.85 * DQ.us"arcsec" # to treat radians as unitful quantities
    xfactor = 1 / (1 + x / xscale)
    yfactor = 1 / (1 + y / yscale)

    return xfactor * yfactor
end

# ╔═╡ 6708c2f5-a291-486e-8e66-93acfcedbcb7
function response_func_good_U(x, y)
    xscale = 0.9 * U.u"arcsecondᵃ" # We use dimensionful angles here
    yscale = 0.85 * U.u"arcsecondᵃ" # to treat radians as unitful quantities
    xfactor = 1 / (1 + x / xscale)
    yfactor = 1 / (1 + y / yscale)

    return xfactor * yfactor
end

# ╔═╡ ce1804da-df61-4e4c-aa3b-d950e83b7c13
response_func_good_DQ(1.0, 1.2)

# ╔═╡ 910920d2-085a-498f-8458-e6e96903ab92
response_func_good_U(1.0, 1.2)

# ╔═╡ f7486c08-b74b-471f-bbc2-690439a487f8
response_func_good_DQ(1.0 * DQ.u"arcmin", 1.2 * DQ.u"arcmin")

# ╔═╡ 08c5e6fc-d8b5-4311-8f35-449a830daeec
response_func_good_U(1.0 * U.u"arcminuteᵃ", 1.2 * U.u"arcminuteᵃ")

# ╔═╡ c143ab62-b3f0-4b6a-83aa-76e5dc3548ac
md"""
The funding agency is impressed at the resolution you achieved, and your instrument is saved! You now go on to win the Nobel Prize due to discoveries the instrument makes. And it was all because you used `Quantity` as the input of code you shared.
"""

# ╔═╡ 0c05c385-67bd-4078-8fe2-34a37a13b312
md"""
!!! note
    **DynamicQuantities.jl:** Note the use of `us""` instead of `u""`, both supplied by DynamicQuantities.jl, to achieve this behavior.

    **Unitful.jl:** Note the use of `ᵃ` supplied by DimensionfulAngles.jl instead of the base angle units provided by Unitful.jl to achieve this behavior.
"""

# ╔═╡ 7d4f9a02-81a7-4bdc-a22d-3435192f9f15
let
    sol = details(
        "Example solution",
        md"""
        ```julia
        v_orb(M, r) = sqrt(C_DQ.G * M / r)
        ```
        """
    )

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
TableOfContents(; depth = 4)

# ╔═╡ e0d2d6c4-d363-4bb2-9d12-42c4a52aba3b
function side_by_side(content = nothing)
    return @htl(
        """
        <style>
        @media (min-width: 900px) {
            pluto-cell:has(> pluto-output .sbs-marker) + pluto-cell {
                margin-right: 1rem;
            }
            pluto-cell:has(> pluto-output .sbs-marker) + pluto-cell,
            pluto-cell:has(> pluto-output .sbs-marker) + pluto-cell + pluto-cell {
                display: inline-block;
                vertical-align: top;
                width: calc((100% - 1rem) / 2);
                box-sizing: border-box;
            }
        }
        </style>
        <div class="sbs-marker" style="display:none"></div>
        $(content)
        """
    )
end

# ╔═╡ 60c89d86-942e-4c97-bd7a-ad2f792b1155
md"""
## 1. Galaxy mass

In this first example, we will use objects with units to estimate a hypothetical galaxy's mass, given its half-light radius and radial velocities of stars in the galaxy.

Let's assume that we measured the half-light radius of the galaxy to be 29 pc projected on the sky at the distance of the galaxy. This radius is often called the "effective radius", so we'll store it with the name `Reff`. The easiest way to create an object with units is by multiplying the value with its unit:
""" |> side_by_side

# ╔═╡ a9d14d8e-aefa-48d9-a3b3-ca13e284ab08
md"""
We can access the value and unit of quantities using the `ustrip()` and `dimension()`/`unit()` functions:
""" |> side_by_side

# ╔═╡ cc736aec-c48a-411a-af83-4255309d77e9
md"""
Furthermore, we can convert the radius to any other unit of length using the `uconvert()` function. Here, we convert it to kilometers:
""" |> side_by_side

# ╔═╡ b7f7e8e2-7f26-4d27-8dc5-8204761a9965
md"""
Units are also "callable," meaning we can treat them like functions. Even though the first example with calling a unit is a bit awkward, the `|>` syntax makes it much more natural:
""" |> side_by_side

# ╔═╡ eb31da4c-7011-49ec-a9ac-194eaee10c7a
md"""using the `Normal` distribution from the Distributions.jl package:
""" |> side_by_side

# ╔═╡ f7ceb7c4-6488-43eb-b211-9f7087762a56
md"""
which we can then visualize with Makie.jl using its automatic unit support:
""" |> side_by_side

# ╔═╡ ce5aea87-2614-4e72-8e88-25a8c590b89f
md"""
Now we can calculate the velocity dispersion of the galaxy. This demonstrates how you can perform basic operations like subtraction and division with objects with units, and also use them in standard Julia functions such as `mean()` and `size()`. They retain their units through these operations just as you would expect them to:
""" |> side_by_side

# ╔═╡ bb9d675b-00ce-4d20-b78b-a08835b47137
md"""
Now for the actual mass calculation. If a galaxy is pressure-supported (for example, an elliptical or dwarf spheroidal galaxy), its mass within the stellar extent can be estimated using a straightforward formula: ``M_{1/2} = 4σ^2 R_\text{eff}/G``. There are caveats to the use of this formula for science -- see [Wolf et al. 2010](http://ui.adsabs.harvard.edu/abs/2010MNRAS.406.1220W/abstract) for details.

!!! note
    Constants [are included in DQ](https://juliaphysics.github.io/DynamicQuantities.jl/stable/constants/):

    ```julia
    using DynamicQuantities: Constants as C_DQ
    ```

    For demonstration purposes, we will actually be importing from `SymbolicConstants` instead. In practice, it is much more performant to work in the base SI system provided by `Constants`, and just convert to whatever units are desired in the end.

    For Unitful.jl, a separate package like [PhysicalConstants.jl for Unitful.jl](https://juliaphysics.github.io/PhysicalConstants.jl/stable/) is required:

    ```julia
    using PhysicalConstants: CODATA2018 as C_U
    ```
""" |> side_by_side

# ╔═╡ 0ab125c6-568b-414c-acaa-abcfebb559e2
md"""
We can also easily express the mass in whatever form we like -- solar masses are common in astronomy, or maybe you want the default SI and CGS units:
""" |> side_by_side

# ╔═╡ 5f438b76-7fc7-49ed-8a0a-7f2d5a760414
md"""
Or, if you want the log of the mass, you can just use the builtin `log10` as long as the logarithm's argument is dimensionless.
""" |> side_by_side

# ╔═╡ 5a9dfdf7-aaf8-421d-a451-6a6a5c35e1cc
md"""
Note that this is different than:
""" |> side_by_side

# ╔═╡ a3f42e7e-af38-46b0-b1fb-bdbf77312e16
md"""
emphasizing the importance of being explicit with our units. Similarly, taking the logarithm of something with units is not mathematically well defined, so this will sensibly error as well:
""" |> side_by_side

# ╔═╡ 78740b86-45e9-45f5-b4c2-be6cb65f1368
md"""
### Setting up the data cube

Let's assume that we've mapped the inner part of a molecular cloud in the ``J = 1 - 0`` rotational transition of ``\mathrm{C}^{18}\mathrm{O}`` and are interested in measuring its total mass. The measurement produced a data cube with RA and Dec as spatial coordinates and velocity as the third axis. Each voxel in this data cube represents the brightness temperature of the emission at that position and velocity. Furthermore, we'll assume that we have an independent measurement of distance to the cloud ``d = 250\, \mathrm{pc}`` and that the excitation temperature is known and constant throughout the cloud: ``T_\text{ex} = 25\, \mathrm{K}``:
""" |> side_by_side

# ╔═╡ b86a31df-9de1-4e27-8ca7-ce7f3d2576fa
md"""
We'll generate a synthetic dataset, assuming the cloud follows a Gaussian distribution in each of RA, Dec, and velocity. We start by creating a 100×100×300 array, such that the first coordinate is right ascension, the second is declination, and the third is velocity. In this data cube, the cloud is positioned at the center, with ``\sigma`` and the center in each dimension shown below:
""" |> side_by_side

# ╔═╡ e13cbda5-ddd4-491e-bf2a-d6ed24d2bdbd
md"""
Note in particular that the ``\sigma`` for RA and Dec have different units from the center, but we are able to automatically handle the relevant conversions before computing the exponential:
""" |> side_by_side

# ╔═╡ 00b9d37f-ce7b-4491-b1df-f0963f2598a8
md"""
We will also need to know the size of each pixel:
""" |> side_by_side

# ╔═╡ ccb58b55-b9c8-4229-aff8-92541af123ec
md"""
and the width of each velocity bin:
""" |> side_by_side

# ╔═╡ 7ab2fe80-812f-4a21-933a-e169a17f6c32
md"""
We're interested in the integrated intensity over all of the velocity channels, so let's create a 2D quantity array by summing our data cube along the velocity axis (multiplying by the velocity width of a pixel):
""" |> side_by_side

# ╔═╡ 4ee79c7c-5c4e-457b-b713-103937e50355
md"""
We can plot the 2D quantity using Makie's `heatmap` function:

!!! todo
    Heatmap + colorbar + units support coming in <https://github.com/MakieOrg/Makie.jl/pull/5673> 
""" |> side_by_side

# ╔═╡ 3a3328d5-31de-4deb-bff7-d25c1fcbc4ef
md"""
Here we have given an expression for ``C`` scaled to the values for ``\mathrm{C}^{13}\mathrm{O}`` (``\nu_{13}`` and ``A_{13}``). In order to use this relation for ``\mathrm{C}^{18}\mathrm{O}``, we need to rescale the frequencies ``\nu`` and the Einstein coefficients (``A``). Lastly, ``C`` is in funny mixed units, but that's okay. We'll be able to do our unit handling in the usual way.

First, we look up the wavelength for these emission lines and store them as quantities:
""" |> side_by_side

# ╔═╡ e0a2c745-41cf-4359-bb99-3117fbb507cc
md"""
Since the wavelength and frequency of light are related using the speed of light, we can convert between them. However, doing so just using the `uconvert()` function or its equivalent `|>` fails, as units of length and frequency are not directly convertible:
""" |> side_by_side

# ╔═╡ 7fa76f60-c623-4b25-b185-0cd5e805ef71
md"""
Fortunately, the Unitful.jl ecosystem comes to the rescue by providing a feature called "unit equivalences." Equivalences provide a way to convert between two physically different units that are not normally equivalent, but in a certain context have a one-to-one mapping. For more on equivalencies, see [the documentation for UnitfulEquivalences.jl](https://sostock.github.io/UnitfulEquivalences.jl) (UE).

In this case, passing the `Spectral()` argument to `uconvert()` provides the equivalences necessary to handle conversions between wavelength and frequency:
""" |> side_by_side

# ╔═╡ d0215082-52df-4eae-a409-7c3c7cfc69ed
md"""
Next, we look up Einstein coefficients (in units of s⁻¹), and calculate the ratios in constant ``CC``. Note how the ratios of frequency and Einstein coefficient units are dimensionless, so the unit of ``CC`` is unchanged:
""" |> side_by_side

# ╔═╡ d9e2b828-9c32-4638-bfcd-0c16b221aa43
md"""
Now we move on to calculate the constant ``B``. This is given by the ratio of ``\dfrac{hν}{k_\text{B}T}``, where ``h`` is Planck's constant, ``k_\text{B}`` is the Boltzmann's constant, ``ν`` is the emission frequency, and ``T`` is the excitation temperature. The constants were imported from `DynamicQuantities.Constants`, and the other two values are already calculated, so here we just take the ratio:
""" |> side_by_side

# ╔═╡ a0b64ba9-cbec-404e-bf39-aee03ae407ae
md"""
!!! note
    Note how DynamicQuantities.jl intelligently cancelled the units for us, while still keeping this as a `Quantity` object. For Unitful.jl, a conversion to `NoUnits` is needed. In practice, it is generally prefable to not worry about these intermediate units and just convert everything at the end.
""" |> side_by_side

# ╔═╡ 6aff587f-be6f-4fd8-96df-9c12f3769f32
md"""
At this point we have all the ingredients to calculate the number density of ``\mathrm{CO}`` molecules in this cloud. We already integrated (summed) over the velocity channels above to show the integrated intensity map, but we'll do it again here for clarity. This gives us the column density of ``\mathrm{CO}`` for each spatial pixel in our map:
""" |> side_by_side

# ╔═╡ 63d629a1-379f-44bd-a02d-d252d0594357
md"""
We can then print out the peak column column density:
""" |> side_by_side

# ╔═╡ 072535ef-9a0d-4c1e-aea7-6486d329f66f
md"""
we can calculate the ``\text{H}_2`` column density by multiplying the ``\mathrm{CO}`` column density found previously with this ratio:
""" |> side_by_side

# ╔═╡ 00f69c40-f31c-49e6-a1cc-d84a4cb2b43d
md"""
The peak column column density is then:
""" |> side_by_side

# ╔═╡ c6b35992-973b-49d4-bfa7-855e8fe10924
md"""
That's a peak column density of roughly 50 magnitudes of visual extinction (assuming the conversion between ``N_{\mathrm{H}_2}`` and ``A_V`` from [Bohlin et al. 1978](http://ui.adsabs.harvard.edu/abs/1978ApJ...224..132B/abstract)), which seems reasonable for a molecular cloud.

We obtain the mass column density by multiplying the number column density by the mass of an individual ``\mathrm{H_2}`` molecule:
""" |> side_by_side

# ╔═╡ fbd0d75b-b3af-4eb5-a249-423317d656e5
md"""
A final step in going from the column density to mass is summing up over the area. If we do this in the straightforward way of length × width of a pixel, this area is then in units of deg²:
""" |> side_by_side

# ╔═╡ a427567e-b6dd-44a9-a61f-a725960cc459
md"""
In the small angle approximation, multiplying the pixel area with the square of distance yields the cross-sectional area of the cloud that the pixel covers, in physical units:
""" |> side_by_side

# ╔═╡ cefdeee7-6a4b-4c60-a14c-c1edac53be98
md"""
Angles in both packages are dimensionless by default, so the above quantities can be freely converted to other compatible dimensions:
""" |> side_by_side

# ╔═╡ e81f3b44-bc04-4768-81c2-c2733958f5a9
md"""
Finally, multiplying the column density with the pixel area and summing over all the pixels gives us the cloud mass:
""" |> side_by_side

# ╔═╡ 71cca345-cffc-4aee-a671-08901a92babe
md"""
And now they tell all their friends how terrible the instrument is, because it's supposed to have arcsecond resolution, but your function clearly shows it can only resolve an arcmin at best. But you can solve this by requiring they pass in Quantity objects. The new function could simply be:
""" |> side_by_side

# ╔═╡ 7cfb17fc-81a4-4522-a88b-2fd10001baaa
md"""
And your collaborator now has to pay attention. If they just blindly put in a number, they get an error:
""" |> side_by_side

# ╔═╡ 5590c909-6103-4427-b1a2-4cf35265dc07
md"""
which is their cue to provide the units explicitly:
""" |> side_by_side

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
ComputePipeline = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
DimensionfulAngles = "2d4b8d7a-02d9-40f9-9abe-9c695b77de0d"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
PhysicalConstants = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAstro = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
UnitfulEquivalences = "da9c4bc3-91c8-4f02-8a40-6b990d2a7e0c"

[sources]
CairoMakie = {rev = "ff/color-dim-converts", subdir = "CairoMakie", url = "https://github.com/MakieOrg/Makie.jl"}
ComputePipeline = {rev = "ff/color-dim-converts", subdir = "ComputePipeline", url = "https://github.com/MakieOrg/Makie.jl"}
DynamicQuantities = {rev = "main", url = "https://github.com/JuliaPhysics/DynamicQuantities.jl"}
Makie = {rev = "ff/color-dim-converts", subdir = "Makie", url = "https://github.com/MakieOrg/Makie.jl"}

[compat]
CairoMakie = "~0.16.0"
ComputePipeline = "~0.2.0"
DimensionfulAngles = "~0.2.2"
Distributions = "~0.25.132"
DynamicQuantities = "~1.13.0"
HypertextLiteral = "~1.0.0"
Makie = "~0.25.0"
PhysicalConstants = "~0.2.5"
PlutoUI = "~0.7.83"
StatsBase = "~0.34.13"
Unitful = "~1.29.0"
UnitfulAstro = "~1.2.2"
UnitfulEquivalences = "~0.2.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.2"
manifest_format = "2.1"
project_hash = "ee62d42e02852c5640547f826fb7997f922e8208"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
registries = "General"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
registries = "General"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.45"

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
deps = ["LinearAlgebra"]
git-tree-sha1 = "dee8388a9bec89c8713ad176f96960ae2f98598a"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.4"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptiveArrayPools]]
deps = ["Preferences", "Printf", "Random"]
git-tree-sha1 = "5feaa5055572e809b69b62f06826b1c1c8d3ebf3"
registries = "General"
uuid = "4f381ef7-9af0-4cbe-99d4-cf36d7b0f233"
version = "0.3.11"

    [deps.AdaptiveArrayPools.extensions]
    AdaptiveArrayPoolsCUDAExt = "CUDA"
    AdaptiveArrayPoolsMetalExt = "Metal"

    [deps.AdaptiveArrayPools.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
registries = "General"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
registries = "General"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "94eab0b3ccdcac361188cc661daf69d4433c1818"
registries = "General"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.2.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
registries = "General"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
registries = "General"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "8c290a1b223deaeea9aea44b235d24546da8eb98"
registries = "General"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.4.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
registries = "General"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
registries = "General"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
registries = "General"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
registries = "General"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
registries = "General"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "c815ac9dd39d95b225b1e00ad13eb5dc823402ee"
repo-rev = "ff/color-dim-converts"
repo-subdir = "CairoMakie"
repo-url = "https://github.com/MakieOrg/Makie.jl"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.16.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "0cfb779e24fdced77a05df1dbe99200821a9f563"
registries = "General"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.8+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
registries = "General"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ChainRulesCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "4378c3544d8fcc3c8a521782db0e93e2c32d1838"
registries = "General"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.4"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "da54a6cd93c54950c15adf1d336cfd7d71f51a56"
registries = "General"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.7"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
registries = "General"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
registries = "General"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "a7cf827a79b378a381e31e5c7feca853ec19cb94"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.4"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
registries = "General"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "a633156e74096f90b3b5e833529292b651b5671a"
registries = "General"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.3"

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
registries = "General"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
registries = "General"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "a2a82d239c8f7132c25a810a85a502e0ddbaa8cf"
repo-rev = "ff/color-dim-converts"
repo-subdir = "ComputePipeline"
repo-url = "https://github.com/MakieOrg/Makie.jl"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.2.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
registries = "General"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
registries = "General"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
registries = "General"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
registries = "General"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
registries = "General"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
registries = "General"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "e1e28f37cbdca4b594ad86ecfe2ee8f25796ebff"
registries = "General"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.7.0"

[[deps.DimensionfulAngles]]
deps = ["Unitful", "UnitfulAngles", "UnitfulEquivalences"]
git-tree-sha1 = "bbf9ded2c249cd3b08bd0053f4acc40294c9ef64"
registries = "General"
uuid = "2d4b8d7a-02d9-40f9-9abe-9c695b77de0d"
version = "0.2.2"

[[deps.DispatchDoctor]]
deps = ["MacroTools", "Preferences"]
git-tree-sha1 = "278d0120d7c22ef36bd8d22f349868d6a782a949"
registries = "General"
uuid = "8d63f2c5-f18a-4cf2-ba9d-b3f60fc568c8"
version = "0.4.29"

    [deps.DispatchDoctor.extensions]
    DispatchDoctorChainRulesCoreExt = "ChainRulesCore"
    DispatchDoctorEnzymeCoreExt = "EnzymeCore"

    [deps.DispatchDoctor.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "Roots", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "2ddbc341251ed713d39e503b504593a5d6789551"
registries = "General"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.132"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.DynamicQuantities]]
deps = ["DispatchDoctor", "PrecompileTools", "TestItems", "Tricks"]
git-tree-sha1 = "1bb3a30dcdea0e3f8497c18d30b81ec71bb2790f"
repo-rev = "main"
repo-url = "https://github.com/JuliaPhysics/DynamicQuantities.jl"
uuid = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
version = "1.13.0"

    [deps.DynamicQuantities.extensions]
    DynamicQuantitiesLinearAlgebraExt = "LinearAlgebra"
    DynamicQuantitiesMeasurementsExt = "Measurements"
    DynamicQuantitiesRecursiveArrayToolsExt = "RecursiveArrayTools"
    DynamicQuantitiesSciMLBaseExt = "SciMLBase"
    DynamicQuantitiesScientificTypesExt = "ScientificTypes"
    DynamicQuantitiesUnitfulExt = "Unitful"

    [deps.DynamicQuantities.weakdeps]
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"
    SciMLBase = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
    ScientificTypes = "321657f4-b219-11e9-178b-2701a2544e81"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
registries = "General"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
registries = "General"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
registries = "General"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0f666deaf8099a7d7ac229a74402ada88898be1"
registries = "General"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.9.0+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
registries = "General"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FastInterpolations]]
deps = ["AdaptiveArrayPools", "HelpPlots", "LinearAlgebra", "Preferences", "Printf"]
git-tree-sha1 = "f7bcd496f0a2e4e1174e7e1f667f27fb671fd7df"
registries = "General"
uuid = "9ea80cae-fc13-4c00-8066-6eaedb12f34b"
version = "0.4.20"

    [deps.FastInterpolations.extensions]
    FastInterpolationsChainRulesCoreExt = "ChainRulesCore"
    FastInterpolationsColorVectorSpaceExt = ["ColorTypes", "ColorVectorSpace"]
    FastInterpolationsEnzymeExt = "Enzyme"
    FastInterpolationsForwardDiffExt = "ForwardDiff"
    FastInterpolationsRecipesBaseExt = "RecipesBase"
    FastInterpolationsSymbolicsExt = ["Symbolics", "SymbolicUtils"]

    [deps.FastInterpolations.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ColorTypes = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
    ColorVectorSpace = "c3611d14-8923-5661-9e6a-0046d554d3a4"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6621fef488e496356c9c9625d0562c12a6070819"
registries = "General"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.20.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
registries = "General"
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
registries = "General"
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
git-tree-sha1 = "417289ac840d201005d9ed16789a46defccbe36d"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.18.1"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "e1e9ba3a8f5d71e6ae698f138a1362aeabb42e9f"
registries = "General"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.18.3+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
registries = "General"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d95a7a3a7e0f2fffc6291bd33035830aab01310b"
registries = "General"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.15.0+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
registries = "General"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f3f50c8001bdfa84f312234fa86715b7c9243771"
registries = "General"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.18+0"

[[deps.Gamma]]
deps = ["LogExpFunctions"]
git-tree-sha1 = "becc397f7cfb06e343496ae6ffb04818a851da51"
registries = "General"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.2.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "ec46c5825710fa1a15d468acb2d93cc939a7a5fe"
registries = "General"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.13"

    [deps.GeometryBasics.extensions]
    ExtentsExt = "Extents"
    GeometryBasicsGeoInterfaceExt = "GeoInterface"
    IntervalSetsExt = "IntervalSets"

    [deps.GeometryBasics.weakdeps]
    Extents = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
registries = "General"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
registries = "General"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3efbbc027441271444dcd0c0a46f2d119dc4329"
registries = "General"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "6.1.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "ea9f1fb37ff10597353b56f3ba0ebd525acef6e0"
registries = "General"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.90.1+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
registries = "General"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
registries = "General"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "ef70da5e123a06a29e2d6ddff0f09985bc226491"
registries = "General"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "9e544f0758a2341e7e782fde6f3f67610366661c"
registries = "General"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "100.14006.0+0"

[[deps.HelpPlots]]
git-tree-sha1 = "07d389d81ca31e42d68bd88cfa1e900a63b7d378"
registries = "General"
uuid = "5880b2e8-ad0a-48c9-bacd-30dc9c1357e0"
version = "1.1.0"

    [deps.HelpPlots.extensions]
    HelpPlotsExt = "Plots"

    [deps.HelpPlots.weakdeps]
    Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[[deps.HypergeometricFunctions]]
deps = ["Gamma", "LinearAlgebra"]
git-tree-sha1 = "31bb6c92405c084617facc1d7ed9eb6c402d061e"
registries = "General"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.30"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
registries = "General"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
registries = "General"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
registries = "General"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "f0f005f997dfb8c5fe23920d99458a9619873893"
registries = "General"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.10"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
registries = "General"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
registries = "General"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
registries = "General"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
registries = "General"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "c72458f1962faeb003bf23cbdb75164fe6280906"
registries = "General"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "48922d06068130f87e43edef52382e6a94305ae6"
registries = "General"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.3"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "ff294afb9a15d31d8d7422da138844641a73135f"
registries = "General"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.11"

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
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

    [deps.IntervalSets.weakdeps]
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.InverseFunctions]]
git-tree-sha1 = "5cc6ad3e2bfc75a11eed92632ee341d37f4a4536"
registries = "General"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.18"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "d22e8bdb3f03c300a103ff23c0dd5a654c6701de"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.7"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
registries = "General"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
registries = "General"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
registries = "General"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "9d8bead1c9ed9f2dc4eb5b51947262163b38ffce"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.12.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
registries = "General"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "037babc10853eeb8e585418922246cb97b8e5b74"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+1"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "9eda8292dd3268b3b7ec9df21bbfac24e177ec52"
registries = "General"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.12"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "39bca05343661c347aae0bca57a5994a0bf4f08d"
registries = "General"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.2.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5b100780d4d30d63b4618d7930d48af409c1772"
registries = "General"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "23.1.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "df7566479bd64f20bd16b09960145e70160ffb3b"
registries = "General"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.12"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
registries = "General"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.22.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.7+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
registries = "General"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
registries = "General"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
registries = "General"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b8a91037379e1fd661e9b5ff3d7a770424b1f9f2"
registries = "General"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.4+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
registries = "General"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9e887184f306b28dbfcff88cefbb0ec80f1aa9e0"
registries = "General"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.4+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "28d620bcdb37a3fb02de9f521928ee453690f414"
registries = "General"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.3"

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
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FastInterpolations", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "Latexify", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "69b2d8f8064365e695f7a9bcb60bceafad2eee7a"
repo-rev = "ff/color-dim-converts"
repo-subdir = "Makie"
repo-url = "https://github.com/MakieOrg/Makie.jl"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.25.0"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"
    MakieFFMPEGExt = "FFMPEG_jll"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
    FFMPEG_jll = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
registries = "General"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "aa1078778be5a8e5259ff04fbc3d258b3e78d464"
registries = "General"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.9"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "Requires"]
git-tree-sha1 = "bdcde8ec04ca84aef5b124a17684bf3b302de00e"
registries = "General"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.11.0"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
registries = "General"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.MuladdMacro]]
deps = ["PrecompileTools"]
git-tree-sha1 = "283bf85d4a767481dd924dff0eee1735e95f449e"
registries = "General"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.7"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "676c0f1c1c34224d6869675047a54b42fdc7b8c5"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.5"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
registries = "General"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
registries = "General"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
registries = "General"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "38a93f17e431141c6470bb67a88952a7c4f0e928"
registries = "General"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.34+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
registries = "General"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "OpenJPH_jll", "Zlib_jll", "Zstd_jll", "libdeflate_jll"]
git-tree-sha1 = "eb47260eaf61dca5c0e17bd126f4426edac35975"
registries = "General"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.5.2+0"

[[deps.OpenJPH_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "1b1f3cabb1130670812b848f4623d69d578ea805"
registries = "General"
uuid = "95196869-6d61-5868-8eee-c38c33804e0c"
version = "0.32.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.8+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
registries = "General"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "b079a6affebde167dad8ffe0d670b53917fbcb68"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "123266c25174ef6c8d4718920abc206452cf8de6"
registries = "General"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.41"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "32b657a0d57c310a1a172bfc8c8cf68c5e674323"
registries = "General"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.5"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
registries = "General"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
registries = "General"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1912a9f1b9ca55005b03ba075f8e19993583e237"
registries = "General"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.58.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools"]
git-tree-sha1 = "663e8b48b789916221e0765393b289ca6c88f24e"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "3.0.0"

[[deps.PhysicalConstants]]
deps = ["Measurements", "Roots", "Unitful"]
git-tree-sha1 = "3260db49a2599ed83d649b2c3043f1733c3d2895"
registries = "General"
uuid = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
version = "0.2.5"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
registries = "General"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
registries = "General"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Reexport", "Statistics"]
git-tree-sha1 = "ddcc36fd83fb4cd459a53c891981697bde0d3de6"
registries = "General"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.5.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
registries = "General"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "aef1987d90dcc2e2cbd919ad10b2d25bcce4f418"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
registries = "General"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Printf"]
git-tree-sha1 = "ce5ba3edfdeba107a721748d53b5f5887f415853"
registries = "General"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.13.0"

    [deps.ProgressMeter.extensions]
    ProgressMeterDistributedExt = "Distributed"
    ProgressMeterIJuliaExt = "IJulia"

    [deps.ProgressMeter.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
registries = "General"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
registries = "General"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
registries = "General"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
registries = "General"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
registries = "General"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
registries = "General"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d40b2fe70437b01397d2a4d5b020008da4e7019"
registries = "General"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.2+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "e6ffa7bd4f8d0b8744d18f219304738d952f9a37"
registries = "General"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.11"

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
registries = "General"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "37ac518b3c924f110d82d50173eaf4ebeaf27ad7"
registries = "General"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.8.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
registries = "General"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "57aa595158717ef165e6f5ab639fe2e3178c0a2b"
registries = "General"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
registries = "General"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
registries = "General"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2c79185e5261b159474903598c571ddd9f12ff7b"
registries = "General"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.6"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
registries = "General"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "88ef0d98204515bc53dd4aaea5f7788882d1b66c"
registries = "General"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.9.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
registries = "General"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "8a8a52f887f38d44014e24ee46b8a7311c625f6d"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.24"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
registries = "General"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
registries = "General"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "adb9da019510162e67a4493fc235c23203d8b09e"
registries = "General"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.13"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91a5737baed20ee31f3faea0e51f57461f6a689e"
registries = "General"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.2.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "4daaddbe37aa80b4551bf7c56adf53e7cfcdf35a"
registries = "General"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.4"

    [deps.StructArrays.extensions]
    AdaptExt = "Adapt"
    GPUArraysCoreExt = "GPUArraysCore"
    KernelAbstractionsExt = "KernelAbstractions"
    LinearAlgebraExt = "LinearAlgebra"
    SparseArraysExt = "SparseArrays"
    StaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "044ed182312a8ad0f28437f048f3e1eae9093a80"
registries = "General"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.11.1"

    [deps.StructUtils.extensions]
    StructUtilsDataFramesExt = ["DataFrames"]
    StructUtilsLazilyInitializedFieldsExt = ["LazilyInitializedFields"]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazilyInitializedFields = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
registries = "General"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
registries = "General"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
registries = "General"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TestItems]]
git-tree-sha1 = "a6dd904babc04670c784f81b67957a3545271610"
registries = "General"
uuid = "1c621080-faea-4a02-84b6-bbd5e436b8fe"
version = "1.1.0"

[[deps.TiffImages]]
deps = ["CodecZstd", "ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "9ca5f1f2d42f80df4b8c9f6ab5a64f438bbd9976"
registries = "General"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.9"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
registries = "General"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
registries = "General"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
deps = ["Serialization"]
git-tree-sha1 = "50eda8011f098869a76da69a524ff59e0298711a"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.2"

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
registries = "General"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "1f0f9f401753701a7e4113b5056ca38d33875b55"
registries = "General"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.29.0"

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
registries = "General"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.7.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "fbe44a0ade62ae5ed0240ad314dfdd5482b90b40"
registries = "General"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.2"

[[deps.UnitfulEquivalences]]
deps = ["Unitful"]
git-tree-sha1 = "249cafa9e88c0c28dbd7b55e258a53a01761f974"
registries = "General"
uuid = "da9c4bc3-91c8-4f02-8a40-6b990d2a7e0c"
version = "0.2.1"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
registries = "General"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
registries = "General"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e52eca002a11c30a858185efdfb15311e1c7a6bf"
registries = "General"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "ce91b54a5f90c4f27df6950706234b34e49cf79e"
registries = "General"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.14+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
registries = "General"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
registries = "General"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
registries = "General"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a54cf54677c9f6eac0d97e9d757c233b56bc21b9"
registries = "General"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.13+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
registries = "General"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
registries = "General"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.2+0"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
registries = "General"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdeflate_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7ee221abc4f125ce18e73d48367cde0e16bc93b9"
registries = "General"
uuid = "46979653-d7f6-5232-b59e-dd310c4598de"
version = "1.26.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "1fcd6cf05d23adab0b5ffc496d020537215751dc"
registries = "General"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.59+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "e067c8bae65bb40866552296a2d5b4dc65b8928f"
registries = "General"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.80.702+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "52d3b9475133c3bc8c0a7f90f18bfc8cdb5a443a"
registries = "General"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.70.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.3+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

# ╔═╡ Cell order:
# ╟─c6ad0267-65d1-4372-a538-22acd9b5d02b
# ╟─2670973e-750c-47f0-9c6a-a58aadf87682
# ╟─05b485e7-115a-4dbb-aa73-0ca6ace2f5c0
# ╠═fd88a6c1-0abe-4a5a-9414-bb15730c9d18
# ╟─60c89d86-942e-4c97-bd7a-ad2f792b1155
# ╠═1d2293bd-a236-4b41-a9d7-9c27463b5062
# ╠═426c37ea-83dd-4ef4-8a19-f0e16536c034
# ╟─a9d14d8e-aefa-48d9-a3b3-ca13e284ab08
# ╠═1b239188-bba1-44d1-bc9c-10b60f762e0d
# ╠═f2c441b3-374b-439d-ba8c-4b28c0765d45
# ╟─cc736aec-c48a-411a-af83-4255309d77e9
# ╠═7ba07e6b-34e0-4206-9cf7-5006c43635a6
# ╠═d6e4fd70-8750-457f-a17f-9a044c25c482
# ╟─b7f7e8e2-7f26-4d27-8dc5-8204761a9965
# ╠═2c7a71f6-e618-40db-80b7-1ad5f09e59d5
# ╠═2a11aed4-6516-40dc-94fd-30aa94e76744
# ╟─984fcb2e-049b-4b35-a4b6-5c8c686921a3
# ╟─c5b4f340-c774-4f09-af4c-f326afce5de3
# ╠═2da8169e-a39d-4314-87c0-60b57f80ae96
# ╟─eb31da4c-7011-49ec-a9ac-194eaee10c7a
# ╠═892f19a6-af21-4353-9987-9de795ba7ad7
# ╠═64774c2c-0991-48b2-af4a-ced8fedd115b
# ╟─f7ceb7c4-6488-43eb-b211-9f7087762a56
# ╠═9f14f051-d1e9-4638-8008-840a886b50ab
# ╠═0cd8f84d-a6aa-46c7-9e13-535e8e0cea1f
# ╟─ce5aea87-2614-4e72-8e88-25a8c590b89f
# ╠═2d236d32-faf4-41ef-84de-1d40c02fb238
# ╠═b388f8be-3074-443c-97a9-72eb7098b5e3
# ╟─bb9d675b-00ce-4d20-b78b-a08835b47137
# ╠═a3baf049-5419-4773-8b7f-0ba07a9f1728
# ╠═efa15492-4811-4896-8947-b90191651952
# ╟─0ab125c6-568b-414c-acaa-abcfebb559e2
# ╠═bb2622f9-98a5-48ca-92cb-bbce84d797f1
# ╠═14909ad6-eb42-48ab-afa7-8b40b700fe3b
# ╟─5f438b76-7fc7-49ed-8a0a-7f2d5a760414
# ╠═9e13bde0-eccf-4d21-b603-fd186e87b7d0
# ╠═77a81c41-0034-4a18-bb70-d2a1f9379593
# ╟─5a9dfdf7-aaf8-421d-a451-6a6a5c35e1cc
# ╠═a13de611-145a-40da-9414-8f3f8f85ad98
# ╠═1e24cdfa-f4a9-4f9e-b95d-af83e724cd90
# ╟─a3f42e7e-af38-46b0-b1fb-bdbf77312e16
# ╠═095e7efe-f4c6-4a8a-a029-03ce97cd15bd
# ╠═e4a98c65-70de-42fc-b490-f019dd93d3fb
# ╟─2648b454-2762-4941-b13c-2303ffcd6521
# ╠═48fac714-ab96-4475-ad0b-0c61432bf849
# ╟─9e8eabc7-f1ce-4dbb-bfd5-e6d28793f9d6
# ╟─78740b86-45e9-45f5-b4c2-be6cb65f1368
# ╠═2c76c406-d15c-4271-baa4-ebebf5299429
# ╠═be30f3c4-0c20-4949-bcc8-f3812b39befc
# ╟─b86a31df-9de1-4e27-8ca7-ce7f3d2576fa
# ╠═aab4113d-d7ad-4521-a208-e192ad218cf0
# ╠═af22e5da-ff9b-4278-95d9-491f84a55add
# ╟─e13cbda5-ddd4-491e-bf2a-d6ed24d2bdbd
# ╠═539290af-787d-4deb-928c-50e4e9f28173
# ╠═86d608ec-aca2-4738-9c3e-656ccd562da6
# ╟─e86d6cf7-2274-403a-b1db-017973f33fb7
# ╟─00b9d37f-ce7b-4491-b1df-f0963f2598a8
# ╠═1f6a0d6c-832b-474d-8fbf-a9de6e821d80
# ╠═efa21003-ee7a-4834-b6f7-7625a7820f70
# ╟─ccb58b55-b9c8-4229-aff8-92541af123ec
# ╠═6fa07efa-fd25-4bca-bd18-8ba5c557d146
# ╠═95bf572e-23a0-4554-8ee3-b7932721e850
# ╟─7ab2fe80-812f-4a21-933a-e169a17f6c32
# ╠═6ddeb42d-37a5-48a3-ac82-47e4ec2ca541
# ╠═9537c91e-03a7-44d7-95ae-e60ea8d5496a
# ╟─4d4794fb-804d-4b4a-8c32-a32d88e43e30
# ╟─1fa8010d-9a5d-4297-9665-9fa8795ef5f7
# ╟─4ee79c7c-5c4e-457b-b713-103937e50355
# ╠═f2b222ec-0783-487e-9c52-835976a555b6
# ╠═f2c3768f-1c81-456d-ba81-6a91fc09e81b
# ╟─c749ce2d-17ae-45f4-b721-3f486b1cbc23
# ╟─3a3328d5-31de-4deb-bff7-d25c1fcbc4ef
# ╠═bf618161-ef96-4445-8fc8-25dc5f662242
# ╠═182861a3-2902-4e56-8fff-3bfc182268c7
# ╟─e0a2c745-41cf-4359-bb99-3117fbb507cc
# ╠═321fa5f5-be6a-4e5f-ba7f-4c117a240253
# ╠═2962b22a-5e4a-4dc1-9be6-8e57979ecba5
# ╟─7fa76f60-c623-4b25-b185-0cd5e805ef71
# ╠═ee46b865-dbd4-4939-b241-514941dd138d
# ╠═554321c7-cdcb-415c-8ac2-c63b42b6889a
# ╟─0af912a5-18ad-4974-a0b5-c1f66c8fa37a
# ╟─d0215082-52df-4eae-a409-7c3c7cfc69ed
# ╠═3a24f4aa-b074-4beb-b77d-6778e2fe580a
# ╠═98869c25-2644-47d1-b8cc-05699292f2a8
# ╟─d9e2b828-9c32-4638-bfcd-0c16b221aa43
# ╠═4092a893-818e-49f4-93d0-be7bd652dddc
# ╠═a0967b01-4eb8-4160-af68-a908734dc25e
# ╟─a0b64ba9-cbec-404e-bf39-aee03ae407ae
# ╠═82f11626-722b-46fa-9173-8b5d7a80a190
# ╠═392d3055-2ce0-4db7-8702-f79395f24b3b
# ╟─6aff587f-be6f-4fd8-96df-9c12f3769f32
# ╠═ef498536-8fb8-46e3-9d8d-f7eb3994a3f5
# ╠═84f948ad-f903-46a6-a520-2499864f348f
# ╟─63d629a1-379f-44bd-a02d-d252d0594357
# ╟─25598fcd-761d-40b8-95ac-8b68a00026da
# ╟─336f8626-7b0c-459d-bedb-281d010fcbd2
# ╟─eba2f06f-ebe9-492d-81d2-1cc4fccd5b0a
# ╠═de22c778-1529-4177-85e6-a0178f437a8c
# ╟─072535ef-9a0d-4c1e-aea7-6486d329f66f
# ╠═c262aa56-b9b6-4da4-8810-d6120f0724c6
# ╠═9a98c6bb-d0d5-426e-9137-9fc94fc944be
# ╟─00f69c40-f31c-49e6-a1cc-d84a4cb2b43d
# ╟─1e8ec5c5-3d41-44a7-8dfa-81d981750d9e
# ╟─dd22bc93-28bf-412e-a52b-d4332475daa2
# ╟─c6b35992-973b-49d4-bfa7-855e8fe10924
# ╠═522d88d6-5f60-401b-8786-0236c0859eda
# ╠═014b94b3-a6fa-4627-a24d-fe8c5e275582
# ╟─fbd0d75b-b3af-4eb5-a249-423317d656e5
# ╠═0eefe540-1669-4510-bef7-8fb7f8be68f1
# ╠═085c04f2-0762-4213-b874-14f3ce1f0b49
# ╟─a427567e-b6dd-44a9-a61f-a725960cc459
# ╠═b057812d-1a8e-479c-8192-1c1acdb97d23
# ╠═98f0753f-635e-46e3-8045-ff40ca2e827b
# ╟─cefdeee7-6a4b-4c60-a14c-c1edac53be98
# ╠═479b900b-b1bf-4d89-a53e-b1b20da7f5ad
# ╠═28237d12-6065-4fa3-bd40-6397dae7ecb6
# ╟─9c0862d0-8028-4ef4-a867-18c7932dad91
# ╟─e81f3b44-bc04-4768-81c2-c2733958f5a9
# ╠═88966b7e-b261-45a6-9f12-6b74f24c03e4
# ╠═1aaaebb6-839a-4598-82b9-dc17efdf7623
# ╟─802f5cad-a0c0-4426-94f2-426f89dea7e1
# ╟─2eed0bc8-2306-42bb-9803-33ae131021e6
# ╠═9922bdbe-cf9e-486a-8e8a-28a7ded12d04
# ╟─118aff27-dc6d-4686-b64a-ff840b731030
# ╠═ffc4bdc3-16e3-4b3e-b448-cfdc16428378
# ╟─71cca345-cffc-4aee-a671-08901a92babe
# ╠═1c3faf9f-adcf-4232-acf4-655d828623e3
# ╠═6708c2f5-a291-486e-8e66-93acfcedbcb7
# ╟─7cfb17fc-81a4-4522-a88b-2fd10001baaa
# ╠═ce1804da-df61-4e4c-aa3b-d950e83b7c13
# ╠═910920d2-085a-498f-8458-e6e96903ab92
# ╟─5590c909-6103-4427-b1a2-4cf35265dc07
# ╠═f7486c08-b74b-471f-bbc2-690439a487f8
# ╠═08c5e6fc-d8b5-4311-8f35-449a830daeec
# ╟─c143ab62-b3f0-4b6a-83aa-76e5dc3548ac
# ╟─0c05c385-67bd-4078-8fe2-34a37a13b312
# ╟─7d4f9a02-81a7-4bdc-a22d-3435192f9f15
# ╟─59b4d441-9a74-468f-ad8c-882516a09049
# ╠═bedc8ccd-e6f6-4dd1-a0b6-1889f4b5b658
# ╟─e0d2d6c4-d363-4bb2-9d12-42c4a52aba3b
# ╠═90fd3ea9-e115-41a5-b020-b8ade6cc6398
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

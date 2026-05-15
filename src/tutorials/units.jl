### A Pluto.jl notebook ###
# v0.20.25

#> [frontmatter]
#> image = "/assets/units.png"
#> order = 1
#> title = "Unit handling"
#> layout = "layout.jlhtml"
#> date = "2025-11-19"
#> description = "Work with units in astrophysical calculations."
#> tags = ["units", "plots", "radio astronomy", "data cubes"]

using Markdown
using InteractiveUtils

# ╔═╡ fd88a6c1-0abe-4a5a-9414-bb15730c9d18
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

    Pkg.add(["TOML", "StatsBase", "Distributions", "Random", "DynamicQuantities", "Unitful", "PhysicalConstants", "UnitfulEquivalences", "UnitfulAstro", "DimensionfulAngles", "PlutoUI", "HypertextLiteral"])

    # Statistical analysis
    using StatsBase: mean
    using Distributions: Normal
    using Random: Xoshiro

    # Units - DynamicQuantities
    using DynamicQuantities: DynamicQuantities as DQ
    using DynamicQuantities: SymbolicConstants as C_DQ

    # Units - Unitful
    using Unitful: Unitful as U
    using PhysicalConstants: CODATA2018 as C_U
    using UnitfulEquivalences: UnitfulEquivalences as UE
    import UnitfulAstro
    import DimensionfulAngles

    # Data Viz
    using CairoMakie: Colorbar, stephist, heatmap
    using Makie: Makie as M

    deps_ready = true
end

# ╔═╡ 90fd3ea9-e115-41a5-b020-b8ade6cc6398
begin
    deps_ready

    using TOML: TOML
    using PlutoUI: TableOfContents, details
    using HypertextLiteral: @htl
    using Test: @test
end

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
let
    fig, ax, p = heatmap(
        ras_DQ, decs_DQ, intcloud_DQ,
        axis = (;
            xreversed = true,
            xlabel = "RA",
            ylabel = "Dec",
        ),
        colormap = :cividis,
    )

    Colorbar(fig[1, 2], p; label = "Intensity")

    fig
end

# ╔═╡ f2c3768f-1c81-456d-ba81-6a91fc09e81b
let
    fig, ax, p = heatmap(
        ras_U, decs_U, intcloud_U,
        axis = (;
            xreversed = true,
            xlabel = "RA",
            ylabel = "Dec",
        ),
        colormap = :cividis,
    )

    Colorbar(fig[1, 2], p; label = "Intensity")

    fig
end

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

# ╔═╡ 3dbbae92-27cf-4573-8752-c2c400017812
function frontmatter(path)
    prefix = "#> "
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ dd1fc1c9-c55e-453a-bda7-a2036542cdcb
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
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

$(keywords())

!!! warning "Companion content"
	Content here.
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
# ╟─3dbbae92-27cf-4573-8752-c2c400017812
# ╟─dd1fc1c9-c55e-453a-bda7-a2036542cdcb
# ╠═90fd3ea9-e115-41a5-b020-b8ade6cc6398

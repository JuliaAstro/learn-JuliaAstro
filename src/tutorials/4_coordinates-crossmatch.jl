### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Astronomical Coordinates 4: Cross-matching Catalogs"
#> layout = "layout.jlhtml"
#> date = "2025-12-31"
#> description = "Retrieve the 2MASS catalog from VizieR over TAP, cross-match it against Gaia-selected members of NGC 188, and build color–magnitude diagrams from the matched photometry."
#> tags = ["coordinates", "gaia"]

using Markdown
using InteractiveUtils

# ╔═╡ 3fffcbae-34a3-4ee4-bdce-f7f9eda5fdae
begin
    # Can remove this block after Makie v0.25 is merged:
    # https://github.com/MakieOrg/Makie.jl/pull/5484
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "TOML"),
            Pkg.PackageSpec(; name = "Dates"),
            Pkg.PackageSpec(; name = "CSV"),
            Pkg.PackageSpec(; name = "DataFramesMeta"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; name = "VirtualObservatory"),
            Pkg.PackageSpec(; name = "NearestNeighbors"),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "Makie",
                rev = "ff/breaking-0.25",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "ComputePipeline",
                rev = "ff/breaking-0.25",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "CairoMakie",
                rev = "ff/breaking-0.25",
            ),
            Pkg.PackageSpec(;
                rev = "altaz",
                url = "https://github.com/activexray/SkyCoords.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaPhysics/DynamicQuantities.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/AstroTime.jl",
            ),
        ]
    )

    using CSV
    using DataFramesMeta
    using CairoMakie
    using DynamicQuantities
    using SkyCoords
    using NearestNeighbors: NearestNeighbors # activates SkyCoords.match
    using AstroTime
    using Dates: Dates
    using VirtualObservatory: execute, TAPService
    using LinearAlgebra: norm

    using DynamicQuantities.Units: °, yr
    using DynamicQuantities.Constants: pc
end

# ╔═╡ b37fc349-2c11-48b9-bea7-457e7a8fe650
begin
    using TOML: TOML
    using PlutoUI: TableOfContents
end

# ╔═╡ 209789ee-24f3-4c7e-9cf5-d60bcebf284f
md"""
## Summary
In the previous tutorials in this series, we introduced many of the key concepts underlying how to represent and transform astronomical coordinates in Julia, including how to work with both position and velocity data.

In this tutorial, we will explore how these tools can be used to cross-match two catalogs that contain overlapping sources that may have been observed at different times. You may find it helpful to keep the SkyCoords.jl documentation ([docs](https://juliaastro.org/SkyCoords/stable/)) open alongside this tutorial for reference or additional reading.

*Note: This is the 4th and final tutorial in a series of tutorials about coordinates in Julia. If you are new to this series, you may want to start from the beginning or an earlier tutorial.*
"""

# ╔═╡ 83a841a3-6881-447a-8c7a-8feb17f7de6d
md"""
### Packages 📦

We start by importing some general packages that we will need below:
"""

# ╔═╡ ec52515b-ba67-40e3-835e-47638fb80b2c
md"""
## Cross-matching and comparing catalogs

In this tutorial, we are going to return to the set of data that we downloaded from the *Gaia* archive back in the first tutorial of this series.

Let's recap what we did there: we defined a coordinate object for the center of the open cluster NGC 188, we queried the *Gaia* DR2 catalog to select stars that are close (on the sky) to the center of the cluster, and we used the parallax values from *Gaia* to select stars that are near NGC 188 in 3D position. Here, we will briefly reproduce those selections — with an additional cut on proper motion, using the velocity concepts from the previous tutorial — so that we can start with a catalog of sources that are likely members of NGC 188.
"""

# ╔═╡ 5a3c5a82-089f-484d-aa84-7315a09f399a
md"""
### Selecting candidate members of NGC 188

We start by re-loading the *Gaia* results table from the first tutorial, keeping only sources with a measured parallax larger than 0.25 mas:
"""

# ╔═╡ 391261c2-4140-4dd4-904e-228a73c95ba2
ngc188_table = @rsubset(
    dropmissing(
        CSV.read(joinpath(@__DIR__, "data", "gaia_results.csv"), DataFrame),
        [:parallax, :pmra, :pmdec],
    ),
    :parallax > 0.25, # mas
)

# ╔═╡ f9ed5b60-7c24-406a-a2d7-776e1bbe64de
md"""
This time around, our bundle for the cluster center holds its 3D position *and* its mean proper motion (taken from the literature), following the pattern we used for velocity data in the previous tutorial:
"""

# ╔═╡ 1ed20caf-6cb9-4247-8212-bca346151e06
ngc188_center = (;
    coord = ICRSCoords(12.11°, 85.26°),
    distance = 1.96e3pc |> us"Constants.pc",
    pm_ra_cosdec = -2.3087 * us"mas/yr",
    pm_dec = -0.9565 * us"mas/yr",
)

# ╔═╡ c7cb9c5f-edac-4157-ac4c-d8c19b0d1736
md"""
As in the first tutorial, we convert the parallaxes to distances and combine them with the sky positions to get the 3D position of every star, along with the 3D position of the cluster center:
"""

# ╔═╡ 12552458-a19b-4e2b-826d-becfa8caeadc
parallax_to_distance(ϖ_mas) = 1000 / ϖ_mas * pc

# ╔═╡ 76fc23e5-9c51-4367-9089-4c56625e0132
ngc188_coords_3d = let
    directions = ICRSCoords.(ngc188_table.ra * °, ngc188_table.dec * °)
    distances = parallax_to_distance.(ngc188_table.parallax)
    [cartesian(c).vec .* d for (c, d) in zip(directions, distances)]
end

# ╔═╡ 790f67a9-05b8-428a-93fe-6611fd824478
ngc188_center_3d = cartesian(ngc188_center.coord).vec .* ngc188_center.distance

# ╔═╡ 210cc51e-13eb-46ce-a129-88b86c0b191d
md"""
Cluster members should be close to the cluster center in 3D position *and* should be moving with the cluster, i.e., have a similar proper motion. We compute both the 3D separation from the cluster center and the total proper motion difference for every star:
"""

# ╔═╡ 8773d609-a612-47a2-bb65-202326d9e4a9
sep3d = norm.(ngc188_coords_3d .- Ref(ngc188_center_3d)) .|> us"Constants.pc"

# ╔═╡ abdbff93-d07b-4154-919b-ad2eef4bf4ab
pm_diff = hypot.(
    ngc188_table.pmra .- ustrip(us"mas/yr", ngc188_center.pm_ra_cosdec),
    ngc188_table.pmdec .- ustrip(us"mas/yr", ngc188_center.pm_dec),
) .* us"mas/yr"

# ╔═╡ a1606ee0-bc46-4f80-a6b8-70c1bb25a81d
md"""
Requiring both to be small selects our candidate cluster members:
"""

# ╔═╡ 05c6b210-98dd-4695-8597-bae1d2cbcaf6
ngc188_members_mask = (sep3d .< 50pc) .&& (pm_diff .< 1.5 * us"mas/yr")

# ╔═╡ 428e8762-3fac-4da6-93dd-ecaea5e7ea2a
ngc188_members = ngc188_table[ngc188_members_mask, :]

# ╔═╡ 7bbdb149-8a56-4110-927c-07f28d10ac6b
ngc188_members_coords = ICRSCoords.(ngc188_members.ra * °, ngc188_members.dec * °)

# ╔═╡ 6bc1f44d-80bb-4ae6-922b-6833fda53905
nrow(ngc188_members)

# ╔═╡ 8d013a10-6cf6-427a-aa71-4db0c9297aa7
md"""
From the selections above, the table `ngc188_members` and the coordinate vector `ngc188_members_coords` contain a couple hundred sources that, based on their 3D positions and proper motions, are consistent with being members of the open cluster NGC 188.
"""

# ╔═╡ 49d653f0-0d54-48a6-9ca2-c710caa1847f
md"""
### Querying the 2MASS catalog from VizieR

Let's assume that we now want to cross-match our catalog of candidate members of NGC 188 — here, based on *Gaia* data — to some other catalog. In this tutorial, we will demonstrate how to manually cross-match these *Gaia* sources with the 2MASS photometric catalog to retrieve infrared magnitudes for these stars, and then we will plot a color–magnitude diagram. To do this, we first need to query the 2MASS catalog to retrieve all sources in a region around the center of NGC 188, as we did for *Gaia*. Here, we will also take into account the fact that the *Gaia* data release 2 reference epoch is J2015.5, whereas the 2MASS coordinates are reported at their time of observation (in the late 1990's).

!!! note
    Some data archives, like the *Gaia* science archive, support running cross-matches at the database level and even support epoch propagation. If you need to perform a large cross-match, it will be much more efficient to use those services!

The main 2MASS photometric catalog lives in the [VizieR](https://vizier.cds.unistra.fr/) catalog collection under the identifier `II/246`. Like the *Gaia* archive and SIMBAD, VizieR exposes a TAP service — available in VirtualObservatory.jl as the built-in `:vizier` service — so we can express our cone search around the cluster center in ADQL just like before. We also request the `Date` column, which records the observation date of each source:
"""

# ╔═╡ 7609ca7f-dcea-4c89-a36c-22b6ee60014d
tmass_query = """
SELECT RAJ2000, DEJ2000, Jmag, Hmag, Kmag, Date
FROM "II/246/out"
WHERE 1 = CONTAINS(
    POINT('ICRS', RAJ2000, DEJ2000),
    CIRCLE('ICRS', $(rad2deg(ngc188_center.coord.ra)), $(rad2deg(ngc188_center.coord.dec)), 0.5)
)
"""

# ╔═╡ 427ec338-b963-45d1-abc1-0358ba46a59b
md"""
Running this query requires an internet connection, so we have included the results file next to this notebook. The cell below only submits the query (and saves the results to `data/2MASS_results.csv`) if that file is not already present:
"""

# ╔═╡ abb9d655-97a5-45e1-813f-13379ca565cf
tmass_table = let
    f = joinpath(@__DIR__, "data", "2MASS_results.csv")
    if !isfile(f)
        CSV.write(f, DataFrame(execute(TAPService(:vizier), tmass_query; unitful = false)))
    end
    CSV.read(f, DataFrame)
end

# ╔═╡ daf88313-b529-4013-a0ca-d814ccc25bda
md"""
As with the *Gaia* results table, we can now create a vector of coordinate objects to represent all of the sources returned from our query to the 2MASS catalog. Let's look at the column names in this table by displaying the first few rows:
"""

# ╔═╡ de46927b-4b84-4526-b8e1-283fd3329a9a
first(tmass_table, 3)

# ╔═╡ e5b374fc-1c32-4ab8-827f-c321bd299e15
md"""
From looking at the column names, the two relevant sky coordinate columns are `RAJ2000` for `ra` and `DEJ2000` for `dec`:
"""

# ╔═╡ 093efdd5-4abb-4759-935b-f41a17c4f17a
tmass_coords = ICRSCoords.(tmass_table.RAJ2000 * °, tmass_table.DEJ2000 * °)

# ╔═╡ 3d20cf17-0c2e-4d73-9cbf-088cff2d7511
length(tmass_coords)

# ╔═╡ 5902eb19-7c3d-4039-99bd-9cc2c3dd65dd
md"""
### Accounting for the difference in epochs

Note also that the table contains the `Date` column that specifies the epoch of the coordinates. Are all of these epochs the same?
"""

# ╔═╡ 65809c36-b524-4e34-a0dc-827fb570e8c1
tmass_date = (only ∘ unique)(tmass_table.Date) |> Dates.DateTime

# ╔═╡ 8de5f017-6c00-4c14-9e0b-b77dbfc70c76
md"""
It looks like all of the sources in our 2MASS table have the same epoch. CSV.jl already parsed this column into `Dates.Date` objects, so all that is left is to express this date as a Julian epoch year (the same time format we used for the *Gaia* reference epoch J2015.5 in the previous tutorial) — which is exactly what `jyear` from AstroTime.jl ([docs](https://juliaastro.org/AstroTime.jl/stable/)) computes for us, just like astropy's `Time(...).jyear`. Note that we interpret the calendar date directly as the timestamp of the epoch (via `TAIEpoch`) rather than converting it with `from_utc`, which would shift the clock reading by the accumulated leap seconds — at day precision the time scale label is immaterial, and this matches astropy's convention:
"""

# ╔═╡ 28269af8-2039-4f5c-9b14-c300fc644c39
tmass_epoch = tmass_date |> TAIEpoch

# ╔═╡ 4dd25f46-1179-46ee-ba97-1a831ae449c4
md"""
We now want to cross-match our *Gaia*-selected candidate members of NGC 188, `ngc188_members_coords`, with this table of photometry from 2MASS. However, as noted previously, the *Gaia* coordinates are given at a different epoch, J2015.5, which is nearly ~16 years after the 2MASS epoch of the data we downloaded (1999-10-19, or roughly J1999.8). We will therefore first use the `apply_space_motion` helper function that we wrote in the previous tutorial to transform the *Gaia* positions back to the 2MASS epoch before we do the cross-match, broadcasting it over the members with each star's own *Gaia* proper motion:
"""

# ╔═╡ eb6f3a12-7ae1-4d8e-9cec-361ee23d6cea
"""
    apply_space_motion(c::ICRSCoords, pm_ra_cosdec, pm_dec, Δt)

Linearly propagate the sky position `c` along its proper motion over the time interval `Δt`, returning the new position as an `ICRSCoords` object.
"""
function apply_space_motion(c::ICRSCoords, pm_ra_cosdec, pm_dec, Δt)
    Δra = pm_ra_cosdec * Δt / cos(c.dec)
    Δdec = pm_dec * Δt
    return ICRSCoords(c.ra + ustrip(u"rad", Δra), c.dec + ustrip(u"rad", Δdec))
end

# ╔═╡ 67bcf154-eb68-4c21-8163-c5f36a68b24b
ngc188_members_coords_1999 = apply_space_motion.(
    ngc188_members_coords,
    ngc188_members.pmra .* us"mas/yr",
    ngc188_members.pmdec .* us"mas/yr",
    # Switch from AstroTime.jl units to DynamicQuantities.jl units
    value(jyear(tmass_epoch) - 2015.5 * years) * yr,
)

# ╔═╡ dd4bb2ce-56a9-42b5-be76-37c76ebbfb30
md"""
The vector `ngc188_members_coords_1999` now contains the coordinates of our *Gaia*-selected members of NGC 188, as we think they would appear if observed on 1999-10-19.
"""

# ╔═╡ 3ada813d-9792-4b98-b492-9682a24d63ce
md"""
### Cross-matching the catalogs

In astropy, the cross-match itself is done with the `SkyCoord.match_to_catalog_sky` method. SkyCoords.jl provides the equivalent as `SkyCoords.match`, an efficient KD-tree-based matcher that lives in a [package extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)): it becomes available as soon as NearestNeighbors.jl ([docs](https://github.com/KristofferC/NearestNeighbors.jl)) is loaded alongside SkyCoords.jl, which is why we added it to our packages above. (The function is deliberately not exported, so it never clashes with `Base.match` — we call it as `SkyCoords.match`.)
"""

# ╔═╡ 6a63b5b9-d185-440b-bc61-af35ba800b0a
md"""
Note that argument order matters with this function, and it is the *reverse* of astropy's method: the reference catalog to search in comes first, followed by the coordinates to find matches for. Here we match our *Gaia* members against 2MASS, so `tmass_coords` is the reference and the returned indices point into it, along with the on-sky separation (in radians) of each match:
"""

# ╔═╡ 7e96dabc-0a41-447b-aed3-863ab2004af8
idx_gaia, sep2d_gaia = let
    i, s = SkyCoords.match(tmass_coords, ngc188_members_coords_1999)
    i, s * u"rad" |> us"arcsec"
end

# ╔═╡ f85ae247-59f1-4c4f-806d-c2e22706a2aa
md"""
Let's now look at the distribution of separations (in arcseconds) for all of the cross-matched sources:
"""

# ╔═╡ 258b2649-fe7e-49de-a80c-3972446e72e1
stephist(
    sep2d_gaia;
    bins = 10 .^ range(-2, 2, 64),
    axis = (;
        xlabel = "separation",
        ylabel = "counts",
        xscale = log10,
        yscale = log10,
    )
)

# ╔═╡ aae217a6-820e-466a-b510-dbb14fba8929
md"""
From this, it looks like the large majority of sources in our *Gaia* NGC 188 member list cross-match to another source within an arcsecond or so, so these all seem like they are correctly matched to a 2MASS source!
"""

# ╔═╡ 7d826fca-0a45-4505-bda9-fe6cb757283b
count(sep2d_gaia .< 2 * us"arcsec"), length(sep2d_gaia)

# ╔═╡ 0c884706-9718-4305-9d28-629d0dcaae30
md"""
### Color–magnitude diagrams

With our cross-match done, we can now make *Gaia*+2MASS color–magnitude diagrams of our candidate NGC 188 members using the information returned by the cross-match:
"""

# ╔═╡ 181f383a-fc85-4651-97f9-7546fe44a6b2
begin
    Jmag = tmass_table.Jmag[idx_gaia] # note that we use the index array returned above
    Gmag = ngc188_members.phot_g_mean_mag
    Bmag = ngc188_members.phot_bp_mean_mag
end

# ╔═╡ be7e1852-3531-4060-bb86-35bbe79c0cef
let
    fig = Figure(size = (900, 450))

    ax1 = Axis(
        fig[1, 1];
        xlabel = L"G - J",
        ylabel = L"G",
        yreversed = true, # backwards because magnitudes!
        limits = (0, 3, 10, 19),
    )
    scatter!(ax1, Gmag .- Jmag, Gmag; color = (:black, 0.5))

    ax2 = Axis(
        fig[1, 2];
        xlabel = L"G_\mathrm{BP} - G",
        ylabel = L"J",
        yreversed = true, # backwards because magnitudes!
        limits = (0.2, 1, 8, 17),
    )
    scatter!(ax2, Bmag .- Gmag, Jmag; color = (:black, 0.5))

    fig
end

# ╔═╡ 2dc2ec0b-fdfc-421a-bbc3-1096b0e073a9
md"""
Those both look like color–magnitude diagrams of a main sequence + red giant branch of an intermediate-age stellar cluster, so it looks like our selection and cross-matching has worked!

This concludes our series of tutorials about astronomical coordinates in Julia: we have represented and formatted sky positions with SkyCoords.jl and AstroAngles.jl, transformed them between frames, attached distances and velocities, queried the *Gaia*, SIMBAD, and VizieR archives over TAP with VirtualObservatory.jl, propagated positions between epochs, and cross-matched catalogs observed years apart. For more on what you can do with these coordinate objects, see the [SkyCoords.jl documentation](https://juliaastro.org/SkyCoords/stable/).
"""

# ╔═╡ a49724fc-50d6-4900-9cba-ca81e67e5935
md"""
# Notebook setup 🔧
"""

# ╔═╡ 31d9e987-b49b-4131-93b0-4a212c149f5d
TableOfContents(; title = "On this page", depth = 4)

# ╔═╡ 29856300-b8bd-43d4-aec8-5f43d6ad2a43
function frontmatter(path)
    prefix = "#> "
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ c2805f81-d01c-455e-870a-aa3ff07682b6
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ a128cd53-182c-4e66-a3a4-dceb9fc0cf18
md"""
# Astronomical Coordinates 4: Cross-matching Catalogs

This notebook is modified from <https://learn.astropy.org/tutorials/4_Coordinates-Crossmatch.html>

_Original authors: Adrian Price-Whelan_

!!! tip "Learning Goals"
    - Demonstrate how to retrieve a catalog from VizieR over TAP with VirtualObservatory.jl
    - Show how to perform positional cross-matches between catalogs of sky coordinates

$(keywords())
"""

# ╔═╡ Cell order:
# ╟─a128cd53-182c-4e66-a3a4-dceb9fc0cf18
# ╟─209789ee-24f3-4c7e-9cf5-d60bcebf284f
# ╟─83a841a3-6881-447a-8c7a-8feb17f7de6d
# ╠═3fffcbae-34a3-4ee4-bdce-f7f9eda5fdae
# ╟─ec52515b-ba67-40e3-835e-47638fb80b2c
# ╟─5a3c5a82-089f-484d-aa84-7315a09f399a
# ╠═391261c2-4140-4dd4-904e-228a73c95ba2
# ╟─f9ed5b60-7c24-406a-a2d7-776e1bbe64de
# ╠═1ed20caf-6cb9-4247-8212-bca346151e06
# ╟─c7cb9c5f-edac-4157-ac4c-d8c19b0d1736
# ╠═12552458-a19b-4e2b-826d-becfa8caeadc
# ╠═76fc23e5-9c51-4367-9089-4c56625e0132
# ╠═790f67a9-05b8-428a-93fe-6611fd824478
# ╟─210cc51e-13eb-46ce-a129-88b86c0b191d
# ╠═8773d609-a612-47a2-bb65-202326d9e4a9
# ╠═abdbff93-d07b-4154-919b-ad2eef4bf4ab
# ╟─a1606ee0-bc46-4f80-a6b8-70c1bb25a81d
# ╠═05c6b210-98dd-4695-8597-bae1d2cbcaf6
# ╠═428e8762-3fac-4da6-93dd-ecaea5e7ea2a
# ╠═7bbdb149-8a56-4110-927c-07f28d10ac6b
# ╠═6bc1f44d-80bb-4ae6-922b-6833fda53905
# ╟─8d013a10-6cf6-427a-aa71-4db0c9297aa7
# ╟─49d653f0-0d54-48a6-9ca2-c710caa1847f
# ╠═7609ca7f-dcea-4c89-a36c-22b6ee60014d
# ╟─427ec338-b963-45d1-abc1-0358ba46a59b
# ╠═abb9d655-97a5-45e1-813f-13379ca565cf
# ╟─daf88313-b529-4013-a0ca-d814ccc25bda
# ╠═de46927b-4b84-4526-b8e1-283fd3329a9a
# ╟─e5b374fc-1c32-4ab8-827f-c321bd299e15
# ╠═093efdd5-4abb-4759-935b-f41a17c4f17a
# ╠═3d20cf17-0c2e-4d73-9cbf-088cff2d7511
# ╟─5902eb19-7c3d-4039-99bd-9cc2c3dd65dd
# ╠═65809c36-b524-4e34-a0dc-827fb570e8c1
# ╟─8de5f017-6c00-4c14-9e0b-b77dbfc70c76
# ╠═28269af8-2039-4f5c-9b14-c300fc644c39
# ╟─4dd25f46-1179-46ee-ba97-1a831ae449c4
# ╠═eb6f3a12-7ae1-4d8e-9cec-361ee23d6cea
# ╠═67bcf154-eb68-4c21-8163-c5f36a68b24b
# ╟─dd4bb2ce-56a9-42b5-be76-37c76ebbfb30
# ╟─3ada813d-9792-4b98-b492-9682a24d63ce
# ╟─6a63b5b9-d185-440b-bc61-af35ba800b0a
# ╠═7e96dabc-0a41-447b-aed3-863ab2004af8
# ╟─f85ae247-59f1-4c4f-806d-c2e22706a2aa
# ╠═258b2649-fe7e-49de-a80c-3972446e72e1
# ╟─aae217a6-820e-466a-b510-dbb14fba8929
# ╠═7d826fca-0a45-4505-bda9-fe6cb757283b
# ╟─0c884706-9718-4305-9d28-629d0dcaae30
# ╠═181f383a-fc85-4651-97f9-7546fe44a6b2
# ╠═be7e1852-3531-4060-bb86-35bbe79c0cef
# ╟─2dc2ec0b-fdfc-421a-bbc3-1096b0e073a9
# ╟─a49724fc-50d6-4900-9cba-ca81e67e5935
# ╠═31d9e987-b49b-4131-93b0-4a212c149f5d
# ╟─29856300-b8bd-43d4-aec8-5f43d6ad2a43
# ╟─c2805f81-d01c-455e-870a-aa3ff07682b6
# ╠═b37fc349-2c11-48b9-bea7-457e7a8fe650

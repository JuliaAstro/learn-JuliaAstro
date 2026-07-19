### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> title = "Astronomical Coordinates 2: Transforming Coordinate Systems and Representations"
#> layout = "layout.jlhtml"
#> date = "2025-12-31"
#> description = "TODO"
#> tags = ["coordinates"]

using Markdown
using InteractiveUtils

# ╔═╡ 6f72fec9-eaf8-4831-8f59-49c4cc153f02
begin
    # Can remove this block after Makie v0.25 is merged:
    # https://github.com/MakieOrg/Makie.jl/pull/5484
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "Downloads"),
            Pkg.PackageSpec(; name = "TOML"),
            Pkg.PackageSpec(; name = "CSV"),
            Pkg.PackageSpec(; name = "DataFramesMeta"),
            Pkg.PackageSpec(; name = "PlutoUI"),
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
                rev = "makie-v0.25",
                url = "https://github.com/icweaver/GeoMakie.jl",
            ),
            Pkg.PackageSpec(;
                rev = "makie-v0.25",
                url = "https://github.com/icweaver/GeoInterface.jl",
            ),
            Pkg.PackageSpec(;
                rev = "makie-v0.25",
                url = "https://github.com/icweaver/GeoJSON.jl",
            ),
            Pkg.PackageSpec(;
                rev = "makie",
                url = "https://github.com/JuliaAstro/AstroImages.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/AstroAngles.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/FITSFiles.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/FITSWCS.jl",
            ),
            Pkg.PackageSpec(;
                rev = "altaz",
                url = "https://github.com/activexray/SkyCoords.jl",
            ),
            Pkg.PackageSpec(;
                rev = "makie-0.25",
                url = "https://github.com/icweaver/DimensionalData.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaPhysics/DynamicQuantities.jl",
            ),
            Pkg.PackageSpec(;
                rev = "sofa-v2-migration",
                url = "https://github.com/JuliaAstro/SOFA.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/AstroTime.jl",
            ),
        ]
    )

    using Downloads: download
    using CSV
    using DataFramesMeta
    using CairoMakie
    using GeoMakie
    using DynamicQuantities
    # using Unitful, UnitfulAstro
    # using Unitful: °
    # using UnitfulAstro: pc
    using SkyCoords
    using AstroAngles
    using SOFA
    using AstroTime

    using DynamicQuantities.Units: °
    using DynamicQuantities.Constants: pc
end

# ╔═╡ edf44d69-a40d-43ed-a4de-ceee627ffaab
using Dates

# ╔═╡ 56aaebbc-8166-4e52-a43e-93453137ad75
begin
    using TOML: TOML
    using PlutoUI: TableOfContents
end

# ╔═╡ d502528d-1b8b-46c0-9e46-5d3196cd3656
md"""
## Summary
In the previous tutorial in this series, we showed how astronomical coordinates in the ICRS or equatorial coordinate system can be represented in Julia using the SkyCoords.jl package ([docs](https://juliaastro.org/SkyCoords/)). There are many other coordinate systems that are commonly used in astronomical research. For example, the Galactic coordinate system is often used in radio astronomy and Galactic science, the "horizontal" or altitude-azimuth frame is often used for observatory-specific observation planning, and Ecliptic coordinates are often used for solar system science or space mission footprints. All of these coordinate frames (and others!) are supported by SkyCoords.jl. As we will see below, subtypes of the `AbstractSkyCoords` supertype in SkyCoords.jl are designed to make transforming between these systems a straightforward task.

In this tutorial, we will explore how the SkyCoords.jl package can be used to transform astronomical coordinates between different coordinate systems or frames.

!!! todo
    Create Part 1 of this tutorial series.
"""

# ╔═╡ f81fc70e-d1f3-410f-a6ce-cb6cd5cd3ca2
md"""
### Packages 📦

We start by importing some general packages that we will need below:
"""

# ╔═╡ 13b7b907-1005-4bce-9b0c-1787a8867f84
md"""
## Key Concepts: Component Formats, Representations, and Frames

Usage of the term "coordinates" is overloaded in astronomy and is often used interchangeably when referring to i) data formats (e.g., sexagesimal vs. decimal), ii) representations (e.g., Cartesian vs. spherical), and iii) frames (e.g., equatorial vs. galactic). In SkyCoords.jl, we have tried to formalize these three concepts and have made them a core part of the way we interact with objects in this package. Here we will give an overview of these different concepts as we build up to demonstrating how to transform between different astronomical reference frames or systems.
"""

# ╔═╡ cc989039-910a-4956-b725-cbe9592e0e22
md"""
### 1. Coordinate Component Formats

In our previous tutorial, we showed that it is possible to pass in coordinate component data to our `AbstractSkyCoords` subtype as strings or as unitful objects in a variety of formats and units with AstroAngles.jl ([docs](https://juliaastro.org/AstroAngles/stable/)). We also saw that the coordinate components of these objects can be re-formatted. For example, we can change the coordinate format by changing the component units, or converting the data to a string:
"""

# ╔═╡ aed06f42-4f0b-4b08-bb48-c120d129e54f
c = ICRSCoords(15.9932°, -10.52351344°)

# ╔═╡ 4d41090a-5774-451a-9d0b-846bcb3048e2
c.ra |> rad2ha # Hour angle

# ╔═╡ 1ff0b6f3-7748-4a79-b59d-645f8d68bb0b
format_angle(rad2hms(c.ra); delim = ["h", "m", "s"])

# ╔═╡ 104d045a-31fb-468a-80a3-d708c0b17085
format_angle(rad2dms(c.dec); delim = ["d", "m", "s"])

# ╔═╡ 11cc01df-b5a0-4deb-ae86-d1ae7d0f9e25
format_angle(rad2dms(c.dec); digits = 5)

# ╔═╡ 712fd752-8b85-43f9-a305-9f2e3ac23a2d
md"""
!!! todo
    Maybe add a convenience function extension for `format_angle(::AbstractSkyCoords; kwargs...)` in SkyCoords.jl at some point
"""

# ╔═╡ ceeca820-5bea-48dd-8fb3-884e0516650b
md"""
See the previous tutorial for more examples. See the documentation of [Printf.jl](https://docs.julialang.org/en/v1/stdlib/Printf/) and [Format.jl](https://github.com/JuliaString/Format.jl) for more fine-grain string formatting.
"""

# ╔═╡ 23e5a876-506f-4283-98b6-17b2af055850
md"""
### 2. Coordinate Representations

In the previous tutorial, we only worked with coordinate data in spherical representations (longitude/latitude), but SkyCoords.jl also supports other coordinate representations like Cartesian and cylindrical. To retrieve the coordinate data in a different representation, we can use the associated function. For example, to convert into Cartesian coordinates, we can use the `cartesian` function:
"""

# ╔═╡ eed3c3ed-8326-423f-8d4d-b72b73b56e9d
cartesian(c)

# ╔═╡ 5781172f-cb7b-4337-b2aa-40d7206cbafa
md"""
In the `AbstractSkyCoords` object `c` that we defined earlier, we only specified sky positions (i.e., no distance data), so the units of the Cartesian components that are returned above are dimensionless and are interpreted as being on the surface of the (dimensionless) unit sphere. If we instead pass in a distance to the constructor using the `distance` keyword argument, we get the 3D position with positional units. For example:
"""

# ╔═╡ cb25c15e-af7f-4987-9939-8b0d56dd1ca6
# c2 = ICRSCoords(ra = 15.9932°, dec = -10.52351344°, distance = 127.4pc)

# ╔═╡ ff03249b-5d0f-44fe-8205-0e0a24170937
c2 = cartesian(c).vec * 127.4 .* pc # Current workaround

# ╔═╡ 2457add6-4c77-42e2-a2bb-ba9a3a3a3000
c2 .|> us"Constants.pc"

# ╔═╡ 3a13ae4f-714e-408d-993f-c3a0859496dc
md"""
Or, we could represent this data with cylindrical components:
"""

# ╔═╡ e6868fd1-9bcf-4f8a-81ca-c5df66bc5e1d
# Implement cylindrical(c2)

# ╔═╡ 722812b9-4984-4f26-b193-786a42242cf9
md"""
To summarize, using `cartesian`/`spherical`/`cylindrical` is a convenient way to retrieve your coordinate data in a different representation.
"""

# ╔═╡ dfc16a40-f501-4a76-b014-0721e0f1b7e4
md"""
!!! todo
    - Implement cylindrical coords
    - Document `cartesian` and `spherical` [#84](https://github.com/JuliaAstro/SkyCoords.jl/issues/84)
    - Add `distance` kwarg to `AbstractSkyCoords` constructors
"""

# ╔═╡ 92ae5f5b-1a81-44d9-b4b6-14f53e2a5817
md"""
### 3. Transforming Between Coordinate Frames

The third key concept to keep in mind when thinking about astronomical coordinate data is the reference frame or coordinate system that the data are in. In the previous tutorial, and so far here, we have worked with the International Celestial Reference System (ICRS; [some important definitions and context about the ICRS is given here](https://arxiv.org/abs/astro-ph/0602086)). The ICRS is the fundamental coordinate system used in most modern astronomical contexts and is generally what people mean when they refer to “equatorial” or “J2000” or “RA/Dec” coordinates (but there are some important caveats if you are working with older data). As noted above, however, there are many other coordinate systems used in different astronomical, solar, or solar system contexts.
"""

# ╔═╡ ba2d2e02-6541-401c-8746-5f7ec7bb230b
md"""
Some other common coordinate systems are defined as a rotation away from the ICRS that is defined to make science applications easier to interpret. One example here is the Galactic coordinate system, which is rotated with respect to the ICRS to approximately align the Galactic plane with latitude=0. As an example of the SkyCoords.jl transformation machinery, we will load in a subset of a catalog of positions and distances to a set of open clusters in the Milky Way from Cantat-Gaudin et al. 2018 ([Table 1 in this catalog](http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/A%2bA/618/A93/table1)). We have pre-selected the 474 clusters within 2 kpc of the sun and provide the catalog as a data file next to this notebook. This catalog provides sky position (columns `RAJ2000` and `DEJ2000` in the original catalog) and distance estimates (column `dmode` in the original catalog), which we have renamed in the table we provide to column names `"ra"`, `"dec"`, and `"distance"`. We will start by loading the catalog as a `DataFrames.DataFrame` ([docs](https://dataframes.juliadata.org/stable/)):
"""

# ╔═╡ cc5f05a7-c80f-4646-865f-e0bb8f3591c1
filename = let
    url = "https://raw.githubusercontent.com/astropy-learn/tutorial--astropy-coordinates/02ce4ab85a818dc8a6244e1a4285eea92e5649c7/Cantat-Gaudin-open-clusters.ecsv"
    f = joinpath(@__DIR__, "data", "Cantat-Gaudin-open-clusters.ecsv")
    isfile(f) ? f : download(url, f)
end;

# ╔═╡ 7653593c-531b-4749-a2dc-91078faa11f3
tbl = CSV.read(filename, DataFrame; comment = "#")

# ╔═╡ f031ba93-4ab8-4b26-acb9-3c9c162a3a8b
md"""
We next create a vector of `ICRSCoords` to represent the positions of all of the open clusters in this catalog:
"""

# ╔═╡ 75ca894d-df36-4b8e-a41c-c83eeaa6de2b
open_cluster_c = ICRSCoords.(tbl.ra * °, tbl.dec * °)

# ╔═╡ 82dc476a-d1b8-45c2-a104-ca650f5f1c75
length(open_cluster_c)

# ╔═╡ 703baf6d-19c2-4fee-89a3-f1e5e5e43967
md"""
Here are the first few coordinate entries:
"""

# ╔═╡ e000650c-fdc2-4dfc-8d1e-cc92fde93a4f
open_cluster_c[begin:4]

# ╔═╡ 2e96000f-d39a-4052-92ec-a76f5de67ddd
md"""
Let’s now visualize the sky positions of all of these clusters, colored by their distances. We plot these in an all-sky spherical projection (e.g., aitoff) using GeoMakie.jl ([docs](https://geo.makie.org/stable/)), with longitude increasing to the left as is typically done for plotting astronomical objects on the sky:
"""

# ╔═╡ 149ac51b-a571-4ab6-b213-decce18c06d1
# !!!! NO!!
begin
    Base.rad2deg(c::ICRSCoords) = (rad2deg(c.ra), rad2deg(c.dec))
    Base.rad2deg(c::GalCoords) = (rad2deg(c.l), rad2deg(c.b))
end

# ╔═╡ 098d72bc-45c9-4e65-8a09-c496ffca64b3
wrap180(x) = mod(x + 180, 360) - 180

# ╔═╡ 397fb01b-5120-48ee-b295-6105a574a36a
with_theme(fontsize = 12) do
    fig = Figure()

    ax = GeoAxis(
        fig[1, 1];
        dest = "+proj=aitoff",
        xlabel = "RA [deg]",
        ylabel = "Dec [deg]",
        xreversed = true,
        xticklabelsvisible = false,
        # yticklabelsvisible = false,
        # yaxisposition = :left,
    )

    # Manually label x-axis for now
    for ra in filter(!=(180), 0:30:330)
        text!(
            ax, wrap180(ra), 0;
            text = "$(ra)°",
            align = (:center, :top),
            # offset = (0, -4),
            # fontsize = 12,
        )
    end

    plt = scatter!(
        ax, rad2deg.(open_cluster_c);
        color = tbl.distance,
    )

    Colorbar(
        fig[1, 2], plt;
        label = "distance [pc]",
        ticks = 0:250:2000,
    )

    # TODO: Workaround until ax.yaxisposition = :left is working again
    colgap!(fig.layout, -20)

    # Remove vertical gaps
    rowsize!(fig.layout, 1, Aspect(1, 0.5))
    resize_to_layout!(fig)

    fig
end

# ╔═╡ e65c1300-bb4a-475e-bdb7-c6fe2f5bed0f
md"""
!!! warning
    Missing axis tick labels issue tracked here:
    - axis tick labels: <https://github.com/MakieOrg/GeoMakie.jl/issues/134>
    - axis labels: <https://github.com/MakieOrg/GeoMakie.jl/issues/221>
"""

# ╔═╡ f7657425-a75f-4628-898f-d6e2320df3e1
md"""
The majority of these open clusters are relatively close to the Galactic midplane, which is why they form a fairly narrow “band” on the sky in ICRS coordinates. If we transform these positions to Galactic coordinates, we would therefore expect the points to appear around the latitude ``b = 0`` line.
"""

# ╔═╡ af45cb56-eb28-4224-9c43-97675b370e3e
md"""
To transform our coordinates from ICRS to Galactic (or any other coordinate system), we can use the corresponding constructor, e.g., [`GalCoords`](https://juliaastro.org/SkyCoords/stable/api/#SkyCoords.GalCoords):
"""

# ╔═╡ 70f3d3b6-be88-4e25-9dc9-f21e574f4222
open_cluster_gal = GalCoords.(open_cluster_c)

# ╔═╡ 43986706-9ab5-403f-8f8b-e51c112f3689
md"""
which is just syntactic sugar for `convert.(GalCoords, open_cluster_c)`.
"""

# ╔═╡ e361cc88-7c1a-46f2-baca-ac1d64fb1a35
md"""
Comparing this to the original `ICRCoords`, note that the names of the longitude and latitude components have changed from `ra` to `l` and from `dec` to `b`, per convention. We can therefore access the new Galactic longitude/latitude data using these new attribute names:
"""

# ╔═╡ 468384b2-725e-4a93-9a8c-aacaa5a1bbc1
map(open_cluster_gal[begin:4]) do c
    format_angle(rad2dms(c.l); delim = ["°", "′", "″"])
end

# ╔═╡ ac07f24c-ba3b-4fa0-a26d-4cd743a7a079
map(open_cluster_gal[begin:4]) do c
    format_angle(rad2dms(c.b); delim = ["°", "′", "″"])
end

# ╔═╡ ec546913-3da1-4604-ae59-ad8251656798
md"""
!!! tip
    We can use any delimeters we prefer. In the example above, we use `\degree`, `\prime`, and `\pprime` which can be entered by tab completing them in the notebook/REPL.
"""

# ╔═╡ 7ab09154-c5ee-4310-8fad-5213965afed9
md"""
With this new AbstractSkyCoords object (in the Galactic frame), let's re-make a sky plot to visualize the sky positions of the open clusters in Galactic coordinates:
"""

# ╔═╡ 7abac62e-39ab-467e-95d2-e4f7d8eb4371
with_theme(fontsize = 12) do
    fig = Figure()

    ax = GeoAxis(
        fig[1, 1];
        dest = "+proj=aitoff",
        xlabel = "RA [deg]",
        ylabel = "Dec [deg]",
        xreversed = true,
        xticklabelsvisible = false,
        # yticklabelsvisible = false,
        yaxisposition = :left,
        limits = (nothing, nothing, -90, 90),
    )

    # Manually label x-axis for now
    for ra in filter(!=(180), 0:30:330)
        text!(
            ax, wrap180(ra), 0;
            text = "$(ra)°",
            align = (:center, :top),
            offset = (0, -45),
            # fontsize = 12,
        )
    end

    plt = scatter!(
        ax, rad2deg.(open_cluster_gal);
        color = tbl.distance,
    )

    Colorbar(
        fig[1, 2], plt;
        label = "distance [pc]",
        ticks = 0:250:2000,
    )

    # ylims!(-90u"°", 90u"°") # TODO: currently doesn't work

    # Workaround until ax.yaxisposition = :left is working again
    # colgap!(fig.layout, -20)
    # ax.yaxisposition = :left

    # Remove vertical gaps
    rowsize!(fig.layout, 1, Aspect(1, 0.5))
    resize_to_layout!(fig)

    fig
end

# ╔═╡ afe3b38d-7381-403b-8dd0-de56e4c9309c
md"""
As we hoped and expected, in the Galactic coordinate frame, the open clusters predominantly appear at low galactic latitudes!
"""

# ╔═╡ ef9592cc-673a-488f-9e8f-fe9314f6d8b6
md"""
## Transforming to More Complex Coordinate Frames: Computing the Altitude of a Target at an Observatory
"""

# ╔═╡ 605fd406-3edc-476a-ac66-8a3d5e5b6661
md"""
!!! warning "Under construction"
    This section is not available until <https://github.com/JuliaAstro/SkyCoords.jl/issues/35> is closed.
"""

# ╔═╡ e7b5928b-35b3-438d-8b07-7d79899c6df6
md"""
To determine whether a target is observable from a given observatory on Earth or to find out what targets are observable from a city or place on Earth at some time, we sometimes need to convert a coordinate or set of coordinates to a frame that is local to an on-earth observer. The most common choice for such a frame is “horizontal” or “Altitude-Azimuth” coordinates. In this frame, the sky coordinates of a source can be specified as an altitude from the horizon and an azimuth angle at a specified time. This coordinate frame is supported in SkyCoords.jl through the [`AltAzCoords`](https://juliaastro.org/SkyCoords/stable/api/#SkyCoords.AltAzCoords) coordinate frame.
"""

# ╔═╡ 9faa3737-b931-43c9-85cb-c70e63e507ce
md"""
The `AltAzCoords` frame is different from the previously-demonstrated `GalCoords` frame in that it requires additional metadata to define the frame instance. Since the Galactic frame is close to being a 3D rotation away from the ICRS frame, and that rotation matrix is fixed, we could transform to Galactic by converting with no arguments (see the example above where we used `GalCoords(<ICRSCoords>)`. In order to specify an instance of the `AltAzCoords` frame, we have to (at minimum) pass in (1) a location on Earth, and (2) the time we are requesting the frame at.
"""

# ╔═╡ e3cae13e-1c09-4ae0-9b77-f7db3b942075
md"""
We will use Kitt Peak National Observatory (in AZ, USA) Kitt as our site by storing its latitude and longitude in an [`Observer`](https://juliaastro.org/SkyCoords/stable/api/#SkyCoords.Observer) object:
"""

# ╔═╡ f1a88ead-193f-4ec0-a431-544ea9fc4b0c
obs_location = Observer(deg2rad(31.96333333), deg2rad(-111.6))

# ╔═╡ fa0f8370-dff2-4ecc-9cdc-d0409cbbcfe9
md"""
As an example, we will now compute the altitude of a few of the open clusters from our catalog above over the course of a night. We have an object (`obs_location`) to represent our location on Earth, so now we need to create a set of times to compute the `AltAzCoords` frame for. `AltAzCoords` expects time information to be passed in UTC JD. Let's pretend we had an observing run coming up on Dec 18, 2020, and we would like to compute the altitude/azimuth coordinates for our open clusters over that whole night:

!!! todo
    Support units, AstroTime.jl, etc.
"""

# ╔═╡ 922f4393-d3ca-4d2a-aec4-bcba885e4f6b
# 1AM UTC = 6PM local time (AZ mountain time), roughly the start of a night
obs_date = from_utc("2020-12-18T01:00:00")

# ╔═╡ 38754f06-c9a0-41ee-9a03-0b1e4882b9e4
# Compute the alt/az over a 14 hour period, starting at 6PM local time,
# with 256 equally spaced time points:
times = obs_date .+ range(0, 14, 256) .* hours

# ╔═╡ c26b86b5-4662-4285-bfb1-f381f0aa9957
md"""
Now we use our location, `obs_location`, and this grid of times, `times`, to transform our ICRS positions to their corresponding alt/az. Let's do this for the first open cluster in our catalog:
"""

# ╔═╡ fa106a43-e57e-4d11-9f94-a89eaa53ed71
# Should probably make an AstroTime extension in SkyCoords.jl or SOFA.jl for this
jds = (julian ∘ to_utc).(AstroDates.DateTime, times)

# ╔═╡ 721f0a73-ac9a-4e47-bab0-7eefd81e3ad1
altaz(c, jd; obs = obs_location) = AltAzCoords(c, obs, jd)

# ╔═╡ f7cba3b2-95bd-46f1-ad03-530a06c9239f
md"""
!!! note
    `AltAzCoords` accepts even more parameters about the atmosphere, which can be used to correct for atmospheric refraction. But here we leave those additional parameters set to their defaults, which ignores refraction.
"""

# ╔═╡ 9ef4da1b-215c-4678-adfb-59df132afec5
md"""
Let's now plot the altitude of this open cluster over the course of the night:
"""

# ╔═╡ a7b9e31f-f529-48cf-8ae6-cc0891040d8f
alts_cluster1 = [altaz(open_cluster_c[1], jd).alt for jd in jds] * u"rad"

# ╔═╡ a11fb0a5-ecbd-4818-92d7-cec37318e6f9
lines(
    to_utc.(Dates.DateTime, times), alts_cluster1;
    axis = (;
        xlabel = "Date/Time [UTC]",
        ylabel = "Altitude",
        dim2_conversion = Makie.DQConversion(us"°"),
    )
)

# ╔═╡ fb0c9eba-20cb-475b-a044-d199095ce325
n_clusters, n_times = length(open_cluster_c), length(times)

# ╔═╡ e857a4d0-de14-44d1-aeef-14df00a542a8
md"""
Here we can see that this open cluster reaches a high altitude above the horizon from Kitt Peak, and so it looks like it would be observable from this site. The above curve only shows the altitude trajectory for the first open cluster in our catalog, but we would like to compute the equivalent for all of the open clusters in the catalog. We have $(n_clusters) open clusters and we want to evaluate the alt/az coordinates of these clusters at $(n_times) different times:
"""

# ╔═╡ a5dc2919-6a35-4495-b341-510766d5946c
"""
We therefore want to produce a two-dimensional coordinate object that is indexed along one axis by the open cluster index, and along other axis by the time index. We _could_ do this by broadcasting over two dimensions:

```julia
altaz_grid = altaz.(open_cluster_c, permutedims(jds))
```

but this would be pretty slow because we would be computing things like precession-nutation, Earth rotation, polar motion, and refraction constants for each of our $(n_clusters) × $(n_times) = $(n_clusters * n_times) total points. There are only **$(n_times)** distinct time points to consider, so there is a more efficient route that we can go.
""" |> Markdown.parse

# ╔═╡ f67881ba-1f84-4380-9c87-d365e0a82ef7
md"""
The strategy is to use SOFA.jl's two-stage idiom: i) build an astrometry context once per time with `apco13`, then ii) run cheap per-star transforms (`atciqz` --> `atioq`) against it. We implement this below:
"""

# ╔═╡ 8080d735-b994-474d-9d8f-8d8e985d3c03
# quick: one apco13 per TIME, then cheap per-star transforms
function altaz_fast(cl, obs, jds)
    out = Matrix{AltAzCoords{Float64}}(undef, length(cl), length(jds))
    for (j, jd) in pairs(jds)
        r = SOFA.apco13(
            float(jd), 0.0, 0.0,
            float(obs.longitude), float(obs.latitude), float(obs.altitude),
            0.0, 0.0, 0.0, 0.0, 0.0, 0.55,
        )
        astrom = r[1]
        for (i, c) in pairs(cl)
            icrs = convert(ICRSCoords, c)
            ri, di = SOFA.atciqz(icrs.ra, icrs.dec, astrom)
            o = SOFA.atioq(ri, di, astrom)
            out[i, j] = AltAzCoords(π / 2 - o[2], o[1])
        end
    end
    return out
end

# ╔═╡ a3aabdfa-6a90-4b2b-9b43-179166c1fac1
md"""
!!! todo
    Upstream a more general version of this to the SOFA.jl extension in SkyCoords.jl
"""

# ╔═╡ 0d43b6fc-530d-4f44-8ba5-4d0de7cec52d
altaz_grid = altaz_fast(open_cluster_c, obs_location, jds)

# ╔═╡ 2857dc98-1780-4689-a172-9d9b08f906e7
md"""
Let’s now over-plot the trajectories for the first 10 open clusters in the catalog:
"""

# ╔═╡ ce47bdd9-dc4d-47ee-9430-c11cbcab3fff
# This is an alternative to `alts_cluster1`, which handled the
# unit conversion to degrees within Makie.jl
alts_grid = getproperty.(altaz_grid, :alt) * u"rad" |> us"°"

# ╔═╡ f537cb9e-a465-48a1-80a4-a8b79f29e150
series(
    to_utc.(Dates.DateTime, times),
    alts_grid[begin:10, :];
    color = :Spectral,
    axis = (;
        xlabel = "Date/Time [UTC]",
        ylabel = "Altitude",
    ),
)

# ╔═╡ 5bc81327-4ed8-4c78-9856-7402d370a044
md"""
From this, we can see that only two of the clusters in this batch seem to be easily observable.
"""

# ╔═╡ 4d8cf50e-0735-4114-ac72-b50ea2565647
md"""
## Potential Caveats
"""

# ╔═╡ eef0b9d5-32e4-4d2e-9116-7b44ce0e1567
md"""
Transformations between some reference frames require knowing more information about a source. For example, the transformation from ICRS to Galactic coordinates (as demonstrated above) amounts to a 3D rotation, but no change of origin. This transformation is therefore supported for any spherical position (with or without distance information). However, some transformations include a change of origin, and therefore require that the source data (i.e., the SkyCoord object) has defined distance information. For example, for a SkyCoord with only sky position, we can transform from the ICRS to the FK5 coordinate system:
"""

# ╔═╡ 765aee74-6bb5-40c2-859d-c55ef6860bfa


# ╔═╡ 0a9854b8-d5e1-4b5f-8d22-328a69e9da30
md"""
However, we would NOT be able to transform this position to the Galactocentric frame (docs), because this transformation involves a shift of origin from the solar system barycenter to the Galactic center:
"""

# ╔═╡ 1a853d21-a144-4486-9956-6d9d5f71620f
md"""
# Exercises
"""

# ╔═╡ 82f44845-5bd4-48ed-b04e-a0cf3c6ca240
md"""
Przybylski’s star or HD101065 is in the southern constellation of Centaurus with a right ascension of 174.4040348 degrees and a declination of -46.70953633 degrees. Create a SkyCoord object of its sky position.
"""

# ╔═╡ a2970366-b716-4574-b93f-49d18372e108
md"""
If the distance to Przybylski’s star is 108.4 pc, retrieve a 3D Cartesian representation. (Hint: we did this earlier in the tutorial and it may help to create a new 3D SkyCoord object.)
"""

# ╔═╡ d3f76adc-ba32-4862-a9c1-35d1a20e3474
md"""
Imagine it is May 2018, and you would like to take an observation of HD 101065 from Greenwich Royal Observatory. Use SkyCoords to figure out if you can observe the star that month. You can use any time and date of that month for your timeframe.
"""

# ╔═╡ 03b572a4-e51b-41e7-bc38-247644e41ebd
md"""
# Notebook setup 🔧
"""

# ╔═╡ 4540e3ef-9db5-4fa7-bf6b-d5f8dd4244c4
TableOfContents(; title = "On this page", depth = 4)

# ╔═╡ f404cd49-186b-493e-bee8-fb08c88f3f88
function frontmatter(path)
    prefix = "#> "
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ e24cb37e-acfd-441c-9839-40649611c1c7
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ dd98b24e-610e-11ef-1180-ef02be7d7cac
md"""
# Astronomical Coordinates 2: Transforming Coordinate Systems and Representations

This notebook is modified from <https://learn.astropy.org/tutorials/2_Coordinates-Transforms.html>

_Original authors: Adrian Price-Whelan, Saima Siddiqui, Zihao Chen, Luthien Liu_

!!! tip "Learning Goals"
    - Introduce key concepts in SkyCoords.jl: coordinate component formats, representations, and frames
    - Demonstrate how to work with coordinate representations, for example, to change from Cartesian to Cylindrical coordinates
    - Introduce coordinate frame transformations and demonstrate transforming from ICRS coordinates to Galactic and Altitude-Azimuth coordinates

$(keywords())
"""

# ╔═╡ Cell order:
# ╟─dd98b24e-610e-11ef-1180-ef02be7d7cac
# ╟─d502528d-1b8b-46c0-9e46-5d3196cd3656
# ╟─f81fc70e-d1f3-410f-a6ce-cb6cd5cd3ca2
# ╠═6f72fec9-eaf8-4831-8f59-49c4cc153f02
# ╟─13b7b907-1005-4bce-9b0c-1787a8867f84
# ╟─cc989039-910a-4956-b725-cbe9592e0e22
# ╠═aed06f42-4f0b-4b08-bb48-c120d129e54f
# ╠═4d41090a-5774-451a-9d0b-846bcb3048e2
# ╠═1ff0b6f3-7748-4a79-b59d-645f8d68bb0b
# ╠═104d045a-31fb-468a-80a3-d708c0b17085
# ╠═11cc01df-b5a0-4deb-ae86-d1ae7d0f9e25
# ╟─712fd752-8b85-43f9-a305-9f2e3ac23a2d
# ╟─ceeca820-5bea-48dd-8fb3-884e0516650b
# ╟─23e5a876-506f-4283-98b6-17b2af055850
# ╠═eed3c3ed-8326-423f-8d4d-b72b73b56e9d
# ╟─5781172f-cb7b-4337-b2aa-40d7206cbafa
# ╠═cb25c15e-af7f-4987-9939-8b0d56dd1ca6
# ╠═ff03249b-5d0f-44fe-8205-0e0a24170937
# ╠═2457add6-4c77-42e2-a2bb-ba9a3a3a3000
# ╟─3a13ae4f-714e-408d-993f-c3a0859496dc
# ╠═e6868fd1-9bcf-4f8a-81ca-c5df66bc5e1d
# ╟─722812b9-4984-4f26-b193-786a42242cf9
# ╟─dfc16a40-f501-4a76-b014-0721e0f1b7e4
# ╟─92ae5f5b-1a81-44d9-b4b6-14f53e2a5817
# ╟─ba2d2e02-6541-401c-8746-5f7ec7bb230b
# ╠═cc5f05a7-c80f-4646-865f-e0bb8f3591c1
# ╠═7653593c-531b-4749-a2dc-91078faa11f3
# ╟─f031ba93-4ab8-4b26-acb9-3c9c162a3a8b
# ╠═75ca894d-df36-4b8e-a41c-c83eeaa6de2b
# ╠═82dc476a-d1b8-45c2-a104-ca650f5f1c75
# ╟─703baf6d-19c2-4fee-89a3-f1e5e5e43967
# ╠═e000650c-fdc2-4dfc-8d1e-cc92fde93a4f
# ╟─2e96000f-d39a-4052-92ec-a76f5de67ddd
# ╠═149ac51b-a571-4ab6-b213-decce18c06d1
# ╠═397fb01b-5120-48ee-b295-6105a574a36a
# ╟─098d72bc-45c9-4e65-8a09-c496ffca64b3
# ╟─e65c1300-bb4a-475e-bdb7-c6fe2f5bed0f
# ╟─f7657425-a75f-4628-898f-d6e2320df3e1
# ╟─af45cb56-eb28-4224-9c43-97675b370e3e
# ╠═70f3d3b6-be88-4e25-9dc9-f21e574f4222
# ╟─43986706-9ab5-403f-8f8b-e51c112f3689
# ╟─e361cc88-7c1a-46f2-baca-ac1d64fb1a35
# ╠═468384b2-725e-4a93-9a8c-aacaa5a1bbc1
# ╠═ac07f24c-ba3b-4fa0-a26d-4cd743a7a079
# ╟─ec546913-3da1-4604-ae59-ad8251656798
# ╟─7ab09154-c5ee-4310-8fad-5213965afed9
# ╠═7abac62e-39ab-467e-95d2-e4f7d8eb4371
# ╟─afe3b38d-7381-403b-8dd0-de56e4c9309c
# ╟─ef9592cc-673a-488f-9e8f-fe9314f6d8b6
# ╟─605fd406-3edc-476a-ac66-8a3d5e5b6661
# ╟─e7b5928b-35b3-438d-8b07-7d79899c6df6
# ╟─9faa3737-b931-43c9-85cb-c70e63e507ce
# ╟─e3cae13e-1c09-4ae0-9b77-f7db3b942075
# ╠═f1a88ead-193f-4ec0-a431-544ea9fc4b0c
# ╟─fa0f8370-dff2-4ecc-9cdc-d0409cbbcfe9
# ╠═922f4393-d3ca-4d2a-aec4-bcba885e4f6b
# ╠═38754f06-c9a0-41ee-9a03-0b1e4882b9e4
# ╟─c26b86b5-4662-4285-bfb1-f381f0aa9957
# ╠═fa106a43-e57e-4d11-9f94-a89eaa53ed71
# ╠═721f0a73-ac9a-4e47-bab0-7eefd81e3ad1
# ╟─f7cba3b2-95bd-46f1-ad03-530a06c9239f
# ╟─9ef4da1b-215c-4678-adfb-59df132afec5
# ╠═a7b9e31f-f529-48cf-8ae6-cc0891040d8f
# ╠═edf44d69-a40d-43ed-a4de-ceee627ffaab
# ╠═a11fb0a5-ecbd-4818-92d7-cec37318e6f9
# ╟─e857a4d0-de14-44d1-aeef-14df00a542a8
# ╠═fb0c9eba-20cb-475b-a044-d199095ce325
# ╟─a5dc2919-6a35-4495-b341-510766d5946c
# ╟─f67881ba-1f84-4380-9c87-d365e0a82ef7
# ╠═8080d735-b994-474d-9d8f-8d8e985d3c03
# ╟─a3aabdfa-6a90-4b2b-9b43-179166c1fac1
# ╠═0d43b6fc-530d-4f44-8ba5-4d0de7cec52d
# ╟─2857dc98-1780-4689-a172-9d9b08f906e7
# ╠═ce47bdd9-dc4d-47ee-9430-c11cbcab3fff
# ╠═f537cb9e-a465-48a1-80a4-a8b79f29e150
# ╟─5bc81327-4ed8-4c78-9856-7402d370a044
# ╟─4d8cf50e-0735-4114-ac72-b50ea2565647
# ╟─eef0b9d5-32e4-4d2e-9116-7b44ce0e1567
# ╠═765aee74-6bb5-40c2-859d-c55ef6860bfa
# ╟─0a9854b8-d5e1-4b5f-8d22-328a69e9da30
# ╟─1a853d21-a144-4486-9956-6d9d5f71620f
# ╟─82f44845-5bd4-48ed-b04e-a0cf3c6ca240
# ╟─a2970366-b716-4574-b93f-49d18372e108
# ╟─d3f76adc-ba32-4862-a9c1-35d1a20e3474
# ╟─03b572a4-e51b-41e7-bc38-247644e41ebd
# ╠═4540e3ef-9db5-4fa7-bf6b-d5f8dd4244c4
# ╟─f404cd49-186b-493e-bee8-fb08c88f3f88
# ╟─e24cb37e-acfd-441c-9839-40649611c1c7
# ╠═56aaebbc-8166-4e52-a43e-93453137ad75

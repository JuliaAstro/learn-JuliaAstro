### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Astronomical Coordinates 1: Getting Started with SkyCoords.jl"
#> layout = "layout.jlhtml"
#> date = "2025-12-31"
#> description = "TODO"
#> tags = ["coordinates", "gaia"]

using Markdown
using InteractiveUtils

# ╔═╡ 62074dd9-df2c-4d8e-8d19-15204b7ddb21
begin
    # Can remove this block after Makie v0.25 is merged:
    # https://github.com/MakieOrg/Makie.jl/pull/5484
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "TOML"),
            Pkg.PackageSpec(; name = "CSV"),
            Pkg.PackageSpec(; name = "DataFramesMeta"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; name = "AstroAngles"),
            Pkg.PackageSpec(; name = "VirtualObservatory"),
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
        ]
    )

    using CSV
    using DataFramesMeta
    using CairoMakie
    using DynamicQuantities
    using SkyCoords
    using AstroAngles
    using VirtualObservatory: execute, TAPService
    using LinearAlgebra: norm

    using DynamicQuantities.Units: °
    using DynamicQuantities.Constants: pc
end

# ╔═╡ eeb03783-3d93-454f-82b8-7e55e0d803d2
begin
    using TOML: TOML
    using PlutoUI: TableOfContents
end

# ╔═╡ f2fae56d-00fc-4aea-a441-82578bba40a9
md"""
## Summary
Astronomers use a wide variety of coordinate systems and formats to represent sky coordinates of celestial objects. For example, you may have seen terms like "right ascension" and "declination" or "galactic latitude and longitude," and you may have seen angular coordinate components represented as `0h39m15.9s`, `00:39:15.9`, or `9.81625°`. The JuliaAstro packages SkyCoords.jl ([docs](https://juliaastro.org/SkyCoords/stable/)) and AstroAngles.jl ([docs](https://juliaastro.org/AstroAngles/stable/)) provide tools for representing the coordinates of objects and transforming them between different systems.

In this tutorial, we will explore how these packages can be used to work with astronomical coordinates in Julia. You may find it helpful to keep the documentation linked above open alongside this tutorial for reference or additional reading. In the text below, you may also see some links that look like ([docs](https://juliaastro.org/SkyCoords/stable/)). These links will take you to parts of the documentation that are directly relevant to the cells from which they link.

*Note: This is the 1st tutorial in a series of tutorials about coordinates in Julia.*
"""

# ╔═╡ a03f5575-1791-44ef-8722-86c525f64eb9
md"""
### Packages 📦

We start by importing some general packages that we will need below:
"""

# ╔═╡ d768eaa3-1535-457a-9455-ab65dbd59a9c
md"""
## Representing On-sky Positions with SkyCoords.jl

In Julia, the most common way of representing and working with sky coordinates is to use objects from the SkyCoords.jl package ([docs](https://juliaastro.org/SkyCoords/stable/)). A sky coordinate can be created directly from angles, as demonstrated below.

To get started, let's assume that we want to create a sky coordinate object for the center of the open star cluster NGC 188 so that later we can query and retrieve stars that might be members of the cluster. Let's also assume, for now, that we already know the sky coordinates of the cluster to be `(12.11, 85.26)` degrees in the ICRS coordinate frame. The ICRS — sometimes referred to as "equatorial" or "J2000" coordinates ([more info about the ICRS](https://arxiv.org/abs/astro-ph/0602086)) — is currently the most common astronomical coordinate frame for stellar or extragalactic astronomy. Since we already know the ICRS position of NGC 188 (see above), we can create an [`ICRSCoords`](https://juliaastro.org/SkyCoords/stable/api/#SkyCoords.ICRSCoords) object for the cluster:
"""

# ╔═╡ 3d35a242-375e-4929-83f5-2f7ab5f34490
ICRSCoords(12.11°, 85.26°)

# ╔═╡ 84d17fc3-f1f1-4d28-b818-e0872a0f58cc
md"""
In astropy, the coordinate frame is specified with a keyword argument (with ICRS as the default). In SkyCoords.jl, the frame *is* the type — every subtype of `AbstractSkyCoords` (`ICRSCoords`, `GalCoords`, `FK5Coords`, ...) names its frame explicitly — so there is no default frame to remember. As we will see later on in this series, there are many other supported coordinate frames.

In the above initialization, we passed in degree quantities from DynamicQuantities.jl ([docs](https://symbolicml.org/DynamicQuantities.jl/stable/)) to specify the angular components of our sky coordinate. Under the hood, `ICRSCoords` stores these components as unitless radians. We can also create the same object from string-formatted coordinates using the string macros from AstroAngles.jl ([docs](https://juliaastro.org/AstroAngles/stable/)). For example, if we have sexagesimal sky coordinate data where the representation of the data includes specifications of the units (the "hms" for "hour minute second", and the "dms" for "degrees minute second"):
"""

# ╔═╡ 98dd7dfe-ba52-48e4-bd58-eab412609ca6
ICRSCoords(hms"00h48m26.4s", dms"85d15m36s")

# ╔═╡ 00fdbf0e-3eb1-492f-8327-c94b578eba53
md"""
Some string representations do not explicitly define units, so it is sometimes necessary to specify these ourselves. With AstroAngles.jl, this is set by our choice of string macro: `hms"..."` interprets the string as an hour angle and `dms"..."` as degrees, regardless of the delimiters used:
"""

# ╔═╡ 4722d06a-5002-40a6-9b56-4d41ebd0a23d
ICRSCoords(hms"00:48:26.4", dms"85:15:36")

# ╔═╡ 3ecdb939-e8b2-4f62-82ef-3ae139db6512
md"""
For more information and examples on initializing these coordinate objects and parsing angle strings, see the [SkyCoords.jl](https://juliaastro.org/SkyCoords/stable/) and [AstroAngles.jl](https://juliaastro.org/AstroAngles/stable/) documentation.
"""

# ╔═╡ a8d6aff0-73e6-484a-8e1b-5e9633d43943
md"""
For the initializations demonstrated above, we assumed that we already had the coordinate component values ready. If you do not know the coordinate values and the object you are interested in is in [SIMBAD](https://simbad.u-strasbg.fr/simbad/), you can also automatically look up and load coordinate values from the name of the object. SkyCoords.jl does not ship a convenience function for this (yet!), but VirtualObservatory.jl ([docs](https://github.com/JuliaAPlavin/VirtualObservatory.jl)) makes rolling our own a one-query affair: SIMBAD exposes a Table Access Protocol (TAP) service — available in VirtualObservatory.jl as the built-in `:simbad` service — that we can query by object identifier ([more on SIMBAD TAP](https://simbad.cds.unistra.fr/simbad/sim-tap/)):

!!! todo
    Add a `from_name`-style lookup convenience function to SkyCoords.jl
"""

# ╔═╡ 81db32db-afb8-428b-8acf-b174b88fdcf9
"""
    from_name(name)

Look up the ICRS coordinates of `name` in [SIMBAD](https://simbad.u-strasbg.fr/simbad/) and return them as an `ICRSCoords` object.
"""
function from_name(name)
    result = execute(
        TAPService(:simbad),
        "SELECT ra, dec FROM basic JOIN ident ON oidref = oid WHERE id = '$name'";
        unitful = false,
    )
    isempty(result) && error("Could not resolve \"$name\" with SIMBAD")
    obj = only(result)
    return ICRSCoords(deg2rad(obj.ra), deg2rad(obj.dec))
end

# ╔═╡ 52c2de4e-96ee-407d-b7f1-0d95dae52f96
ngc188_center = from_name("NGC 188")

# ╔═╡ aa2eb71a-4033-4e8e-9a7c-7c709ba14793
md"""
The coordinate object that we defined has various ways of accessing the information contained within it. All `AbstractSkyCoords` objects have fields that allow you to retrieve the coordinate component data, but the component names will change depending on the coordinate frame of the object you have. In our examples we have created `ICRSCoords` objects, so the component names are lower-case abbreviations of Right Ascension, `.ra`, and Declination, `.dec`:
"""

# ╔═╡ 67fff99b-d4ea-4acb-9348-06c518e8cfe1
ngc188_center.ra, ngc188_center.dec

# ╔═╡ 51040904-d714-471f-817d-df42a681e19a
md"""
In astropy, these components are returned as specialized `Angle`/`Longitude`/`Latitude` objects. SkyCoords.jl takes a deliberately lightweight approach: the components are stored as plain `Float64` radians, and the surrounding package ecosystem provides the extra functionality — AstroAngles.jl for angle conversions and formatting, and DynamicQuantities.jl for attaching physical units:
"""

# ╔═╡ 4631a5e1-0214-422c-9a05-491e3ca9a92d
typeof.((ngc188_center.ra, ngc188_center.dec))

# ╔═╡ 705479be-b05c-49c3-92a8-f07e5672eab5
md"""
With the AstroAngles.jl helper functions, we can retrieve the coordinate components in different units, for example as an hour angle, in radians (identity), or in degrees:
"""

# ╔═╡ 4e83bd32-9b7a-468f-8f54-d1d42eb7dbd4
ngc188_center.ra .|> (rad2ha, identity, rad2deg)

# ╔═╡ e6be2d5c-57e0-4ced-bc1b-9c8a5bf9279e
md"""
!!! tip
    We can also work with these components as unitful quantities using DynamicQuantities.jl:
"""

# ╔═╡ 8216ff6d-5b25-43cd-8ea1-fe68bd3e8559
ngc188_center.ra * u"rad" |> us"°"

# ╔═╡ 645a72a5-0ac4-4967-8812-b0c0b5840073
md"""
We can also format the values into strings with specified delimiters and precision, for example:
"""

# ╔═╡ 4ed7215f-80ab-410c-9ca0-135066f62d40
format_angle(rad2hms(ngc188_center.ra); delim = ['ʰ', 'ᵐ', 'ˢ'])

# ╔═╡ e6868532-6497-4666-9d25-394608642da6
format_angle(rad2dms(ngc188_center.dec); delim = ['°', '′', '″'])

# ╔═╡ 98b1ba7f-c758-449c-b436-1dfe7cc05c41
md"""
## Querying the *Gaia* Archive to Retrieve Coordinates of Stars in NGC 188

Now that we have a coordinate object for the center of NGC 188, we can query many different astronomical databases directly from Julia. Here, we will use `ngc188_center` to select sources from the *Gaia* Data Release 2 catalog around the position of the center of NGC 188 to look for stars that might be members of the star cluster. Like SIMBAD above, the *Gaia* archive speaks TAP ([docs](https://www.cosmos.esa.int/web/gaia-users/archive/programmatic-access)) — VirtualObservatory.jl ships it as the built-in `:gaia` service — so we can express our cone search in ADQL (an SQL dialect for astronomy) and let `execute` handle the rest.

In the query below, the `CONTAINS`/`CIRCLE` condition selects sources within 0.5 degrees of the cluster center, `TOP 10000` limits the number of returned rows, and `ORDER BY dist` keeps the sources closest to the center if that limit kicks in. Note how we can interpolate the components of `ngc188_center` directly into the query string:
"""

# ╔═╡ 942cb618-814e-4728-9b59-6c73af5f84e1
adql_query = """
SELECT TOP 10000
    source_id, ra, dec, parallax, parallax_error, pmra, pmdec,
    radial_velocity, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
    DISTANCE(
        POINT('ICRS', ra, dec),
        POINT('ICRS', $(rad2deg(ngc188_center.ra)), $(rad2deg(ngc188_center.dec)))
    ) AS dist
FROM gaiadr2.gaia_source
WHERE 1 = CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', $(rad2deg(ngc188_center.ra)), $(rad2deg(ngc188_center.dec)), 0.5)
)
ORDER BY dist ASC
"""

# ╔═╡ 1448586d-f924-46e9-a846-798cf67cd4cf
md"""
Running this query requires an internet connection, so we have included the results file next to this notebook. The cell below only submits the query (and saves the results to `data/gaia_results.csv`) if that file is not already present:
"""

# ╔═╡ a0c9e9cb-279d-4872-9a3a-e992c0771af3
ngc188_table_all = let
    f = joinpath(@__DIR__, "data", "gaia_results.csv")
    if !isfile(f)
        CSV.write(f, DataFrame(execute(TAPService(:gaia), adql_query; unitful = false)))
    end
    CSV.read(f, DataFrame)
end

# ╔═╡ 934d74a4-4889-4496-821f-867a6e7ba9a0
# only keep stars brighter than G=19 magnitude
ngc188_table = @rsubset(ngc188_table_all, :phot_g_mean_mag < 19)

# ╔═╡ 5b052489-6bf2-45db-8c16-7d03270e95a8
nrow(ngc188_table)

# ╔═╡ a2c88449-5936-4d69-aaf5-bd17beaa0ab4
md"""
The returned `DataFrames.DataFrame` ([docs](https://dataframes.juliadata.org/stable/)) now contains about 5000 stars from *Gaia* DR2 around the coordinate position of the center of NGC 188. Let's now construct coordinate objects with the results table. In the *Gaia* data archive, the ICRS coordinates of a source are given as column names `"ra"` and `"dec"`:
"""

# ╔═╡ 34c9fdc6-49ba-49e8-8e84-eb59f2e3f380
ngc188_table.ra

# ╔═╡ 387c0990-d13a-4235-ae75-348529b418f1
ngc188_table.dec

# ╔═╡ 4b09453d-3eaa-4494-a251-2dcd70eae008
md"""
Note that the *Gaia* archive provides these columns in units of degrees ([Gaia DR2 data model](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html)), which CSV.jl reads in as plain `Float64` values. Note also that these columns contain many (nearly 5000!) coordinate values. We can attach the degree unit and broadcast the `ICRSCoords` constructor over the columns to get a single vector of coordinate objects to represent all of these coordinates:
"""

# ╔═╡ 92ede967-ba73-4714-ac5f-011ecc420663
ngc188_gaia_coords = ICRSCoords.(ngc188_table.ra * °, ngc188_table.dec * °)

# ╔═╡ 59eaf664-cd9d-4a36-858d-e5bdcc70f80c
md"""
### Exercises
"""

# ╔═╡ 6334ceba-1b0f-4c03-8d76-f6a2dfb298bb
md"""
Create a sky coordinate object for the center of the open cluster the Pleiades (either by looking up the coordinates and passing them in, or by using the convenience function we defined above):
"""

# ╔═╡ cb099a3f-3273-4fd5-8fff-e8d3eed16e47
pleiades_center = from_name("Pleiades")

# ╔═╡ 92aba2e6-1280-4f1e-98dc-fc559e5dee50
md"""
Print a string with the RA/Dec of the center of NGC 188 in the form 'HH:MM:SS.S DD:MM:SS.S'. Check your answer against [SIMBAD](http://simbad.u-strasbg.fr/simbad/), which will show you sexagesimal coordinates for the object.

(Hint: `format_angle` along with `rad2hms` and `rad2dms` might be useful)
"""

# ╔═╡ 4467ad5b-2a92-41d0-89ae-7c584fb10702
string(
    format_angle(
        rad2hms(ngc188_center.ra);
        delim = ":", digits = 1
    ),
    " ",
    format_angle(
        rad2dms(ngc188_center.dec);
        delim = ":", digits = 1
    ),
)

# ╔═╡ 35e5c8ed-4159-4be2-921d-82e8c403989f
md"""
Using a single function call, compute the angular separation between each resulting star from our *Gaia* query and the coordinates of the cluster center for NGC 188.

(Hint: [`separation`](https://juliaastro.org/SkyCoords/stable/api/#SkyCoords.separation) might be useful)
"""

# ╔═╡ db73e59e-526b-4b07-bcb8-bdea29b98c67
separation.(ngc188_gaia_coords, Ref(ngc188_center)) * u"rad" |> us"°"

# ╔═╡ 35d602ff-c4be-4f11-95dd-891f88637905
md"""
## More Than Just Sky Positions: Including Distance Information

So far, we have only used our coordinate objects to represent angular sky positions (i.e., `ra` and `dec` only). It is sometimes useful to include distance information with the sky coordinates of a source, thereby specifying the full 3D position of an object. SkyCoords.jl does not (yet) accept a `distance` keyword argument like astropy's `SkyCoord` does, but we can compose the same information ourselves: the `cartesian` function converts a coordinate object to a direction vector on the (dimensionless) unit sphere, and scaling that vector by a distance gives the full 3D position. So, if we knew that the distance to NGC 188 is 1.96 kpc, we could compute its 3D Cartesian position like so:

!!! todo
    - Add `distance` kwarg to `AbstractSkyCoords` constructors
    - Document `cartesian` and `spherical` [#84](https://github.com/JuliaAstro/SkyCoords.jl/issues/84)
"""

# ╔═╡ 9986f577-2a9a-4ade-9181-1207de773145
ngc188_center_3d = cartesian(ngc188_center).vec * 1.96e3 .* pc .|> us"Constants.pc"

# ╔═╡ 1229bdcb-e8b2-4d1f-ba38-79852c75413d
md"""
With the table of *Gaia* data we retrieved above for stars around NGC 188, `ngc188_table`, we also have parallax measurements for each star. For a precisely-measured parallax ``\varpi``, the distance ``d`` to a star can be obtained approximately as ``d \approx 1/\varpi``. This only really works if the parallax error is small relative to the parallax ([see discussion in this paper](https://arxiv.org/abs/1507.02105)), so if we want to use these parallaxes to get distances we first have to filter out stars that have low signal-to-noise parallaxes. Some sources in our table have no parallax measurement at all (the *Gaia* archive returns an empty value, which CSV.jl reads in as `missing`), so we drop those rows first:
"""

# ╔═╡ d1fd3e52-af06-4e4d-8581-b16d0e646626
ngc188_table_3d = @rsubset(
    dropmissing(ngc188_table, [:parallax, :parallax_error]),
    :parallax / :parallax_error > 10,
)

# ╔═╡ b52488ce-b09f-4f86-8286-67d5773e1f6a
nrow(ngc188_table_3d)

# ╔═╡ cbe574c4-8430-45d1-88b9-5b32f09cdba0
md"""
The above selection keeps stars that have a ~10-sigma parallax measurement, but this is an arbitrary selection threshold that you may want to tune or remove in your own use cases. This selection removed over half of the stars in our original table, but for the remaining stars we can be confident that converting the parallax measurements to distances is mostly safe.

astropy provides a specialized `Distance` class for handling common transformations of different distance representations. In Julia, this is just a short function away: the catalog of stars we queried from *Gaia* contains parallax information in milliarcsecond units, and a parallax ``\varpi`` in milliarcseconds corresponds to a distance ``d = 1000 / \varpi`` in parsecs:

!!! todo
    Add support for `mas` (milliarcsecond) angle units and a parallax-based distance convenience
"""

# ╔═╡ d74b9c88-ba7b-42bd-a6c9-cfb8d3e2fc3b
parallax_to_distance(ϖ_mas) = 1000 / ϖ_mas * pc

# ╔═╡ 454c3e7a-8f67-48f2-ae1f-8ce5211ed916
parallax_to_distance(1) |> us"Constants.pc"

# ╔═╡ 9f93b81e-8f66-4937-8bdd-2c753a7f1940
md"""
so we can compute the distances to all of our high parallax signal-to-noise stars directly from the table column:
"""

# ╔═╡ 00c6c05c-fd5d-4a37-9caf-640901b03212
gaia_dist = parallax_to_distance.(ngc188_table_3d.parallax) .|> us"Constants.pc"

# ╔═╡ d0d4a5a0-dfca-43b2-84aa-889ed378523e
md"""
We can then combine these distances with the sky positions to compute the 3D positions of all of the *Gaia* stars, just like we did for the cluster center:
"""

# ╔═╡ 46adb2c5-0efe-40e1-a5be-98812842c40c
ngc188_coords_3d = let
    directions = ICRSCoords.(ngc188_table_3d.ra * °, ngc188_table_3d.dec * °)
    [cartesian(c).vec .* d for (c, d) in zip(directions, gaia_dist)]
end

# ╔═╡ 37359964-af1b-40ae-8ac1-df279793440b
md"""
Let's now use Makie.jl to plot the sky positions of all of these sources, colored by distance to emphasize the cluster stars:
"""

# ╔═╡ ae40d35e-8a52-4061-ad16-bfdf6d1fb411
let
    fig = Figure(size = (650, 520))

    ax = Axis(
        fig[1, 1];
        xlabel = "RA [deg]",
        ylabel = "Dec [deg]",
        title = "Gaia DR2 sources near NGC 188",
    )

    plt = scatter!(
        ax, ngc188_table_3d.ra, ngc188_table_3d.dec;
        color = ustrip.(us"Constants.pc", gaia_dist) ./ 1_000,
        colorrange = (1.5, 2.5),
        colormap = :twilight,
        markersize = 5,
    )

    Colorbar(fig[1, 2], plt; label = "distance [kpc]")

    fig
end

# ╔═╡ 9e47a105-2853-4658-ac7f-e2f9bc044ccc
md"""
Now that we have 3D position information for both the cluster center and for the stars we queried from *Gaia*, we can compute the 3D separation (distance) between all of the *Gaia* sources and the cluster center with LinearAlgebra's `norm`:
"""

# ╔═╡ 98de4251-e7ab-489c-8759-5b83dc0857c6
sep3d = norm.(ngc188_coords_3d .- Ref(ngc188_center_3d)) .|> us"Constants.pc"

# ╔═╡ 7e5b4232-88da-402d-9313-98a37526a55b
md"""
### Exercises

Using the 3D separation values, define a boolean mask to select candidate members of the cluster. Select all stars within 50 pc of the cluster center. How many candidate members of NGC 188 do we have, based on their 3D positions?
"""

# ╔═╡ 38894ce9-6d1a-45da-8958-c94095e7c677
ngc188_3d_mask = sep3d .< 50pc

# ╔═╡ 388d0d5a-7e60-49b9-ad27-3013e8244cad
count(ngc188_3d_mask)

# ╔═╡ 1c478bec-a8f5-4959-9822-b643956779bc
md"""
In this tutorial, we have introduced SkyCoords.jl and AstroAngles.jl as a way to store, represent, and format astronomical sky coordinates. We used coordinate objects to parse and change coordinate representations and units. We also demonstrated how to use a coordinate object to query an astronomical database, the *Gaia* science archive, directly over TAP/ADQL with VirtualObservatory.jl. We then created a vector of coordinate objects from the queried data to represent the sky coordinates of many objects. Finally, we introduced the concept of combining Cartesian direction vectors with distances to represent the 3D position of an object or set of objects.
"""

# ╔═╡ d1eb67e4-a49c-4272-abbc-eea03e18c1ff
md"""
# Notebook setup 🔧
"""

# ╔═╡ 624fbeba-1ec9-4809-ac62-e042ac1efe67
TableOfContents(; title = "On this page", depth = 4)

# ╔═╡ 3d759c8a-ec5c-432a-9543-88a762864fa1
function frontmatter(path)
    prefix = "#> "
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ 6e367762-b07d-4c6f-92b8-090bf1aac8f3
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ f31d8083-9809-4fcb-9063-b10ae54c93ee
md"""
# Astronomical Coordinates 1: Getting Started with SkyCoords.jl

This notebook is modified from <https://learn.astropy.org/tutorials/1_Coordinates-Intro.html>

_Original authors: Adrian Price-Whelan_

!!! tip "Learning Goals"
    - Create sky coordinate objects in SkyCoords.jl using coordinate data and object names
    - Use a sky coordinate object to query the *Gaia* archive directly from Julia
    - Output coordinate data in different string representations with AstroAngles.jl
    - Demonstrate working with 3D sky coordinates (including distance information for objects)

$(keywords())
"""

# ╔═╡ Cell order:
# ╟─f31d8083-9809-4fcb-9063-b10ae54c93ee
# ╟─f2fae56d-00fc-4aea-a441-82578bba40a9
# ╟─a03f5575-1791-44ef-8722-86c525f64eb9
# ╠═62074dd9-df2c-4d8e-8d19-15204b7ddb21
# ╟─d768eaa3-1535-457a-9455-ab65dbd59a9c
# ╠═3d35a242-375e-4929-83f5-2f7ab5f34490
# ╟─84d17fc3-f1f1-4d28-b818-e0872a0f58cc
# ╠═98dd7dfe-ba52-48e4-bd58-eab412609ca6
# ╟─00fdbf0e-3eb1-492f-8327-c94b578eba53
# ╠═4722d06a-5002-40a6-9b56-4d41ebd0a23d
# ╟─3ecdb939-e8b2-4f62-82ef-3ae139db6512
# ╟─a8d6aff0-73e6-484a-8e1b-5e9633d43943
# ╠═81db32db-afb8-428b-8acf-b174b88fdcf9
# ╠═52c2de4e-96ee-407d-b7f1-0d95dae52f96
# ╟─aa2eb71a-4033-4e8e-9a7c-7c709ba14793
# ╠═67fff99b-d4ea-4acb-9348-06c518e8cfe1
# ╟─51040904-d714-471f-817d-df42a681e19a
# ╠═4631a5e1-0214-422c-9a05-491e3ca9a92d
# ╟─705479be-b05c-49c3-92a8-f07e5672eab5
# ╠═4e83bd32-9b7a-468f-8f54-d1d42eb7dbd4
# ╟─e6be2d5c-57e0-4ced-bc1b-9c8a5bf9279e
# ╠═8216ff6d-5b25-43cd-8ea1-fe68bd3e8559
# ╟─645a72a5-0ac4-4967-8812-b0c0b5840073
# ╠═4ed7215f-80ab-410c-9ca0-135066f62d40
# ╠═e6868532-6497-4666-9d25-394608642da6
# ╟─98b1ba7f-c758-449c-b436-1dfe7cc05c41
# ╠═942cb618-814e-4728-9b59-6c73af5f84e1
# ╟─1448586d-f924-46e9-a846-798cf67cd4cf
# ╠═a0c9e9cb-279d-4872-9a3a-e992c0771af3
# ╠═934d74a4-4889-4496-821f-867a6e7ba9a0
# ╠═5b052489-6bf2-45db-8c16-7d03270e95a8
# ╟─a2c88449-5936-4d69-aaf5-bd17beaa0ab4
# ╠═34c9fdc6-49ba-49e8-8e84-eb59f2e3f380
# ╠═387c0990-d13a-4235-ae75-348529b418f1
# ╟─4b09453d-3eaa-4494-a251-2dcd70eae008
# ╠═92ede967-ba73-4714-ac5f-011ecc420663
# ╟─59eaf664-cd9d-4a36-858d-e5bdcc70f80c
# ╟─6334ceba-1b0f-4c03-8d76-f6a2dfb298bb
# ╠═cb099a3f-3273-4fd5-8fff-e8d3eed16e47
# ╟─92aba2e6-1280-4f1e-98dc-fc559e5dee50
# ╠═4467ad5b-2a92-41d0-89ae-7c584fb10702
# ╟─35e5c8ed-4159-4be2-921d-82e8c403989f
# ╠═db73e59e-526b-4b07-bcb8-bdea29b98c67
# ╟─35d602ff-c4be-4f11-95dd-891f88637905
# ╠═9986f577-2a9a-4ade-9181-1207de773145
# ╟─1229bdcb-e8b2-4d1f-ba38-79852c75413d
# ╠═d1fd3e52-af06-4e4d-8581-b16d0e646626
# ╠═b52488ce-b09f-4f86-8286-67d5773e1f6a
# ╟─cbe574c4-8430-45d1-88b9-5b32f09cdba0
# ╠═d74b9c88-ba7b-42bd-a6c9-cfb8d3e2fc3b
# ╠═454c3e7a-8f67-48f2-ae1f-8ce5211ed916
# ╟─9f93b81e-8f66-4937-8bdd-2c753a7f1940
# ╠═00c6c05c-fd5d-4a37-9caf-640901b03212
# ╟─d0d4a5a0-dfca-43b2-84aa-889ed378523e
# ╠═46adb2c5-0efe-40e1-a5be-98812842c40c
# ╟─37359964-af1b-40ae-8ac1-df279793440b
# ╠═ae40d35e-8a52-4061-ad16-bfdf6d1fb411
# ╟─9e47a105-2853-4658-ac7f-e2f9bc044ccc
# ╠═98de4251-e7ab-489c-8759-5b83dc0857c6
# ╟─7e5b4232-88da-402d-9313-98a37526a55b
# ╠═38894ce9-6d1a-45da-8958-c94095e7c677
# ╠═388d0d5a-7e60-49b9-ad27-3013e8244cad
# ╟─1c478bec-a8f5-4959-9822-b643956779bc
# ╟─d1eb67e4-a49c-4272-abbc-eea03e18c1ff
# ╠═624fbeba-1ec9-4809-ac62-e042ac1efe67
# ╟─3d759c8a-ec5c-432a-9543-88a762864fa1
# ╟─6e367762-b07d-4c6f-92b8-090bf1aac8f3
# ╠═eeb03783-3d93-454f-82b8-7e55e0d803d2

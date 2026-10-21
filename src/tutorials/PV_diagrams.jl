### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> tags = ["TODO"]
#> title = "Extracting and Plotting Position-Velocity Diagrams"
#> description = "TODO"
#> date = "2025-12-31"
#> layout = "layout.jlhtml"

using Markdown
using InteractiveUtils

# ╔═╡ 8c5cedbf-c0ca-476b-96f8-a00b56ec3d96
begin
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "DynamicQuantities"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; name = "Downloads"),
            Pkg.PackageSpec(; name = "AstroAngles"),
            Pkg.PackageSpec(; name = "Interpolations"),
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
                rev = "docs-spectral",
                url = "https://github.com/JuliaAstro/AstroImages.jl",
            ),
            Pkg.PackageSpec(;
                rev = "makie-v0.25",
                url = "https://github.com/JuliaAstro/SkyCoords.jl/",
            ),
        ]
    )

    using AstroImages
    using CairoMakie
    using DynamicQuantities, SkyCoords, AstroAngles
    using Interpolations: interpolate, extrapolate, linear_interpolation, BSpline, Linear, NoInterp
    using Downloads: download

    deps_ready = true
end

# ╔═╡ 7ccb9766-6913-4259-bbfd-71aab15d9d0c
begin
    deps_ready
    using PlutoUI: TableOfContents
end

# ╔═╡ e6a87850-4abc-4e03-ac37-d5f6650bc375
using PlutoUI: details

# ╔═╡ cd1a40ec-c93f-11ef-0d74-d1e6690e3edd
md"""
# Extracting and Plotting Position-Velocity Diagrams
"""

# ╔═╡ 2b23fb29-28fd-4676-be6c-2cf7a8185ead
md"""
Tutorial modified from <https://learn.astropy.org/tutorials/PVDiagramPlotting.html>

Original authors: Adam Ginsburg, Eric Koch
"""

# ╔═╡ 796d025b-c08f-41bc-8a05-ff503d5aea99
md"""
!!! tip "Learning Goals"
    - Extract a position-velocity diagram from a spectral cube using both pixel and sky coordinates using AstroImages.jl and Interpolations.jl
    - Display the position-velocity diagram with appropriately labeled coordinates
    - Display the extraction path on the plots
"""

# ╔═╡ 01227a97-8417-4bb5-ad9e-54d99f38b3f8
md"""
!!! note "Keywords"
    pv-diagram, spectral cube, AstroImages.jl, DynamicQuantities.jl, radio astronomy, coordinates
"""

# ╔═╡ 8d49b740-1ee0-4e5d-97dd-af807d16e162
md"""
!!! warning "Summary"
    In this tutorial, we will extract position-velocity (PV) diagrams from a spectral cube and plot them.

    Position-velocity diagrams are often used in radio astronomy for analysis of rotating objects, like protostellar disks and galaxies, to measure rotation curves and determine the contained mass. They are also used in studies of atomic and molecular clouds to show where overlapping emission may point at interactions between distinct clouds. Both radio and optical position-velocity diagrams are used to study outflows, jets, and winds; in the optical, two-dimensional spectra obtained from long-slit spectrographs naturally produce the equivalent of a position-velocity diagram.
"""

# ╔═╡ bf8d1b32-8a04-443a-ac74-498e23ac6253
md"""
### Header material
"""

# ╔═╡ 016676c7-5225-48dd-a5c0-dcf19d41f68c
md"""
We import tools from several packages up front:
"""

# ╔═╡ 4953b21e-5e73-448f-8e2c-596b1ad27088
md"""
Retrieve and open a cube from astropy-data:
"""

# ╔═╡ 4b3a33de-7a53-4148-98ba-cff72d432b2e
cubefilename = download("http://www.astropy.org/astropy-data/l1448/l1448_13co.fits")

# ╔═╡ a72296e3-1c40-4545-b8a3-302ba8509e90
cube = load(cubefilename)

# ╔═╡ b458426d-27b6-419b-a072-8fe803bc5e9e
header(cube)

# ╔═╡ 3b7d4a0f-b172-4cc2-852c-397bdee4c09a
md"""
We show a single channel from the cube to visualize the data spatially. We use pixel units to start; we'll show celestial coordinates later:
"""

# ╔═╡ 015841ff-747e-421f-90ca-894ed00da6e4
cube[:, :, 26] # Default AstroImages view in a notebook (equivalent `imview`)

# ╔═╡ 11aadb2a-8671-4178-b3af-6c25df76ad69
# Automatic DimensionalData.jl axis support in plain Makie.jl
image(cube[:, :, 26], colormap = :viridis, axis = (; aspect = DataAspect()))

# ╔═╡ 593aeb81-1693-4e78-b59e-bec2bf44463c
# Full WCS axis-aware support with AstroImages.jl + Makie.jl
implot(cube[:, :, 26])

# ╔═╡ 2d86ed17-9ab9-455d-87cb-84a9f60a5ddc
# And composite AstroImages.jl plots + Makie.jl + colorbar
implotview(cube[:, :, 26])

# ╔═╡ f233606d-f868-4d69-bd3c-8e1a9c10a111
md"""
## PV Extraction from Pixel Coordinates
"""

# ╔═╡ 7fda8de2-b5b7-40c5-84c2-66f2311c04b4
md"""
First we create an _extraction path_. This `path` is the two-dimensional spatial coordinates through which we are cutting out our PV diagram. The entries are pairs of pixel coordinates, (x,y):
"""

# ╔═╡ ae81fb1a-beaf-4847-b69f-2c822d2d4056
path = [(21, 21), (41, 41), (61, 21)]

# ╔═╡ 69cfd5de-5d35-497f-b379-502f32c4d7f2
md"""
Then, we can overplot it on our original figure with the relevant conversions to WCS handled automatically:
"""

# ╔═╡ 18be1471-1cb8-4971-b6de-34e06995f72b
let
    fig, ax, p = implot(cube[:, :, 26])
    lines!(ax, path, color = "red", linewidth = 2)
    fig
end

# ╔═╡ 257ac597-57ca-48f0-8101-0f3b3321ac1b
md"""
In this case, we've selected arbitrary coordinates for our demonstration, but for your use case, you might pick the axis of an outflow, the major axis of an inclined disk, or a spiral following an arm of a galaxy.
"""

# ╔═╡ ee86b837-d696-41a1-9a35-35de954cd446
md"""
We can then extract the PV diagram with the `extract_pv` function, specifying the separation between points in pixels that we would like to sample by passing the `spacing` keyword:
"""

# ╔═╡ 8c9358f3-be8e-49e3-b031-4fd04dcb289d
md"""
and plot it. pvdiagram is a plain (offset ``\times`` velocity) matrix, so we plot it against offset and velocity axes built from the cube's WCS:

!!! todo
    Upstream this
"""

# ╔═╡ a8c1bc08-78fb-486b-9b38-246e54ecf2b9
md"""
Changing units to the more commonly used km/s and more readable arcminutes can be done with [DynamicQuantities.jl](https://github.com/JuliaPhysics/DynamicQuantities.jl), by attaching units to each axis and stripping them off again in the units we want to display:
"""

# ╔═╡ 619c1f68-e5a7-409f-ac69-d8a09082ba41
md"""
We can put all this together next to a map of the peak intensity along each line of sight, which we will use for display from here on:
"""

# ╔═╡ fe6b0f2d-5cae-4474-9a13-818511726c80
peak = dropdims(maximum(cube; dims = 3); dims = 3)

# ╔═╡ 207d1ea0-cb19-4355-865e-50abd0feb5c2
md"""
We wrap this figure in a function so that we can reuse it in the sections below:
"""

# ╔═╡ 7d361163-62f9-4d14-b085-f457043aeba9
"""
    pvfigure(path, pvdiagram; spacing = 1, aspect = nothing)

Show the peak intensity map with the pixel `path` drawn on it, next to `pvdiagram`
in units of arcmin and km/s.
"""
function pvfigure(path, pvdiagram; spacing = 1, aspect = nothing)
    velocities = [
        pixel_to_world(cube, [1, 1, k])[3] for k in axes(pvdiagram, 2)
    ] * u"m/s"
    offsets = (axes(pvdiagram, 1) .- 0.5) .* spacing .* abs(cube["CDELT2"]) * u"deg"

    fig = Figure(size = (900, 400))

    # Collapsing the cube leaves a reference velocity behind, so we set our own title
    ax_map, _ = implot(fig[1, 1], peak; cmap = :viridis, axis = (; title = "Peak intensity"))
    lines!(ax_map, path, color = "red", linewidth = 2)

    ax_pv = Axis(fig[1, 2]; ylabel = "Velocity", xlabel = "Offset", aspect)
    hm = heatmap!(ax_pv, offsets .|> us"arcmin", velocities .|> us"km/s", pvdiagram, colormap = :viridis)
    Colorbar(fig[1, 3], hm, label = "Brightness temperature [K]")

    return fig
end

# ╔═╡ bf5a0a11-8332-4d1d-a3df-a9e8b395d520
md"""
## PV Extraction from Sky Coordinates
"""

# ╔═╡ 3fbfcc8f-0abd-4501-9805-e572bbfc428b
md"""
We can also make paths by supplying [SkyCoords.jl](https://juliaastro.org/SkyCoords/stable/) coordinates instead:
"""

# ╔═╡ 9df09b72-d577-4962-ac6d-e2982ae38298
const ° = u"deg"

# ╔═╡ 19affd33-ba8d-4ef8-833c-32312e3336f7
skypath = [
    FK5Coords{2000}(ha2deg(3.4)°, 30.5°),
    FK5Coords{2000}(ha2deg(3.43)°, 30.75°),
    FK5Coords{2000}(ha2deg(3.42)°, 30.5°),
]

# ╔═╡ 811d9490-dba9-4377-bb0c-3de58051b90c
md"""
To plot it, we first convert the path to the frame of the cube's WCS. Its header has no `RADESYS` keyword, so this defaults to ICRS:
"""

# ╔═╡ 5d16537a-c733-4557-8ce6-872a88142d79
skypath_icrs = convert.(ICRSCoords, skypath)

# ╔═╡ 82a87989-039f-47ec-83e7-8d243e6f5a38
md"""
AstroImages.jl's `world_transform` then maps these sky coordinates (in degrees) onto the image for us:
"""

# ╔═╡ 6f19ad03-7895-4ab2-a028-c03ef13874fc
let
    fig, ax, p = implot(cube[:, :, 26])
    radec = [rad2deg.((c.ra, c.dec)) for c in skypath_icrs]
    lines!(
        ax, radec;
        color = "red", linewidth = 2, transformation = world_transform(p)
    )
    fig
end

# ╔═╡ 95398ade-ed1f-4e59-bdf5-4813283064fa
md"""
To extract the PV diagram, we convert the path to pixel coordinates with the cube's WCS:
"""

# ╔═╡ 07bd38a1-e2af-4652-8b25-662e35f1da93
skypath_pixels = map(skypath_icrs) do c
    x, y, _ = world_to_pixel(cube, [rad2deg(c.ra), rad2deg(c.dec), 0]) # Any velocity will do
    return (x, y)
end

# ╔═╡ cab0d9c5-19d0-42ee-841e-d42ed2c57fdc
md"""
and pass them to `extract_pv` like before:
"""

# ╔═╡ 8bfb5636-301a-410e-a335-44a433a3d1d0
md"""
Putting it all together again:
"""

# ╔═╡ 0518e78f-3e8c-43bc-b571-7507d25b0866
md"""
We can also change the aspect ratio of the PV diagram. `Figure(size = ...)` controls the overall figure size, and the `aspect` attribute of the PV diagram's `Axis` controls its shape. For example, we can apply a square aspect with proper colorbar spacing with the following:
"""

# ╔═╡ c8d91b03-964f-4c26-9799-df4f3c668e99
md"""
See the [GridLayout > Aspect](https://docs.makie.org/stable/reference/blocks/gridlayout#Aspect) section of the Makie.jl docs for more.
"""

# ╔═╡ ed6647b9-4c83-4fbf-918c-46d40eec6d38
md"""
## PV Extraction with Spatial Averaging
"""

# ╔═╡ d7d3f889-ca88-4a6e-906d-f9880c8c2992
md"""
pvextractor.Path allows you to specify a width to average over, which specifies a spatial range around the path to average over.

In other words, it turns our path into a series of consecutive rectangular regions.

!!! todo
    Maybe pull in Photometry.jl for this
"""

# ╔═╡ 243a86ee-6576-4591-8f06-c1d10e551e0f


# ╔═╡ cee07c57-b54b-42ff-9e54-b9671a85b921


# ╔═╡ 6b71ca84-d5b4-4aad-8399-ea59986ded3f
md"""
We can plot this path as a set of patches to show where we averaged. The default spacing is 1 pixel,so we plot 1-pixel chunks.
"""

# ╔═╡ dfb786d6-5ba3-401f-8faf-a23e656dbc89
# See TODO in "Notebook setup" section

# ╔═╡ 38c3da60-70f1-4b35-b41d-78a076b507ba
md"""
We can also have more widely spaced chunks.

Note that the spacing given to extract_pv_slice affects the shape of the output PV diagram, so we also change the aspect ratio:
"""

# ╔═╡ d41a3253-4369-4a8d-ae5b-5acf698224da
md"""
## Saving
"""

# ╔═╡ b8c5e55c-cdae-466e-aaaa-666d770e6101
md"""
Finally, we can save the extracted PV diagram as a FITS file:
"""

# ╔═╡ 66e6e8ba-41c6-47fc-bd39-0fb35711e4ad
# See TODO in "Notebook setup" section

# ╔═╡ f3705288-c780-4f9a-abf3-3fcc5b35715a
md"""
We can also save the figure as a png or pdf:
"""

# ╔═╡ 877fc6d0-6fc9-4fa8-a681-20a36553bf78
md"""
# Notebook setup
"""

# ╔═╡ e7a62f39-e053-4cc0-afe8-eb06869a872c
TableOfContents()

# ╔═╡ e30ac0e5-40b6-4ee7-ab8e-8a44eb6c3a04
md"""
!!! todo
    Upstream path sampling and PV extraction, e.g. into AstroImages.jl next to its axis-aligned position-velocity slices. These helpers stand in for pvextractor's `Path.sample_points` and `extract_pv_slice`, but interpolate linearly, while pvextractor defaults to cubic interpolation (`order = 3`).

    $(details("Rough plan for an AstroImages.jl PR", md\"\"\"
    1. `src/pv.jl` (new)
        - API: `pvslice(cube, path; spacing = 1)` 
        - returns an AstroImage with a PV header, not a bare matrix. The
        header would carry:
        - `CTYPE1 = 'OFFSET'` and `CDELT1 = spacing` `\times` pixel scale.
        - The cube's spectral axis keywords moved to axis 2 (CTYPE3/CRVAL3/CDELT3/CRPIX3/CUNIT3 become …2, plus SPECSYS).
        - That needs the header and WCS code already in core: header, copyheader, world_to_pixel, and the parent-pixel handling in `src/wcs.jl`.
        - What this then gives us here:
            - implotview(pv) draws WCS-labelled offset/velocity axes.
            - save("pv.fits", pv) writes a valid FITS file through io.jl.
        - Sky-coordinate paths fit here too: accept world-coordinate tuples and convert them with the existing `world_to_pixel(img, ...)`. That avoids a SkyCoords dependency; SkyCoords users pass rad2deg.(lonlat(c)).
    
    1. Interpolation: ext/AstroImagesInterpolationsExt.jl (new)
        - Follow the pattern the package already uses for Makie (implot) and, on the reproject branch, Reproject
        (ext/AstroImagesReprojectExt.jl):
        - Declare `function pvslice end` in core.
    
    1. `src/wcs-grid.jl`: a small fix needed for nice PV plots
        - `isangular` treats any axis with `CUNIT = 'deg'` as angular.
        - So `wcslabels` would label an OFFSET axis in degrees-minutes-seconds, like a declination.
        - Fix: use sexagesimal labels only for celestial `CTYPEs`, and give `OFFSET` plain decimal (or arcmin) ticks.
        - Optionally teach `ctype_label` to print "Offset".
    
    1. Docs: `docs/src/manual/spec.md` § "Position–Velocity Slices"
        - That section currently shows only the straight cut `HIcube[Y = 45]`.
        - Add an "Along an arbitrary path" subsection there.
    
    1. Tests: new `test/pv.jl` that checks for
        - A linear-ramp cube, where bilinear sampling is exact.
        - A path running along a row `(Y = k)`, which must reproduce the existing `cube[Y = k]` slice.
    
    1. Later: width averaging
        - `ext/AstroImagesPhotometryExt.jl` (new), wrapping Photometry.jl's rotated `RectangularAperture` for `pvslice(...; width)`.
    "\"\"))
    """

# ╔═╡ 3996ac46-9104-4745-ac70-0916f50527b4
begin
    """
        sample_path(path; spacing = 1)

    Midpoints of `spacing`-pixel steps along the polyline through `path`,
    like pvextractor's `Path.sample_points`.
    """
    function sample_path(path; spacing = 1)
        xs, ys = float.(first.(path)), float.(last.(path))
        d = [0; cumsum(hypot.(diff(xs), diff(ys)))] # Distance along the path
        n = floor(Int, d[end] / spacing)
        edges = range(0, n * spacing; length = n + 1)
        xe, ye = linear_interpolation(d, xs).(edges), linear_interpolation(d, ys).(edges)
        return (xe[1:(end - 1)] .+ xe[2:end]) ./ 2, (ye[1:(end - 1)] .+ ye[2:end]) ./ 2
    end

    """
        extract_pv(cube, path; spacing = 1)

    Bilinearly sample every channel of `cube` along `path`, returning an
    (offset × velocity) matrix.
    """
    function extract_pv(cube, path; spacing = 1)
        xs, ys = sample_path(path; spacing)
        itp = extrapolate(interpolate(parent(cube), (BSpline(Linear()), BSpline(Linear()), NoInterp())), NaN)
        return [itp(x, y, k) for (x, y) in zip(xs, ys), k in axes(parent(cube), 3)]
    end
end

# ╔═╡ a1b785c7-daa4-47a8-a062-a4aa8bf4acb2
pvdiagram = extract_pv(cube, path; spacing = 1)

# ╔═╡ 2d2e1bf5-5a91-4142-9a96-361db70e2b88
let
    velocities = [pixel_to_world(cube, [1, 1, k])[3] for k in axes(pvdiagram, 2)] # m/s
    offsets = (axes(pvdiagram, 1) .- 0.5) .* abs(cube["CDELT2"]) # deg, for spacing = 1 pixel

    fig = Figure()
    ax = Axis(fig[1, 1], ylabel = "velocity [m/s]", xlabel = "offset [°]")
    hm = heatmap!(ax, offsets, velocities, pvdiagram, colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = "Brightness temperature [K]")
    fig
end

# ╔═╡ a92ecb45-1dfb-48d8-819e-670d343e7db4
let
    velocities = [
        pixel_to_world(cube, [1, 1, k])[3] for k in axes(pvdiagram, 2)
    ] * u"m/s"
    offsets = (axes(pvdiagram, 1) .- 0.5) .* abs(cube["CDELT2"]) * u"deg"

    fig = Figure()
    ax = Axis(fig[1, 1], ylabel = "Velocity", xlabel = "Offset")
    hm = heatmap!(ax, offsets .|> us"arcmin", velocities .|> us"km/s", pvdiagram, colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = "Brightness temperature [K]")
    fig
end

# ╔═╡ ec25be6f-6037-47f3-8a71-e06d2710ecaf
pvfigure(path, pvdiagram)

# ╔═╡ 7314d89b-7f63-4db4-a263-f47d98b9d6f8
pvdiagram2 = extract_pv(cube, skypath_pixels)

# ╔═╡ 8032f5d1-ed0b-4e1c-b0c8-02da7b86a2de
pvfigure(skypath_pixels, pvdiagram2)

# ╔═╡ 47e84d41-df7e-481b-8953-e52ba92bbb8e
let
    fig = pvfigure(skypath_pixels, pvdiagram2; aspect = AxisAspect(1))
    colsize!(fig.layout, 2, Aspect(1, 1))
    fig
end

# ╔═╡ Cell order:
# ╟─cd1a40ec-c93f-11ef-0d74-d1e6690e3edd
# ╟─2b23fb29-28fd-4676-be6c-2cf7a8185ead
# ╟─796d025b-c08f-41bc-8a05-ff503d5aea99
# ╟─01227a97-8417-4bb5-ad9e-54d99f38b3f8
# ╟─8d49b740-1ee0-4e5d-97dd-af807d16e162
# ╟─bf8d1b32-8a04-443a-ac74-498e23ac6253
# ╟─016676c7-5225-48dd-a5c0-dcf19d41f68c
# ╠═8c5cedbf-c0ca-476b-96f8-a00b56ec3d96
# ╟─4953b21e-5e73-448f-8e2c-596b1ad27088
# ╠═4b3a33de-7a53-4148-98ba-cff72d432b2e
# ╠═a72296e3-1c40-4545-b8a3-302ba8509e90
# ╠═b458426d-27b6-419b-a072-8fe803bc5e9e
# ╟─3b7d4a0f-b172-4cc2-852c-397bdee4c09a
# ╠═015841ff-747e-421f-90ca-894ed00da6e4
# ╠═11aadb2a-8671-4178-b3af-6c25df76ad69
# ╠═593aeb81-1693-4e78-b59e-bec2bf44463c
# ╠═2d86ed17-9ab9-455d-87cb-84a9f60a5ddc
# ╟─f233606d-f868-4d69-bd3c-8e1a9c10a111
# ╟─7fda8de2-b5b7-40c5-84c2-66f2311c04b4
# ╠═ae81fb1a-beaf-4847-b69f-2c822d2d4056
# ╟─69cfd5de-5d35-497f-b379-502f32c4d7f2
# ╠═18be1471-1cb8-4971-b6de-34e06995f72b
# ╟─257ac597-57ca-48f0-8101-0f3b3321ac1b
# ╟─ee86b837-d696-41a1-9a35-35de954cd446
# ╠═a1b785c7-daa4-47a8-a062-a4aa8bf4acb2
# ╟─8c9358f3-be8e-49e3-b031-4fd04dcb289d
# ╠═2d2e1bf5-5a91-4142-9a96-361db70e2b88
# ╟─a8c1bc08-78fb-486b-9b38-246e54ecf2b9
# ╠═a92ecb45-1dfb-48d8-819e-670d343e7db4
# ╟─619c1f68-e5a7-409f-ac69-d8a09082ba41
# ╠═fe6b0f2d-5cae-4474-9a13-818511726c80
# ╟─207d1ea0-cb19-4355-865e-50abd0feb5c2
# ╠═7d361163-62f9-4d14-b085-f457043aeba9
# ╠═ec25be6f-6037-47f3-8a71-e06d2710ecaf
# ╟─bf5a0a11-8332-4d1d-a3df-a9e8b395d520
# ╟─3fbfcc8f-0abd-4501-9805-e572bbfc428b
# ╠═9df09b72-d577-4962-ac6d-e2982ae38298
# ╠═19affd33-ba8d-4ef8-833c-32312e3336f7
# ╟─811d9490-dba9-4377-bb0c-3de58051b90c
# ╠═5d16537a-c733-4557-8ce6-872a88142d79
# ╟─82a87989-039f-47ec-83e7-8d243e6f5a38
# ╠═6f19ad03-7895-4ab2-a028-c03ef13874fc
# ╟─95398ade-ed1f-4e59-bdf5-4813283064fa
# ╠═07bd38a1-e2af-4652-8b25-662e35f1da93
# ╟─cab0d9c5-19d0-42ee-841e-d42ed2c57fdc
# ╠═7314d89b-7f63-4db4-a263-f47d98b9d6f8
# ╟─8bfb5636-301a-410e-a335-44a433a3d1d0
# ╠═8032f5d1-ed0b-4e1c-b0c8-02da7b86a2de
# ╟─0518e78f-3e8c-43bc-b571-7507d25b0866
# ╠═47e84d41-df7e-481b-8953-e52ba92bbb8e
# ╟─c8d91b03-964f-4c26-9799-df4f3c668e99
# ╟─ed6647b9-4c83-4fbf-918c-46d40eec6d38
# ╟─d7d3f889-ca88-4a6e-906d-f9880c8c2992
# ╠═243a86ee-6576-4591-8f06-c1d10e551e0f
# ╠═cee07c57-b54b-42ff-9e54-b9671a85b921
# ╟─6b71ca84-d5b4-4aad-8399-ea59986ded3f
# ╠═dfb786d6-5ba3-401f-8faf-a23e656dbc89
# ╟─38c3da60-70f1-4b35-b41d-78a076b507ba
# ╟─d41a3253-4369-4a8d-ae5b-5acf698224da
# ╟─b8c5e55c-cdae-466e-aaaa-666d770e6101
# ╠═66e6e8ba-41c6-47fc-bd39-0fb35711e4ad
# ╟─f3705288-c780-4f9a-abf3-3fcc5b35715a
# ╟─877fc6d0-6fc9-4fa8-a681-20a36553bf78
# ╠═e7a62f39-e053-4cc0-afe8-eb06869a872c
# ╠═7ccb9766-6913-4259-bbfd-71aab15d9d0c
# ╟─e30ac0e5-40b6-4ee7-ab8e-8a44eb6c3a04
# ╠═e6a87850-4abc-4e03-ac37-d5f6650bc375
# ╠═3996ac46-9104-4745-ac70-0916f50527b4

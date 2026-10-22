# The sidebar lists the tags under these headings, in this order (each group sorted A–Z).
# A tag used by a tutorial but missing here is listed under "Other".
tag_groups = [
    "Dimensions" => ["coordinates", "units", "WCS"],
    "Data" => ["ASDF", "data cubes", "file I/O", "FITS", "query", "tables"],
    "Domains" => ["astrostatistics", "cosmology", "extinction", "observational astronomy", "physics", "radio astronomy"],
    "Methods" => ["Bayesian", "MCMC", "models", "optimization", "photometry", "spectroscopy"],
    "Data Viz" => [
        "colorbars", "contours", "corner plots", "errorbars", "heatmaps", "histograms",
        "images", "line plots", "scatter plots", "stairs", "uncertainty bands",
    ],
    "Surveys & Missions" => ["Gaia", "SDSS"],
]

# In a `let`, so these don't overwrite the page's own frontmatter variables (e.g. `tags`) in the template module.
let
    # Every notebook under `tutorials/` (layouts only see pages through `collections`, so it needs at least one tag),
    # sorted A–Z by its frontmatter `title`.
    tutorials = unique(p -> p.url, [p for c in values(collections) for p in c.pages if startswith(p.input.relative_path, "tutorials/")])
    sort!(tutorials; by = p -> lowercase(get(p.output.frontmatter, "title", basename(p.input.relative_path))))

    used = Set(tag for p in tutorials for tag in get(p.output.frontmatter, "tags", String[]))
    grouped = Set(tag for (_, tags) in tag_groups for tag in tags)
    groups = [heading => sort!(filter(in(used), tags); by = lowercase) for (heading, tags) in tag_groups]
    push!(groups, "Other" => sort!(collect(setdiff(used, grouped)); by = lowercase))

    Dict(
        "tutorials" => tutorials,
        # Only tags some tutorial uses, so no sidebar link leads to an empty gallery.
        "tag_groups" => filter(group -> !isempty(last(group)), groups),
    )
end

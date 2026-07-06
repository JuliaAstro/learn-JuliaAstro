### A Pluto.jl notebook ###
# v1.0.1

#> [frontmatter]
#> image = "https://roman-docs.stsci.edu/_/0A8064E70182A738F6DB006B4CF16B4C/1763466662461/assets/img/roman-logo-white.png"
#> order = 3
#> title = "ASDF files"
#> layout = "layout.jlhtml"
#> date = "2025-12-14"
#> description = "View and manipulate data from ASDF files."
#> tags = ["file I/O", "ASDF"]

using Markdown
using InteractiveUtils

# ╔═╡ eb241a66-4965-4f11-a67b-4d344b71c4d1
md"""
# Working with ASDF files

!!! tip "Learning goals"
    - Load a simple ASDF file
    - Edit its contents
    - Save modified data to a new file

$(keywords())

!!! warning "Companion content"
    - <https://learn.juliaastro.org/tutorials/fileio-fits_tables/>
    - <https://learn.juliaastro.org/tutorials/fileio-fits_images/>
"""

# ╔═╡ 3a4c5c6c-00eb-4f73-be50-81848e5b5e0e
md"""
## Summary

We will show how to use the Julia implementation of the [ASDF spec](https://www.asdf-format.org/en/latest/overview.html) to load and manipulate ASDF files.

!!! note
    For more examples, see the [ASDF.jl documentation](juliaastro.org/ASDF/stable/).
"""

# ╔═╡ 0bd0cd34-8cc9-4edb-abd3-d9539c4972cd
md"""
### Packages 📦
"""

# ╔═╡ 959c66f6-d8df-11f0-8c7e-1150993431da
begin
    using ASDF: load, save
    using Downloads: download
end

# ╔═╡ 031f2756-f239-40aa-98f3-9718e6c493f0
md"""
## Load
"""

# ╔═╡ 3711d189-47f3-44d5-bf5d-d0886bb93ed1
url = "https://github.com/JuliaAstro/ASDF.jl/blob/main/test/data/asdf-1.6.0/scalars.asdf?raw=true"

# ╔═╡ 81007f0c-85d1-462d-b753-1d4619a6c9ae
fname = download(url)

# ╔═╡ b69dd6af-448e-4a71-aad2-3c2053c45843
af = (load ∘ download)(url)

# ╔═╡ 95dbe901-b4d4-4226-8087-8300639753b8
md"""
## Modify
"""

# ╔═╡ c1718f79-9bcc-4e90-ba37-80c8f4545d33
af["char"] = 'c'

# ╔═╡ 90f502cc-bef1-40bb-8730-bf0c4bf9d3df
af

# ╔═╡ 8331ac41-8508-42a0-8ee9-7d1d665afd0d
md"""
## Save
"""

# ╔═╡ 786d762f-074e-4f84-88ec-784ca2b9f12e
save("my-file.asdf", af)

# ╔═╡ 209b4d0d-7222-43b8-9a50-736f60ba5dad
md"""
# Notebook setup 🔧
"""

# ╔═╡ 1007d00d-16de-497d-b584-9154af0ca155
TableOfContents(; depth = 4)

# ╔═╡ 5dbb58fd-a635-424c-9ddf-6ac833b971bb
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ 410167f1-5f1f-4806-8121-2f2d480d5730
begin
    using Pluto: frontmatter
    using PlutoUI: TableOfContents
    using Test: @test
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ASDF = "686f71d1-807d-59a4-a860-28280ea06d7b"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
Pluto = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
"""

# ╔═╡ Cell order:
# ╟─eb241a66-4965-4f11-a67b-4d344b71c4d1
# ╟─3a4c5c6c-00eb-4f73-be50-81848e5b5e0e
# ╟─0bd0cd34-8cc9-4edb-abd3-d9539c4972cd
# ╠═959c66f6-d8df-11f0-8c7e-1150993431da
# ╟─031f2756-f239-40aa-98f3-9718e6c493f0
# ╠═3711d189-47f3-44d5-bf5d-d0886bb93ed1
# ╠═81007f0c-85d1-462d-b753-1d4619a6c9ae
# ╠═b69dd6af-448e-4a71-aad2-3c2053c45843
# ╟─95dbe901-b4d4-4226-8087-8300639753b8
# ╠═c1718f79-9bcc-4e90-ba37-80c8f4545d33
# ╠═90f502cc-bef1-40bb-8730-bf0c4bf9d3df
# ╟─8331ac41-8508-42a0-8ee9-7d1d665afd0d
# ╠═786d762f-074e-4f84-88ec-784ca2b9f12e
# ╟─209b4d0d-7222-43b8-9a50-736f60ba5dad
# ╠═1007d00d-16de-497d-b584-9154af0ca155
# ╟─5dbb58fd-a635-424c-9ddf-6ac833b971bb
# ╠═410167f1-5f1f-4806-8121-2f2d480d5730
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

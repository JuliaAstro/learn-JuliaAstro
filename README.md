# learn.JuliaAstro

> [!NOTE]
> **In progress:**
> - Makie.jl: units in axes labels (and more!) <https://github.com/MakieOrg/Makie.jl/pull/5323>
> - DimensionalData.jl: Update [Makie.jl extension](https://github.com/rafaqz/DimensionalData.jl/blob/main/ext/DimensionalDataMakieExt.jl) for <https://github.com/MakieOrg/Makie.jl/pull/5323>
> - GLM.jl: Weighted confidence intervals <https://github.com/JuliaStats/GLM.jl/pull/487>

> [!WARNING]
> **To be implemented:**
> - DimensionalData.jl: Update [Makie.jl extension](https://github.com/rafaqz/DimensionalData.jl/blob/main/ext/DimensionalDataMakieExt.jl) for <https://github.com/MakieOrg/Makie.jl/pull/5323>
> - DynamicQuantities.jl: `Base.sum` support <https://github.com/JuliaPhysics/DynamicQuantities.jl/issues/76#issuecomment-3614719247>
> - GLM.jl
>     - Float32 support: <https://github.com/JuliaStats/GLM.jl/issues/260>
>     - Rename `wts` to `weights`: <https://github.com/JuliaStats/GLM.jl/pull/570>

## Development

In root directory:

```julia-repl
julia --proj

using PlutoPages: PlutoPages

PlutoPages.develop(; input_dir = abspath("src"), output_dir = abspath("build"), cache_dir = abspath("_cache/"))
```

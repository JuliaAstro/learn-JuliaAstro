# learn.JuliaAstro

## Development

In top-level of activated project directory:

```julia-repl
julia> using PlutoPages

julia> PlutoPages.develop(; input_dir = abspath("src"), output_dir = abspath("build"), cache_dir = abspath("_cache/"))
```

To update packages, including those on development branches, run:

```julia-repl
julia> Pluto.activate_notebook_environment(path) do Pkg.update() end
```

This may be turned into a cron job at some point.

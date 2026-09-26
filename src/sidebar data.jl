sections = ["file I/O", "units", "coordinates", "models", "cosmology", "dust", "spectroscopy"]

Dict(
    "main" => [uppercase(section) => collections[section].pages for section in sections],
)

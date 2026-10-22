// Filter the homepage gallery (_includes/gallery.jlmd) to the cards that have every tag in `?tag=a,b`.
// A classic (non-module) script placed after the cards, so it runs before the first paint and filtered-out cards never flash.
{
    const selected = (new URLSearchParams(window.location.search).get("tag") ?? "").split(",").filter(Boolean)

    let shown = 0
    for (const card of document.querySelectorAll(".pp-card")) {
        const tags = card.dataset.tags.split(",")
        card.hidden = !selected.every((tag) => tags.includes(tag))
        if (!card.hidden) shown++
    }
    document.querySelector(".pp-gallery-empty").hidden = shown > 0
    // The sidebar's "Clear all tags" link (hidden on every other page).
    document.querySelector(".pp-clear-tags").hidden = selected.length === 0

    // Each tag link in the sidebar is highlighted when its tag is selected,
    // and clicking it adds the tag to, or removes it from, the current selection.
    for (const a of document.querySelectorAll("a[data-tag]")) {
        const tag = a.dataset.tag
        const active = selected.includes(tag)
        const next = active ? selected.filter((t) => t !== tag) : [...selected, tag]
        a.classList.toggle("active", active)
        // Literal commas between the encoded tags (URLSearchParams would write them as `%2C`).
        a.href = next.length ? "?tag=" + next.map(encodeURIComponent).join(",") : "?"
    }
}

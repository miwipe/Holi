# Development

The `holi_dev` branch reorganises documentation and the literature explorer while retaining existing analysis entry points and scientific figure files. Changes to database layouts or biological filtering criteria should be reviewed separately from presentation changes.

## Folder conventions

| Folder | Contents |
| --- | --- |
| `bin/` | Existing analysis and postprocessing scripts; keep entry points stable |
| `envs/` | Conda environments; `bin/*.yml` links preserve older commands |
| `config/` | Portable examples; `*.local.*` files are ignored |
| `literature/R/` | Literature plotting project and exploratory `sandbox.R` |
| `data/literature/` | Canonical dashboard snapshot and provenance |
| `dashboard/` | Browser code, styles, generated snapshot wrapper, pinned vendor assets |
| `dashboards/reference-catalog/` | Existing catalogue dashboard, retained unchanged |
| `docs/wiki/` | Current proposed wiki source pages |
| `docs/wiki-archive/` | Historical wiki snapshot, preserved rather than treated as current instructions |
| `figures/` | Existing original plots, workflow, and artwork |
| `scripts/` | Build and validation utilities |

The old literature directory points to `literature/`, and `literature/sedaDNA_lit_map` points to `R/`. The original reference-catalogue HTML URL redirects to its new folder and preserves query/hash navigation through JavaScript. Use the new paths for future work.

## Preview the literature explorer

Node.js is required for the small build/validation scripts; the dashboard itself uses ordinary HTML, CSS, and JavaScript and needs no package installation.

```bash
node scripts/build-dashboard.mjs
node scripts/validate-dashboard.mjs
python3 -m http.server 8765 --bind 127.0.0.1
```

Open `http://127.0.0.1:8765/`. The Plotly distribution and world geometry are stored locally. Serve over HTTP so the browser can load map geometry; a `file://` opening is not the supported map preview. The root `index.html` remains the GitHub Pages entry point.

`node scripts/build-dashboard.mjs --standalone` also writes `dist/holi-literature.html` with scripts, styles, and snapshot embedded. This is a convenience export, not a claim of a fully offline package: the map geometry, original figure gallery, and reference catalogue still need the repository assets. The normal multi-file dashboard is the maintained version.

## Update the literature snapshot

1. Review the source sheet and proposed DOI corrections with the curator.
2. Export the `Eukaryotes` data tab, not the DOI review tab or a working list.
3. Preserve source values in `data/literature/eukaryotes.snapshot.json` and record capture time, source URL, and row count in `provenance.json`.
4. Run the build and validation commands above. Review count changes and source-table semantics before committing.

`dashboard/snapshot.js` is generated so that loading the initial table does not require a fetch. Edit the JSON and rebuild rather than editing the generated wrapper. The current snapshot is unchanged from the previous explorer; the original capture time is unknown. Static PNGs may predate it.

Live refresh and local table import are explicit browser actions. They do not rewrite the committed snapshot. The dashboard preserves rows without coordinates for statistics and the source table; only the map excludes them. Unknown coordinates are not swapped automatically. Publication identity normalises DOI case and standard DOI URL prefixes without adopting unreviewed corrections.

## Documentation and wiki publication

Edit the pages in `docs/wiki/`, then check local links:

```bash
python3 scripts/validate-docs.py
python3 scripts/prepare-wiki.py /tmp/holi-wiki-export
```

The export rewrites local page links for GitHub's wiki and repository-document links for the main GitHub repository. It records the intended target branch in those repository URLs. Review the exported pages before publishing.

The live wiki is a separate repository, so publication is an explicit maintainer step after review:

```bash
git clone https://github.com/miwipe/Holi.wiki.git /tmp/Holi.wiki
cp /tmp/holi-wiki-export/*.md /tmp/Holi.wiki/
git -C /tmp/Holi.wiki diff
# Review the diff, then commit and push in the wiki repository when approved.
```

The export does not delete existing wiki pages or change the live wiki. Imported historical content remains available in `docs/wiki-archive/`. The `_Sidebar.md` and `_Footer.md` files in the proposed wiki provide navigation and version context.

## Validation boundaries

The validation scripts check data parsing, DOI deduplication, coordinate handling, dashboard assets, local documentation links, compatibility paths, and preservation of the existing plots. R and Bash sources should also pass parsing/syntax checks. These checks do not run the high-memory biological pipeline, install environments, validate database contents, or establish scientific equivalence between R and browser filtering.

Review the desktop and narrow-screen dashboard views, search and empty results, filters/reset, table pagination, and CSV/PNG downloads before deployment. Deploying GitHub Pages or publishing the live wiki is separate from preparing this branch.

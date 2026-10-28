# Literature and figures

The [Google Sheet](https://docs.google.com/spreadsheets/d/13cmBUi4cigUaTKtQeFLFvS0gXT8AeWxWKzHv2UcOBCI/edit#gid=0) is the curated source. The `Eukaryotes` tab supplies the original plots and literature dashboard. Other workbook tabs contain microbial studies, additions, and working notes.

## Scope and counts

The eukaryotic compilation covers ancient environmental DNA analysed with shotgun sequencing or capture enrichment, excluding studies using only metabarcoding. It extends the [sedimentary aeDNA inventory](https://doi.org/10.5281/zenodo.6847521) beyond sediments to additional environmental and archaeological materials.

A table row is a **site record**, not necessarily a unique physical site. A paper can have many rows. Interactive publication counts use case-insensitive DOI identifiers, falling back to citation text only if no DOI is supplied. Papers can appear in more than one method or database category; totals across categories need not equal the headline publication count.

The [DOI review tab](https://docs.google.com/spreadsheets/d/13cmBUi4cigUaTKtQeFLFvS0gXT8AeWxWKzHv2UcOBCI/edit#gid=1498244963) contains identifiers pending curator review. No proposed correction has been applied automatically to the dashboard snapshot or committed figures.

## Original R figures

```bash
Rscript literature/R/generate_figures.R
```

Run from the repository root, or open `literature/R/sedaDNA_litterature.Rproj`. The plotting script resolves the repository's `figures/` directory, reads the live sheet through `googlesheets4`, and may request authentication. It installs any missing R packages when executed.

The seven current export definitions, projections, palettes, and output names are retained. The existing committed PNGs are preserved byte-for-byte in this redesign. Before regenerating them, review the live data and the script's coordinate-cleaning, target-group filtering, deduplication, and automatic coordinate-swap logic.

## Interactive explorer

`index.html` loads assets from `dashboard/`. The canonical embedded snapshot is `data/literature/eukaryotes.snapshot.json`; provenance is recorded alongside it. This snapshot is extracted unchanged from the previous dashboard, with 319 source rows. Its original capture date is unknown; extraction date is not a claim about data freshness.

The map, annual and cumulative charts, reference-database chart, and mapper chart share the active filters. Search and the paginated source table allow records to be inspected and downloaded. PNG buttons export individual charts. The original static figures appear in their own gallery and do not respond to the interactive filters.

Rows with invalid or missing coordinates remain in the table and chart summaries, but are excluded from the map. This coordinate handling is intentionally more conservative than the R script's automatic swap heuristic; counts are therefore not claimed to reproduce every R filtering decision.

A live refresh or local CSV/TSV import replaces the browser session's dataset. Refresh failure retains the previous data and its source label. Local imports are not uploaded. Only an explicit refresh requests the Google Sheet by default. See [the development guide](../DEVELOPMENT.md) for maintaining the snapshot and previewing the explorer.

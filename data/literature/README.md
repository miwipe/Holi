# Literature snapshot

`eukaryotes.snapshot.json` contains the 319 rows extracted unchanged from the pre-redesign `index.html`. The original capture date was not recorded. `provenance.json` records that limitation, the extraction date, source sheet, and row count.

The data retains pending DOI issues. Corrections belong in the curator's review process before a new snapshot is captured. A site record is a source-table row, not necessarily a unique locality.

Rebuild `dashboard/snapshot.js` with `node scripts/build-dashboard.mjs` after an intentional snapshot update. Live loading and browser imports do not update these files.

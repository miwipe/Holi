> This detailed catalogue guide is retained from the original wiki. Catalogue construction is separate from running Holi. Check index names and shard layout against `bin/holi.sh` before substituting a new build. The existing reference-catalogue dashboard is now at `dashboards/reference-catalog/index.html`.

# Building a Catalog of Reference Genomes for Taxonomic Classification of Metagenomic Data

Constructing a database for taxonomic classification in metagenomics is relatively straightforward with Bowtie2: download FASTA files of reference genomes and run `bowtie2-build` (parallelizable with `--threads <num>`). The challenging part is curating the reference genome collection you choose to index, since this directly influences (i) how many reads are attracted to the database, (ii) what you can classify, and (iii) the certainty of those classifications. It is therefore essential that the database provides strong representation of the organismal groups relevant to your study, while also including additional representative references that can attract reads correctly. In practice, the optimal composition depends on the complexity of your metagenome and the balance of sensitivity versus specificity you seek.

**Resource considerations (Bowtie2 indexing)**

Depending on the size and structure of your input FASTA (for example very long contigs, many contigs, or large combined references), `bowtie2-build` may decide to build a large index and will use an alternative strategy that requires additional memory. A practical rule of thumb is that `bowtie2-build` may use approximately 1.5× the size of the uncompressed FASTA in RAM during indexing. Disk requirements can also grow substantially, since Bowtie2 indices and intermediate files can be several-fold larger than the raw FASTA, depending on your build strategy and file management.

**The real challenge: building a representative reference genome catalog**

The more complex and critical part is curating a representative reference genome catalog, because this directly affects taxonomic sensitivity and specificity. First, verify whether the taxa or communities you care about are actually represented in your reference set. Second, avoid systematic biases that can skew results. For example, a database dominated by arctic taxa can cause competitive mapping to preferentially assign reads to arctic organisms even when the input originates from another biome.

In Holi, we emphasize including as many reference genomes as possible across a broad phylogenetic spectrum. This matters because the pipeline relies on competitive mapping to determine the taxonomic rank of each read. Accurate classification requires sufficient genomic resolution to differentiate closely related taxa. Inadequate representation or diversity reduces classification power and increases the likelihood of misclassification or missed taxa. In addition, and equally important, genome-wide statistics and DNA damage estimates should meet criteria accepted by the ancient DNA community. Always remember that extraordinary claims require extraordinary support.

## What is a reference genome catalog in Holi?

In Holi, a reference genome catalog is a curated, versionable set of genome assemblies used as mapping targets. It includes:

* a reproducible list of genome assemblies (by accession),
* the corresponding genome FASTA sequences (used to build `Bowtie2` indices), and
* metadata and taxonomy for transparent inclusion criteria and downstream summarization.

The Holi pipeline is configured to use the NCBI ecosystem, hence the genome repository and the NCBI taxonomy (TaxIDs and the NCBI taxonomic tree). This provides stable identifiers and enables consistent aggregation and reporting across taxonomic ranks, while it comes with a few challenges. It is possible to add a custom set of references and/or taxonomy, as long as you follow the NCBI conventions (for example unique accessions in FASTA headers, and corresponding entries in `acc2taxid`, `names.dmp`, and `nodes.dmp`).

## Why this catalog is eukaryote-only

This catalog is intentionally restricted to eukaryotes. In the `holi_prefilter` workflow, reads are first mapped against a microbial database to remove confidently microbial-mapping reads. We then retain reads that are likely eukaryotic, or that do not map to microbes, and pass only those reads to the eukaryote-focused mapping and classification step.

This is a fast and flexible strategy that helps suppress microbial-like reads that can otherwise appear in eukaryotic genome assemblies (for example due to contamination, horizontally transferred regions, or microbial sequence fragments embedded in eukaryotic genome builds). Prefiltering reduces spurious eukaryotic assignments driven by microbial sequence similarity and improves the specificity of downstream competitive mapping.

### Overview of the catalog build workflow

**The reference catalog build consists of three parts:**

**Metadata acquisition (NCBI Datasets CLI)**
Use `ncbi-datasets-cli` to query NCBI for the latest eukaryotic assemblies and download genome-associated statistics and identifiers. This metadata snapshot becomes the basis for both exploration and the final selection of assemblies to download.

**Exploration and filtering (R, data-driven)**
Use the accompanying R workflow to explore key assembly metrics (for example contig counts, N50, coverage, and gap-related measures) and to decide filtering thresholds.
Important: filters should be set using a data-driven approach. Make informed decisions for the specific NCBI snapshot and scientific objective you are targeting, and avoid applying previously suggested thresholds by default unless you have validated that they remain appropriate.

**Download and packaging for indexing**
From the filtered list of NCBI accessions, download and merge genomes into manageable FASTA bundles for Bowtie2 indexing. As a practical rule of thumb, the Bowtie2 index build typically requires approximately 1.5× RAM relative to the uncompressed combined FASTA being indexed, and indices can require multiple-fold additional disk space. Once the database is created, update `config.yml` to point to the new database paths.

## Step 1: Download NCBI assembly metadata and build taxonomic paths

This step creates the metadata tables used for downstream filtering (Step 2) and, later, for selecting and downloading assemblies (Step 3). We query **NCBI Datasets** for the latest *eukaryotic* genome assemblies (excluding atypical assemblies), then produce:

* a compressed JSONL metadata snapshot (`wgs_eukaryota.metadata_filter.jsonl.gz`)
* a comprehensive TSV with many genome and assembly statistics (`wgs_eukaryota.metadata_filter.full.tsv.gz`)
* a smaller “short” TSV with selected columns that are convenient for filtering (`wgs_eukaryota.metadata_filter.tsv.gz`)
* a `taxon_paths.tsv` file mapping `tax_id` to a semicolon-delimited NCBI lineage path

### 1.1 Install NCBI Datasets CLI

Install NCBI Datasets command-line tools (`datasets` + `dataformat`) into the current conda/mamba env.

```
mamba install -c conda-forge -c bioconda ncbi-datasets-cli
```

### 1.2 Query NCBI for eukaryotic genome assemblies (metadata snapshot)

Query NCBI for the latest genome assemblies under the taxon `eukaryota`, exclude atypical assemblies, output as JSON Lines, then gzip-compress to a metadata file.

```
datasets summary genome taxon --exclude-atypical --assembly-version latest --as-json-lines eukaryota | gzip > wgs_eukaryota.metadata_filter.jsonl.gz
```

### 1.3 Convert JSONL to a full TSV of genome and assembly statistics

Decompress the JSONL, remove the `.organism.species` field (`jq del`), convert JSON to a TSV with a large set of selected fields (`dataformat`), de-duplicate identical lines (`awk`), normalize quoting (replace `"` with `'`) and fix a specific metaSPAdes formatting issue, then gzip-compress the final full TSV.

```
zcat wgs_eukaryota.metadata_filter.jsonl.gz
| jq -c 'del(.organism.species)'
| dataformat tsv genome --fields accession,annotinfo-busco-complete,annotinfo-busco-duplicated,annotinfo-busco-fragmented,annotinfo-busco-lineage,annotinfo-busco-missing,annotinfo-busco-singlecopy,annotinfo-busco-totalcount,annotinfo-featcount-gene-non-coding,annotinfo-featcount-gene-other,annotinfo-featcount-gene-protein-coding,annotinfo-featcount-gene-pseudogene,annotinfo-featcount-gene-total,annotinfo-method,annotinfo-name,annotinfo-pipeline,annotinfo-provider,annotinfo-release-date,annotinfo-report-url,annotinfo-software-version,annotinfo-status,assminfo-assembly-method,assminfo-atypicalis-atypical,assminfo-atypicalwarnings,assminfo-bioproject,assminfo-biosample-accession,assminfo-biosample-age,assminfo-biosample-bioproject-parent-accession,assminfo-biosample-bioproject-parent-accessions,assminfo-biosample-bioproject-title,assminfo-biosample-breed,assminfo-biosample-collected-by,assminfo-biosample-collection-date,assminfo-biosample-cultivar,assminfo-biosample-description-comment,assminfo-biosample-description-organism-common-name,assminfo-biosample-description-organism-infraspecific-breed,assminfo-biosample-description-organism-infraspecific-cultivar,assminfo-biosample-description-organism-infraspecific-ecotype,assminfo-biosample-description-organism-infraspecific-isolate,assminfo-biosample-description-organism-infraspecific-sex,assminfo-biosample-description-organism-infraspecific-strain,assminfo-biosample-description-organism-name,assminfo-biosample-description-organism-pangolin,assminfo-biosample-description-organism-tax-id,assminfo-biosample-description-title,assminfo-biosample-development-stage,assminfo-biosample-ecotype,assminfo-biosample-geo-loc-name,assminfo-biosample-host,assminfo-biosample-host-disease,assminfo-biosample-ifsac-category,assminfo-biosample-isolate,assminfo-biosample-isolate-name-alias,assminfo-biosample-isolation-source,assminfo-biosample-last-updated,assminfo-biosample-lat-lon,assminfo-biosample-models,assminfo-biosample-owner-contact-lab,assminfo-biosample-owner-name,assminfo-biosample-package,assminfo-biosample-project-name,assminfo-biosample-publication-date,assminfo-biosample-sample-name,assminfo-biosample-serotype,assminfo-biosample-serovar,assminfo-biosample-sex,assminfo-biosample-source-type,assminfo-biosample-status-status,assminfo-biosample-status-when,assminfo-biosample-strain,assminfo-biosample-sub-species,assminfo-biosample-submission-date,assminfo-biosample-tissue,assminfo-blast-url,assminfo-description,assminfo-grouping-method,assminfo-level,assminfo-long-name,assminfo-name,assminfo-notes,assminfo-paired-assm-accession,assminfo-paired-assm-changed,assminfo-paired-assm-manual-diff,assminfo-paired-assm-name,assminfo-paired-assm-only-genbank,assminfo-paired-assm-only-refseq,assminfo-paired-assm-status,assminfo-refseq-category,assminfo-release-date,assminfo-sequencing-tech,assminfo-status,assminfo-submitter,assminfo-suppression-reason,assminfo-synonym,assminfo-type,assmstats-contig-l50,assmstats-contig-n50,assmstats-gaps-between-scaffolds-count,assmstats-gc-count,assmstats-gc-percent,assmstats-genome-coverage,assmstats-number-of-component-sequences,assmstats-number-of-contigs,assmstats-number-of-organelles,assmstats-number-of-scaffolds,assmstats-scaffold-l50,assmstats-scaffold-n50,assmstats-total-number-of-chromosomes,assmstats-total-sequence-len,assmstats-total-ungapped-len,checkm-completeness,checkm-completeness-percentile,checkm-contamination,checkm-marker-set,checkm-marker-set-rank,checkm-species-tax-id,checkm-version,current-accession,organism-common-name,organism-infraspecific-breed,organism-infraspecific-cultivar,organism-infraspecific-ecotype,organism-infraspecific-isolate,organism-infraspecific-sex,organism-infraspecific-strain,organism-name,organism-pangolin,organism-tax-id,source_database,type_material-display_text,type_material-label,wgs-contigs-url,wgs-project-accession,wgs-url
| awk '!x[$0]++'
| sed "s/"/'/g"
| sed 's/metaSPAdes\tv/metaSPAdes v/'
| gzip > wgs_eukaryota.metadata_filter.full.tsv.gz
```


### 1.4 Create a smaller metadata TSV for interactive filtering

From the “full” TSV, select a smaller subset of columns by name (fuzzy matching), treat input as tab-separated, and write a compressed TSV output.

```
csvtk --tabs cut --fuzzy-fields -f "Assembly Accession","Assembly BioProject Accession","Assembly BioSample Description Organism Name","Organism Taxonomic ID","Source Database","Assembly BioSample Last updated","Assembly BioSample Owner Name","Assembly BioSample Geographic location ","Assembly BioSample Isolation source ","Assembly Level","Assembly Stats *" --out-file wgs_eukaryota.metadata_filter.tsv.gz wgs_eukaryota.metadata_filter.full.tsv.gz
```

### 1.5 Download NCBI taxonomy dump and configure `taxonkit`

Install `taxonkit`, download the NCBI taxonomy dump, unpack it, and configure the local `taxonkit` database.

```
conda install -c bioconda taxonkit
wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

tar -zxvf taxdump.tar.gz

mkdir -p $HOME/.taxonkit
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit
```

### 1.6 Build taxonomic lineage paths for all assemblies

Extract TaxIDs from the full TSV (column 45 in the current export), drop the header, then compute semicolon-delimited taxonomic paths.

```
zcat wgs_eukaryota.metadata_filter.full.tsv.gz | cut -f45 | tail -n+2 | taxonkit reformat2 - > taxon_paths.tsv
```

**Outputs of Step 1**

* `wgs_eukaryota.metadata_filter.jsonl.gz` (raw NCBI datasets JSONL snapshot)
* `wgs_eukaryota.metadata_filter.full.tsv.gz` (full TSV export with extensive assembly stats)
* `wgs_eukaryota.metadata_filter.tsv.gz` (short TSV used for Step 2 filtering)
* `taxon_paths.tsv` (TaxID to lineage path mapping used to annotate/filter assemblies)

**Notes / assumptions**

* The TaxID column index (`cut -f45`) depends on the exact `dataformat --fields` selection and ordering. If you change the fields list, verify the TaxID column position before running the `cut`.
* `--exclude-atypical` removes assemblies flagged as atypical by NCBI. If you need atypical assemblies for a specific use case, remove this flag and document the rationale.


## Step 2: Explore metadata and define data-driven filters (R)

This step is where you decide *which* assemblies to include in your reference catalog. The goal is to explore assembly-level statistics (for example contig counts, N50, coverage, and gap structure), compare patterns across sources (for example GenBank vs RefSeq), and define filtering thresholds that are appropriate for **your** scientific objective and the specific NCBI snapshot you downloaded in Step 1.

**Important:** Filters should be set using a data-driven approach. Avoid applying any suggested thresholds by default unless you have validated that they remain appropriate for the current catalog snapshot and your target organismal space.

### Recommended ways to run the R workflow

There are two common ways to use the R workflow below, and they serve different purposes:

1. **Interactive exploration (recommended for developing filters)**
   Open the script in RStudio (or another R GUI), run it section-by-section, inspect plots, and iteratively adjust thresholds. This is the best approach when you are creating a new catalog version or changing filtering logic, because you will want to visually confirm what gets filtered and why.

2. **Batch execution (recommended for reproducible re-runs)**
   Once you have agreed on filtering logic and thresholds, run the same script non-interactively (for example on a server) using `Rscript`. In this mode, the script should write plots and tables to disk (for example `figures/` and `tables/`) so that the run is reproducible and auditable.

A practical workflow is therefore:
- first run the script interactively to *choose* thresholds and filters,
- then lock in those thresholds (for example in a config section at the top of the script),
- and finally run it as `Rscript` to generate the final “kept” and “removed” assembly lists used for Step 3.

### Batch run example (optional)

If you want to run non-interactively, consider rendering as an HTML report (recommended) or running as a plain script:

Render to HTML (produces a self-contained report with plots):

```
Rscript -e "rmarkdown::render('Final_script_for_filteringGenbankRefSeq_genomes.Rmd')"
```

Run as a plain R script (recommended only if the script explicitly saves outputs):
```
Rscript Final_script_for_filteringGenbankRefSeq_genomes.R
```

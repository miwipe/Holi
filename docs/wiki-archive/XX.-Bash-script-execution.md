<html><head></head><body><h1>Holi tutorial: from clone to running the pipeline</h1>
<p>This guide explains how to:</p>
<ol>
<li>Clone the Holi repo into your working folder</li>
<li>Create the holi conda environment from <code>holi.yml</code></li>
<li>Create <code>sample.list</code> with your sample names</li>
<li>Run long jobs safely (screen, tmux, nohup, Slurm)</li>
<li>Set up <code>config.yml</code> with correct paths to folders, databases, and taxonomies</li>
<li>Test that all tools are installed and print versions</li>
<li>Run <code>holi_prefilter.sh</code> with different options and understand why you would use them</li>
</ol>
<hr>
<h2>What Holi does (high level)</h2>
<p>The pipeline runs in 8 stages:</p>

Stage | Name | Description | Key output
-- | -- | -- | --
1 | Preprocessing | fastp merge + trim → vsearch duplicate removal → sga low-complexity filter | <sample>.ppm.vs.d4.fq.gz
2 | GTDB mapping | bowtie2 to 7 GTDB chunks → sort → merge | MICROB_OUT/<sample>.gtdb.merged.bam
3 | Microbial split | getRTax or metaDMG LCA → bacterial + eukaryotic FASTQs | EUK_OUT/<sample>.euk.fastq.gz
4 | Euk mapping | bowtie2 to euk DB parts + mito + phyNor + core_nt + plastid | EUK_OUT/<sample>.pla.bam (last DB mapped)
5 | Compress/merge | compress mapping BAMs → sort → merge | EUK_OUT/<sample>.comp.bam
6 | BAM filtering | Unicorn: alnfilt → refstats → bamstats | EUK_OUT/<sample>.comp.filtered.bam
7 | metaDMG | LCA → dfit → aggregate | EUK_OUT/<sample>.sort.comp.filtered.agg.stat.gz
8 | Unicorn taxstats | Per-taxID stats at family, genus, and species rank | EUK_OUT/<sample>.comp.filtered.species.taxstats


<hr>
<h3>7.7 The Unicorn filtering workflow (stage 6)</h3>
<p>Stage 6 runs a three-step Unicorn pipeline:</p>
<ol>
<li><strong>alnfilt</strong> — alignment-level filtering (<code>--minani 93–100</code>, <code>--mode ALLTOP</code>)</li>
<li><strong>refstats</strong> — reference statistics filtering (reads from <code>alnfilt</code> output)</li>
<li><strong>bamstats</strong> — final BAM statistics and output of <code>comp.filtered.bam</code></li>
</ol>
<p>This is the only filtering method available. The output of stage 6 is <code>EUK_OUT/&lt;sample&gt;.comp.filtered.bam</code>.</p>
<hr>
<h3>7.8 Common use cases</h3>
<p>Run everything from scratch:</p>
<pre><code class="language-bash">./holi.sh config.yml sample.list
</code></pre>
<p>Resume after a crash (auto-detection handles this automatically):</p>
<pre><code class="language-bash">./holi.sh config.yml sample.list
</code></pre>
<p>Rerun only metaDMG and taxstats on existing filtered BAMs:</p>
<pre><code class="language-bash">./holi.sh config.yml sample.list --from-stage 7
</code></pre>
<p>Run only preprocessing and GTDB mapping, then stop:</p>
<pre><code class="language-bash">./holi.sh config.yml sample.list --to-stage 2
</code></pre>
<p>Run single-end data, keep all intermediates, skip metaDMG:</p>
<pre><code class="language-bash">./holi.sh config.yml sample.list --single-end --skip-preprocessing-cleanup --skip-gtdb-cleanup --skip-euk-cleanup --to-stage 6
</code></pre>
<hr>
<h2>8) Common pitfalls and fixes</h2>
<h3>Wrong environment or PATH contamination</h3>
<pre><code class="language-bash">echo "$CONDA_DEFAULT_ENV"
which bowtie2 samtools
</code></pre>
<h3>Config keys not found due to indentation</h3>
<p>Keys must start at column 1. <code>  THREADS: 12</code> will not be parsed — use <code>THREADS: 12</code>.</p>
<h3>Options placement</h3>
<p>Options must come after the two required positional arguments:</p>
<pre><code class="language-bash">./holi.sh config.yml sample.list --skip-preprocessing   # correct
./holi.sh --skip-preprocessing config.yml sample.list   # wrong
</code></pre>
<h3>Oversubscribing CPUs</h3>
<p>Total CPU use ≈ <code>THREADSP * THREADS</code>. If jobs are slow or the node is overloaded, reduce <code>THREADSP</code> first.</p>
<h3>FASTQ discovery problems</h3>
<pre><code class="language-bash">sample=$(head -n 1 sample.list)
ls "${INPATH}/${sample}"*R1*.fastq.gz
ls "${INPATH}/${sample}"*R2*.fastq.gz
</code></pre>
<hr>
<h2>9) Quick start command sequence</h2>
<pre><code class="language-bash">git clone https://github.com/miwipe/Holi.git
cd Holi

conda env create -f holi.yml
conda activate holi

ls /path/to/fastq/*R1*.fastq.gz | awk -F'/' '{print $NF}' | cut -f1 -d_ | sort -u &gt; sample.list

vim config.yml

./test_tools_safe.sh

./holi_prefilter.sh config.yml sample.list
</code></pre></body></html>
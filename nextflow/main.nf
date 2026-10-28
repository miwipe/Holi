nextflow.enable.dsl=2

// Shell quoting for staged filenames, including spaces and apostrophes.
def quote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

process VALIDATE {
    tag 'manifests'
    cpus 1
    memory '1 GB'
    publishDir "${params.outdir}/validation", mode: 'copy', overwrite: true
    input:
    path samples, stageAs: 'samples.csv'
    path settings, stageAs: 'settings.json'
    path references, stageAs: 'reference_data'
    path runner, stageAs: 'runner/*'
    output:
    path 'validated.csv', emit: samples
    path 'validation.json', emit: summary
    script:
    """
    python3 runner/prepare.py --input samples.csv --settings settings.json --references reference_data --check-tools
    """
    stub:
    """
    python3 runner/prepare.py --input samples.csv --settings settings.json --references reference_data
    """
}

process PREPROCESS {
    tag sample
    publishDir "${params.outdir}/stage1", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path(reads, stageAs: 'reads/*')
    path settings, stageAs: 'settings.json'
    path references, stageAs: 'reference_data'
    path runner, stageAs: 'runner/*'
    output:
    tuple val(sample), path("${sample}.stage1")
    script:
    """
    python3 runner/stage.py --stage 1 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --reads ${(reads instanceof List ? reads : [reads]).collect { quote(it) }.join(' ')}
    """
    stub:
    """
    python3 runner/stage.py --stage 1 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --reads ${(reads instanceof List ? reads : [reads]).collect { quote(it) }.join(' ')} --stub
    """
}

process MAP_MICROBIAL {
    tag sample
    publishDir "${params.outdir}/stage2", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path(upstream, stageAs: 'upstream')
    path settings, stageAs: 'settings.json'
    path references, stageAs: 'reference_data'
    path runner, stageAs: 'runner/*'
    output:
    tuple val(sample), path(upstream), path("${sample}.stage2")
    script:
    """
    python3 runner/stage.py --stage 2 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream
    """
    stub:
    """
    python3 runner/stage.py --stage 2 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream --stub
    """
}

process SPLIT_READS {
    tag sample
    publishDir "${params.outdir}/stage3", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path(reads_dir, stageAs: 'preprocessed'), path(upstream, stageAs: 'upstream')
    path settings, stageAs: 'settings.json'
    path references, stageAs: 'reference_data'
    path runner, stageAs: 'runner/*'
    output:
    tuple val(sample), path("${sample}.stage3")
    script:
    """
    python3 runner/stage.py --stage 3 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream --reads ${quote(reads_dir.resolve(sample + '.ppm.vs.d4.fq.gz'))}
    """
    stub:
    """
    python3 runner/stage.py --stage 3 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream --reads ${quote(reads_dir.resolve(sample + '.ppm.vs.d4.fq.gz'))} --stub
    """
}

process MAP_EUKARYOTIC {
    tag sample
    publishDir "${params.outdir}/stage4", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path(upstream, stageAs: 'upstream')
    path settings, stageAs: 'settings.json'
    path references, stageAs: 'reference_data'
    path runner, stageAs: 'runner/*'
    output:
    tuple val(sample), path("${sample}.stage4")
    script:
    """
    python3 runner/stage.py --stage 4 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream
    """
    stub:
    """
    python3 runner/stage.py --stage 4 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream --stub
    """
}

process MERGE_ALIGNMENTS {
    tag sample
    publishDir "${params.outdir}/stage5", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path(upstream, stageAs: 'upstream')
    path settings, stageAs: 'settings.json'
    path references, stageAs: 'reference_data'
    path runner, stageAs: 'runner/*'
    output:
    tuple val(sample), path("${sample}.stage5")
    script:
    """
    python3 runner/stage.py --stage 5 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream
    """
    stub:
    """
    python3 runner/stage.py --stage 5 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream --stub
    """
}

process FILTER_ALIGNMENTS {
    tag sample
    publishDir "${params.outdir}/stage6", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path(upstream, stageAs: 'upstream')
    path settings, stageAs: 'settings.json'
    path references, stageAs: 'reference_data'
    path runner, stageAs: 'runner/*'
    output:
    tuple val(sample), path("${sample}.stage6")
    script:
    """
    python3 runner/stage.py --stage 6 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream
    """
    stub:
    """
    python3 runner/stage.py --stage 6 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream --stub
    """
}

process ESTIMATE_DAMAGE {
    tag sample
    publishDir "${params.outdir}/stage7", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path(upstream, stageAs: 'upstream')
    path settings, stageAs: 'settings.json'
    path references, stageAs: 'reference_data'
    path runner, stageAs: 'runner/*'
    output:
    tuple val(sample), path("${sample}.stage7")
    script:
    """
    python3 runner/stage.py --stage 7 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream
    """
    stub:
    """
    python3 runner/stage.py --stage 7 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream --stub
    """
}

process TAXONOMY_STATS {
    tag sample
    publishDir "${params.outdir}/stage8", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path(upstream, stageAs: 'upstream')
    path settings, stageAs: 'settings.json'
    path references, stageAs: 'reference_data'
    path runner, stageAs: 'runner/*'
    output:
    tuple val(sample), path("${sample}.stage8")
    script:
    """
    python3 runner/stage.py --stage 8 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream
    """
    stub:
    """
    python3 runner/stage.py --stage 8 --sample ${quote(sample)} --settings settings.json --references reference_data --cpus ${task.cpus} --memory-mb ${task.memory.toMega()} --upstream upstream --stub
    """
}

workflow {
    if (params.help) {
        log.info 'Holi: nextflow run nextflow/main.nf --input samples.csv --settings settings.json --references /absolute/reference_bundle --outdir results/run01 [-resume]'
    } else {
        if (!params.input || !params.settings || !params.references) {
            error 'Required: --input samples.csv --settings settings.json --references reference_directory (see nextflow/README.md)'
        }
        samples_file = file(params.input, checkIfExists: true)
        settings = file(params.settings, checkIfExists: true)
        references = file(params.references, checkIfExists: true)
        runner = [file("${projectDir}/runner/prepare.py"), file("${projectDir}/runner/stage.py"), file("${projectDir}/runner/defaults.json")]
        VALIDATE(samples_file, settings, references, runner)
        samples = VALIDATE.out.samples.splitCsv(header: true).map { row ->
            def reads = [file(row.fastq_1, checkIfExists: true)]
            if (row.fastq_2) { reads.add(file(row.fastq_2, checkIfExists: true)) }
            tuple(row.sample, reads)
        }
        PREPROCESS(samples, settings, references, runner)
        MAP_MICROBIAL(PREPROCESS.out, settings, references, runner)
        SPLIT_READS(MAP_MICROBIAL.out, settings, references, runner)
        MAP_EUKARYOTIC(SPLIT_READS.out, settings, references, runner)
        MERGE_ALIGNMENTS(MAP_EUKARYOTIC.out, settings, references, runner)
        FILTER_ALIGNMENTS(MERGE_ALIGNMENTS.out, settings, references, runner)
        ESTIMATE_DAMAGE(FILTER_ALIGNMENTS.out, settings, references, runner)
        TAXONOMY_STATS(FILTER_ALIGNMENTS.out, settings, references, runner)
    }
}

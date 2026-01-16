#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    Bovine IgG Repertoire Analysis Pipeline
========================================================================================
    Analyzes antibody repertoire diversity using the Immcantation framework

    Tools used:
    - IgBLAST: V(D)J gene assignment
    - Change-O: Clonal assignment and AIRR formatting
    - Alakazam: Diversity analysis and visualization

    Input: FASTA files from bovine-igg-pipeline (unique sequences)
    Output: Clonotype diversity metrics, CDR3 analysis, gene usage plots
========================================================================================
*/

// --------------------------------------------------------------- //
// WORKFLOW
// --------------------------------------------------------------- //

workflow {

    // Input channel from FASTA files
    ch_fasta = Channel
        .fromPath( params.fasta_input )
        .map { fasta ->
            def sample_id = fasta.baseName.replaceAll(/_unique$/, '').replaceAll(/_nogroup$/, '')
            tuple( sample_id, fasta )
        }

    // Download/prepare germline database if not provided
    if ( params.germline_db ) {
        ch_germline = Channel.fromPath( params.germline_db )
    } else {
        FETCH_GERMLINES()
        ch_germline = FETCH_GERMLINES.out
    }

    // Build IgBLAST database
    BUILD_IGBLAST_DB( ch_germline.collect() )

    // Run IgBLAST annotation
    IGBLAST_ANNOTATION(
        ch_fasta,
        BUILD_IGBLAST_DB.out.db,
        BUILD_IGBLAST_DB.out.aux
    )

    // Convert IgBLAST output to AIRR format
    MAKEDB(
        IGBLAST_ANNOTATION.out,
        ch_germline.collect()
    )

    // Filter for productive sequences
    FILTER_PRODUCTIVE(
        MAKEDB.out
    )

    // Define clones
    DEFINE_CLONES(
        FILTER_PRODUCTIVE.out
    )

    // Collect all samples for combined analysis
    ch_all_clones = DEFINE_CLONES.out.collect()

    // Run diversity analysis
    DIVERSITY_ANALYSIS(
        ch_all_clones
    )

    // Generate summary report
    GENERATE_REPORT(
        DIVERSITY_ANALYSIS.out.plots,
        DIVERSITY_ANALYSIS.out.stats
    )

}

// --------------------------------------------------------------- //
// PARAMETERS
// --------------------------------------------------------------- //

params.fasta_input = "${launchDir}/*_unique.fasta"
params.germline_db = null  // Path to pre-built germline FASTA files
params.species = "bovine"
params.chain = "IG"  // IG for immunoglobulin
params.results = "${launchDir}/results"

// Clone definition parameters
params.clone_threshold = 0.15  // Junction hamming distance threshold
params.clone_method = "hierarchical"  // clustering method

// Analysis parameters
params.min_seq_count = 1  // Minimum sequences per clone
params.subsample = null  // Subsample size for diversity (null = no subsampling)


// --------------------------------------------------------------- //
// PROCESSES
// --------------------------------------------------------------- //

process FETCH_GERMLINES {

    publishDir "${params.results}/germlines", mode: 'copy'

    output:
    path "*.fasta"

    script:
    """
    # Fetch bovine germlines from IMGT using Immcantation's fetch script
    fetch_imgtdb.sh -o . -s bovine

    # Rename files to standard format if needed
    for f in *.fasta; do
        if [[ ! -f "\$f" ]]; then
            echo "Warning: No germline files downloaded. Please check IMGT availability."
            exit 1
        fi
    done
    """
}

process BUILD_IGBLAST_DB {

    publishDir "${params.results}/igblast_db", mode: 'copy'

    input:
    path germlines

    output:
    path "database/*", emit: db
    path "internal_data/*", emit: aux

    script:
    """
    mkdir -p database internal_data

    # Separate V, D, J genes
    cat *IGHV*.fasta *IGKV*.fasta *IGLV*.fasta 2>/dev/null > database/bovine_V.fasta || true
    cat *IGHD*.fasta 2>/dev/null > database/bovine_D.fasta || true
    cat *IGHJ*.fasta *IGKJ*.fasta *IGLJ*.fasta 2>/dev/null > database/bovine_J.fasta || true

    # Build BLAST databases
    cd database
    for f in bovine_*.fasta; do
        if [[ -s "\$f" ]]; then
            makeblastdb -parse_seqids -dbtype nucl -in "\$f"
        fi
    done
    cd ..

    # Create IgBLAST auxiliary data
    # Copy internal data from IgBLAST installation
    cp -r \${IGDATA}/internal_data/* internal_data/ 2>/dev/null || true

    # Create organism-specific files if not present
    if [[ ! -f internal_data/bovine/bovine_gl.aux ]]; then
        mkdir -p internal_data/bovine
        # Create basic aux file for bovine
        echo -e "# Bovine germline auxiliary data" > internal_data/bovine/bovine_gl.aux
    fi
    """
}

process IGBLAST_ANNOTATION {

    tag "${sample_id}"
    publishDir "${params.results}/igblast", mode: 'copy'

    cpus 4

    input:
    tuple val(sample_id), path(fasta)
    path db
    path aux

    output:
    tuple val(sample_id), path("${sample_id}_igblast.fmt7"), path(fasta)

    script:
    """
    # Run IgBLAST
    AssignGenes.py igblast \
        -s ${fasta} \
        -b \${IGDATA} \
        --organism bovine \
        --loci ig \
        --format blast \
        -o ${sample_id}_igblast.fmt7 \
        --nproc ${task.cpus}
    """
}

process MAKEDB {

    tag "${sample_id}"
    publishDir "${params.results}/airr", mode: 'copy'

    input:
    tuple val(sample_id), path(igblast_out), path(fasta)
    path germlines

    output:
    tuple val(sample_id), path("${sample_id}_db-pass.tsv")

    script:
    """
    # Combine germline references
    cat *IGHV*.fasta *IGKV*.fasta *IGLV*.fasta 2>/dev/null > combined_V.fasta || true
    cat *IGHD*.fasta 2>/dev/null > combined_D.fasta || true
    cat *IGHJ*.fasta *IGKJ*.fasta *IGLJ*.fasta 2>/dev/null > combined_J.fasta || true

    # Convert IgBLAST output to AIRR format
    MakeDb.py igblast \
        -i ${igblast_out} \
        -s ${fasta} \
        -r combined_V.fasta combined_D.fasta combined_J.fasta \
        --extended \
        --format airr \
        -o ${sample_id}_db.tsv

    # Rename output
    mv ${sample_id}_db_db-pass.tsv ${sample_id}_db-pass.tsv 2>/dev/null || mv ${sample_id}_db-pass.tsv ${sample_id}_db-pass.tsv
    """
}

process FILTER_PRODUCTIVE {

    tag "${sample_id}"
    publishDir "${params.results}/filtered", mode: 'copy'

    input:
    tuple val(sample_id), path(airr_tsv)

    output:
    tuple val(sample_id), path("${sample_id}_productive.tsv")

    script:
    """
    # Filter for productive sequences (in-frame, no stop codons)
    ParseDb.py select \
        -d ${airr_tsv} \
        -f productive \
        -u T TRUE True true \
        -o ${sample_id}_productive.tsv
    """
}

process DEFINE_CLONES {

    tag "${sample_id}"
    publishDir "${params.results}/clones", mode: 'copy'

    input:
    tuple val(sample_id), path(airr_tsv)

    output:
    path "${sample_id}_clones.tsv"

    script:
    """
    # Define clones based on V gene, J gene, and junction similarity
    DefineClones.py -d ${airr_tsv} \
        --act set \
        --model ham \
        --norm len \
        --dist ${params.clone_threshold} \
        --format airr \
        -o ${sample_id}_clones.tsv
    """
}

process DIVERSITY_ANALYSIS {

    publishDir "${params.results}/diversity", mode: 'copy'

    input:
    path clone_files

    output:
    path "plots/*", emit: plots
    path "stats/*", emit: stats

    script:
    """
    mkdir -p plots stats

    # Run R script for diversity analysis
    Rscript - <<'RSCRIPT'

    library(alakazam)
    library(shazam)
    library(ggplot2)
    library(dplyr)
    library(airr)

    # Read all clone files
    files <- list.files(pattern = "_clones.tsv\$", full.names = TRUE)

    if (length(files) == 0) {
        stop("No clone files found")
    }

    # Combine all data
    db_list <- lapply(files, function(f) {
        df <- read_rearrangement(f)
        df\$sample_id <- gsub("_clones.tsv", "", basename(f))
        return(df)
    })
    db <- bind_rows(db_list)

    # Basic stats
    stats <- db %>%
        group_by(sample_id) %>%
        summarise(
            total_sequences = n(),
            unique_clones = n_distinct(clone_id, na.rm = TRUE),
            productive = sum(productive == TRUE | productive == "T", na.rm = TRUE),
            mean_cdr3_length = mean(nchar(as.character(junction_aa)), na.rm = TRUE),
            median_cdr3_length = median(nchar(as.character(junction_aa)), na.rm = TRUE)
        )
    write.csv(stats, "stats/basic_stats.csv", row.names = FALSE)

    # CDR3 length distribution
    db\$cdr3_length <- nchar(as.character(db\$junction_aa))

    p1 <- ggplot(db, aes(x = cdr3_length, fill = sample_id)) +
        geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
        labs(title = "CDR3 Length Distribution",
             x = "CDR3 Length (amino acids)",
             y = "Count") +
        theme_minimal() +
        theme(legend.position = "bottom")
    ggsave("plots/cdr3_length_distribution.pdf", p1, width = 10, height = 6)
    ggsave("plots/cdr3_length_distribution.png", p1, width = 10, height = 6, dpi = 150)

    # V gene usage
    v_usage <- db %>%
        filter(!is.na(v_call)) %>%
        mutate(v_gene = gsub("\\\\*.*", "", v_call)) %>%
        group_by(sample_id, v_gene) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(sample_id) %>%
        mutate(freq = count / sum(count))

    write.csv(v_usage, "stats/v_gene_usage.csv", row.names = FALSE)

    p2 <- ggplot(v_usage, aes(x = reorder(v_gene, -freq), y = freq, fill = sample_id)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "V Gene Usage",
             x = "V Gene",
             y = "Frequency") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
              legend.position = "bottom")
    ggsave("plots/v_gene_usage.pdf", p2, width = 14, height = 6)
    ggsave("plots/v_gene_usage.png", p2, width = 14, height = 6, dpi = 150)

    # J gene usage
    j_usage <- db %>%
        filter(!is.na(j_call)) %>%
        mutate(j_gene = gsub("\\\\*.*", "", j_call)) %>%
        group_by(sample_id, j_gene) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(sample_id) %>%
        mutate(freq = count / sum(count))

    write.csv(j_usage, "stats/j_gene_usage.csv", row.names = FALSE)

    p3 <- ggplot(j_usage, aes(x = reorder(j_gene, -freq), y = freq, fill = sample_id)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "J Gene Usage",
             x = "J Gene",
             y = "Frequency") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom")
    ggsave("plots/j_gene_usage.pdf", p3, width = 10, height = 6)
    ggsave("plots/j_gene_usage.png", p3, width = 10, height = 6, dpi = 150)

    # Clone size distribution
    clone_sizes <- db %>%
        filter(!is.na(clone_id)) %>%
        group_by(sample_id, clone_id) %>%
        summarise(clone_size = n(), .groups = "drop")

    write.csv(clone_sizes, "stats/clone_sizes.csv", row.names = FALSE)

    p4 <- ggplot(clone_sizes, aes(x = clone_size, fill = sample_id)) +
        geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
        scale_x_log10() +
        labs(title = "Clone Size Distribution",
             x = "Clone Size (log10)",
             y = "Count") +
        theme_minimal() +
        theme(legend.position = "bottom")
    ggsave("plots/clone_size_distribution.pdf", p4, width = 10, height = 6)
    ggsave("plots/clone_size_distribution.png", p4, width = 10, height = 6, dpi = 150)

    # Diversity curves using Alakazam
    if (nrow(db) > 0 && any(!is.na(db\$clone_id))) {

        # Calculate diversity
        div_curve <- tryCatch({
            alphaDiversity(db, group = "sample_id", clone = "clone_id",
                          min_q = 0, max_q = 4, step_q = 0.1,
                          ci = 0.95, nboot = 100)
        }, error = function(e) {
            message("Could not calculate diversity curve: ", e\$message)
            NULL
        })

        if (!is.null(div_curve)) {
            p5 <- plot(div_curve, legend_title = "Sample") +
                labs(title = "Repertoire Diversity (Hill Numbers)") +
                theme_minimal()
            ggsave("plots/diversity_curve.pdf", p5, width = 10, height = 6)
            ggsave("plots/diversity_curve.png", p5, width = 10, height = 6, dpi = 150)

            # Save diversity values
            write.csv(div_curve@diversity, "stats/diversity_values.csv", row.names = FALSE)
        }

        # Rarefaction curve
        rarefaction <- tryCatch({
            estimateAbundance(db, group = "sample_id", clone = "clone_id",
                            ci = 0.95, nboot = 100)
        }, error = function(e) {
            message("Could not calculate rarefaction: ", e\$message)
            NULL
        })

        if (!is.null(rarefaction)) {
            p6 <- plot(rarefaction, legend_title = "Sample") +
                labs(title = "Clonal Abundance Rarefaction") +
                theme_minimal()
            ggsave("plots/rarefaction_curve.pdf", p6, width = 10, height = 6)
            ggsave("plots/rarefaction_curve.png", p6, width = 10, height = 6, dpi = 150)
        }
    }

    # Summary diversity metrics
    diversity_summary <- db %>%
        filter(!is.na(clone_id)) %>%
        group_by(sample_id) %>%
        summarise(
            total_sequences = n(),
            unique_clones = n_distinct(clone_id),
            simpson_index = 1 - sum((table(clone_id)/n())^2),
            shannon_index = -sum((table(clone_id)/n()) * log(table(clone_id)/n())),
            chao1 = n_distinct(clone_id) + (sum(table(clone_id) == 1)^2) / (2 * max(1, sum(table(clone_id) == 2))),
            .groups = "drop"
        )

    write.csv(diversity_summary, "stats/diversity_summary.csv", row.names = FALSE)

    message("Diversity analysis complete!")

    RSCRIPT
    """
}

process GENERATE_REPORT {

    publishDir "${params.results}", mode: 'copy'

    input:
    path plots
    path stats

    output:
    path "repertoire_report.html"

    script:
    """
    # Generate HTML report
    cat > repertoire_report.html << 'HTML'
<!DOCTYPE html>
<html>
<head>
    <title>Bovine IgG Repertoire Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        .stats-table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        .stats-table th, .stats-table td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        .stats-table th { background-color: #3498db; color: white; }
        .stats-table tr:nth-child(even) { background-color: #f2f2f2; }
        .plot-container { margin: 20px 0; text-align: center; }
        .plot-container img { max-width: 100%; border: 1px solid #ddd; }
        .section { margin: 30px 0; }
    </style>
</head>
<body>
    <h1>Bovine IgG Repertoire Analysis Report</h1>
    <p>Generated: $(date)</p>

    <div class="section">
        <h2>Summary Statistics</h2>
        <p>See stats/ directory for detailed CSV files.</p>
    </div>

    <div class="section">
        <h2>CDR3 Length Distribution</h2>
        <div class="plot-container">
            <img src="diversity/plots/cdr3_length_distribution.png" alt="CDR3 Length Distribution">
        </div>
    </div>

    <div class="section">
        <h2>V Gene Usage</h2>
        <div class="plot-container">
            <img src="diversity/plots/v_gene_usage.png" alt="V Gene Usage">
        </div>
    </div>

    <div class="section">
        <h2>J Gene Usage</h2>
        <div class="plot-container">
            <img src="diversity/plots/j_gene_usage.png" alt="J Gene Usage">
        </div>
    </div>

    <div class="section">
        <h2>Clone Size Distribution</h2>
        <div class="plot-container">
            <img src="diversity/plots/clone_size_distribution.png" alt="Clone Size Distribution">
        </div>
    </div>

    <div class="section">
        <h2>Diversity Analysis</h2>
        <div class="plot-container">
            <img src="diversity/plots/diversity_curve.png" alt="Diversity Curve">
        </div>
        <div class="plot-container">
            <img src="diversity/plots/rarefaction_curve.png" alt="Rarefaction Curve">
        </div>
    </div>

    <div class="section">
        <h2>Output Files</h2>
        <ul>
            <li><strong>airr/</strong> - AIRR-formatted sequence annotations</li>
            <li><strong>clones/</strong> - Clone assignments</li>
            <li><strong>diversity/stats/</strong> - CSV files with diversity metrics</li>
            <li><strong>diversity/plots/</strong> - PDF and PNG plots</li>
        </ul>
    </div>

</body>
</html>
HTML
    """
}

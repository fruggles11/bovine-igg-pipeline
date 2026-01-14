#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	// Input channels - heavy and light chain barcodes
	ch_heavy_dir = Channel
		.fromPath( "${params.fastq_dir}/${params.heavy_barcode}/", type: 'dir' )
		.map { dir -> tuple( "heavy", dir ) }

	ch_light_dir = Channel
		.fromPath( "${params.fastq_dir}/${params.light_barcode}/", type: 'dir' )
		.map { dir -> tuple( "light", dir ) }

	ch_input_dirs = ch_heavy_dir.mix( ch_light_dir )

	ch_primers = Channel
		.fromPath( params.primer_table )
		.splitCsv( header: true )
		.map { row -> tuple( row.chain, row.primer_seq, row.differentiating ) }
		.filter { it[2] == "true" }
		.groupTuple( by: 0 )
		.map { chain, primers, diff -> tuple( chain, primers ) }

	// Germline gene files (if available)
	ch_germlines = Channel
		.fromPath( "${params.germline_dir}/*.fasta", checkIfExists: false )
		.collect()
		.ifEmpty( [] )

	// Workflow steps

	// Stage 1: Read Processing
	MERGE_READS(
		ch_input_dirs
	)

	QUALITY_FILTER(
		MERGE_READS.out
			.map { chain, fastq -> tuple( chain, file(fastq), file(fastq).countFastq() ) }
			.filter { it[2] >= params.min_reads }
			.map { chain, fastq, count -> tuple( chain, fastq ) }
	)

	FIND_ADAPTERS(
		QUALITY_FILTER.out
	)

	TRIM_PRIMERS(
		QUALITY_FILTER.out
			.join( ch_primers )
			.join( FIND_ADAPTERS.out )
	)

	// Stage 2: Clustering & Consensus
	CONVERT_TO_FASTA(
		TRIM_PRIMERS.out
	)

	CLUSTER_READS(
		CONVERT_TO_FASTA.out
	)

	// Stage 3: Annotation (conditional on germline files being available)
	if ( !params.skip_annotation ) {
		BUILD_IGBLAST_DB(
			ch_germlines
		)

		ANNOTATE_IGBLAST(
			BUILD_IGBLAST_DB.out,
			CLUSTER_READS.out
		)

		PARSE_ANNOTATIONS(
			ANNOTATE_IGBLAST.out
		)

		// Stage 4: Reporting
		COLLECT_STATS(
			PARSE_ANNOTATIONS.out.collect()
		)
	} else {
		// Skip annotation, just collect consensus stats
		COLLECT_CONSENSUS_STATS(
			CLUSTER_READS.out.collect()
		)
	}

	if ( !params.skip_annotation ) {
		GENERATE_REPORT(
			COLLECT_STATS.out
		)
	} else {
		GENERATE_REPORT(
			COLLECT_CONSENSUS_STATS.out
		)
	}

}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
if ( params.debugmode == true ){
	errorMode = 'terminate'
} else {
	errorMode = 'ignore'
}

params.merged_reads = params.results + "/1_merged_reads"
params.filtered_reads = params.results + "/2_filtered_reads"
params.consensus_seqs = params.results + "/3_consensus_sequences"
params.annotations = params.results + "/4_annotations"
params.reports = params.results + "/5_reports"
// --------------------------------------------------------------- //



// PROCESS SPECIFICATION
// --------------------------------------------------------------- //

process MERGE_READS {

	tag "${chain}"
	publishDir params.merged_reads, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(chain), path(read_dir)

	output:
	tuple val(chain), path("${chain}_chain.fastq.gz")

	script:
	"""
	seqkit scat -j ${task.cpus} -f `realpath ${read_dir}` -o ${chain}_chain.fastq.gz
	"""

}

process QUALITY_FILTER {

	tag "${chain}"
	publishDir params.filtered_reads, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(chain), path(reads)

	output:
	tuple val(chain), path("${chain}_chain_filtered.fastq.gz")

	script:
	"""
	seqkit seq \
	--min-len ${params.min_len} \
	--max-len ${params.max_len} \
	--min-qual ${params.min_qual} \
	--validate-seq \
	--threads ${task.cpus} \
	${reads} \
	-o ${chain}_chain_filtered.fastq.gz
	"""

}

process FIND_ADAPTERS {

	tag "${chain}"

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	tuple val(chain), path(reads)

	output:
	tuple val(chain), path("${chain}_adapters.fasta")

	script:
	"""
	bbmerge.sh in=`realpath ${reads}` outa="${chain}_adapters.fasta" ow qin=33
	"""

}

process TRIM_PRIMERS {

	tag "${chain}"
	publishDir params.filtered_reads, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(chain), path(reads), val(primer_seqs), path(adapters)

	output:
	tuple val(chain), path("${chain}_chain_trimmed.fastq.gz")

	script:
	seq_patterns = primer_seqs
		.toString()
		.replace("[", "")
		.replace("]", "")
		.replace("'", "")
		.replace(" ", "")
		.replace(",", "\n")
	"""
	echo "${seq_patterns}" > primers.fasta
	bbduk.sh in=`realpath ${reads}` out=${chain}_chain_trimmed.fastq.gz \
	ref=primers.fasta,`realpath ${adapters}` \
	ktrim=r k=19 mink=11 hdist=2 \
	minlength=${params.min_len} maxlength=${params.max_len} \
	qin=33 threads=${task.cpus}
	"""

}

process CONVERT_TO_FASTA {

	tag "${chain}"

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	tuple val(chain), path(reads)

	output:
	tuple val(chain), path("${chain}_chain.fasta")

	script:
	"""
	seqkit fq2fa `realpath ${reads}` > ${chain}_chain.fasta
	"""

}

process CLUSTER_READS {

	tag "${chain}"
	publishDir params.consensus_seqs, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(chain), path(fasta)

	output:
	tuple val(chain), path("${chain}_consensus")

	script:
	def amb_flag = params.use_ambiguous ? '-amb' : ''
	"""
	amplicon_sorter.py \
	-i ${fasta} \
	-o ${chain}_consensus \
	-min ${params.min_len} -max ${params.max_len} \
	-sg ${params.similar_genes} \
	-ss ${params.similar_species} \
	-sc ${params.similar_consensus} \
	-ldc ${params.length_diff_consensus} \
	${amb_flag} \
	-ar -maxr ${params.max_reads} -ra -np ${task.cpus}
	"""

}

process BUILD_IGBLAST_DB {

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	path fastas

	output:
	path "bovine_ig_db*"

	script:
	"""
	# Combine all germline fastas
	cat *.fasta > combined_germlines.fasta

	# Format for IgBLAST (remove gaps, standardize headers)
	sed 's/\\./-/g' combined_germlines.fasta | \
	awk '/^>/{print; next}{gsub(/\\./, ""); print}' > bovine_ig_db

	# Build BLAST database
	makeblastdb -parse_seqids -dbtype nucl -in bovine_ig_db
	"""

}

process ANNOTATE_IGBLAST {

	tag "${chain}"
	publishDir params.annotations, mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	input:
	path db_files
	tuple val(chain), path(consensus_dir)

	output:
	tuple val(chain), path("${chain}_igblast.tsv"), path(consensus_dir)

	script:
	"""
	# Find consensus fasta files
	consensus_fasta=\$(find ${consensus_dir} -name "*.fasta" -o -name "*consensus*.fa" | head -1)

	if [ -z "\$consensus_fasta" ]; then
		# If no fasta found, create empty output
		echo "No consensus sequences found" > ${chain}_igblast.tsv
	else
		# Run IgBLAST with custom bovine database
		igblastn \
		-germline_db_V bovine_ig_db \
		-germline_db_J bovine_ig_db \
		-germline_db_D bovine_ig_db \
		-auxiliary_data optional_file/human_gl.aux \
		-query "\$consensus_fasta" \
		-outfmt "7 std qseq sseq" \
		-out ${chain}_igblast.tsv \
		|| echo "IgBLAST completed with warnings" > ${chain}_igblast.tsv
	fi
	"""

}

process PARSE_ANNOTATIONS {

	tag "${chain}"
	publishDir params.annotations, mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	input:
	tuple val(chain), path(igblast_out), path(consensus_dir)

	output:
	tuple val(chain), path("${chain}_annotations.tsv"), path("${chain}_cdr3.fasta")

	script:
	"""
	parse_igblast.py \
	--input ${igblast_out} \
	--consensus_dir ${consensus_dir} \
	--chain ${chain} \
	--output_tsv ${chain}_annotations.tsv \
	--output_cdr3 ${chain}_cdr3.fasta
	"""

}

process COLLECT_STATS {

	publishDir params.reports, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	path annotations

	output:
	path "summary_stats.tsv"

	script:
	"""
	echo -e "chain\tnum_sequences\tnum_unique_v\tnum_unique_j\tavg_cdr3_len" > summary_stats.tsv

	for f in *_annotations.tsv; do
		chain=\$(basename "\$f" _annotations.tsv)
		if [ -s "\$f" ]; then
			num_seqs=\$(tail -n +2 "\$f" | wc -l)
			echo -e "\${chain}\t\${num_seqs}\tNA\tNA\tNA" >> summary_stats.tsv
		fi
	done
	"""

}

process GENERATE_REPORT {

	publishDir params.reports, mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	input:
	path stats

	output:
	path "*.pdf", optional: true
	path "report_summary.txt"

	script:
	"""
	echo "Bovine IgG Repertoire Analysis Report" > report_summary.txt
	echo "=====================================" >> report_summary.txt
	echo "" >> report_summary.txt
	cat ${stats} >> report_summary.txt
	echo "" >> report_summary.txt
	echo "Analysis completed: \$(date)" >> report_summary.txt
	"""

}

process COLLECT_CONSENSUS_STATS {

	publishDir params.reports, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	path consensus_dirs

	output:
	path "summary_stats.tsv"

	script:
	"""
	echo -e "chain\tnum_consensus_sequences\ttotal_reads" > summary_stats.tsv

	for dir in */; do
		chain=\$(basename "\$dir" _consensus)
		if [ -d "\$dir" ]; then
			# Count consensus sequences
			num_seqs=\$(find "\$dir" -name "*.fasta" -exec grep -c "^>" {} + 2>/dev/null | awk -F: '{sum+=\$2} END {print sum}' || echo 0)
			echo -e "\${chain}\t\${num_seqs}\tNA" >> summary_stats.tsv
		fi
	done
	"""

}

// --------------------------------------------------------------- //

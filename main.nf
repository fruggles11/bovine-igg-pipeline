#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	// Auto-detect all barcode directories
	ch_barcodes = Channel
		.fromPath( "${params.fastq_dir}/barcode*/", type: 'dir' )
		.map { dir -> tuple( dir.getName(), dir ) }

	// Load primers for trimming (chain-specific, grouped by chain)
	ch_primers = Channel
		.fromPath( params.primer_table )
		.splitCsv( header: true )
		.map { row -> tuple( row.chain, row.primer_seq, row.differentiating ) }
		.filter { it[2] == "true" }
		.groupTuple( by: 0 )
		.map { chain, primers, diff -> tuple( chain, primers ) }

	// Load classification primers (separated by chain type)
	ch_class_primers = Channel
		.fromPath( params.primer_table )
		.splitCsv( header: true )
		.filter { row -> row.differentiating == "true" }

	heavy_primers = ch_class_primers
		.filter { row -> row.chain == 'heavy' }
		.map { row -> row.primer_seq }
		.collect()

	light_primers = ch_class_primers
		.filter { row -> row.chain == 'light' }
		.map { row -> row.primer_seq }
		.collect()

	// Germline gene files (if available)
	ch_germlines = Channel
		.fromPath( "${params.germline_dir}/*.fasta", checkIfExists: false )
		.collect()
		.ifEmpty( [] )

	// Workflow steps

	// Stage 1: Read Processing - Merge reads per barcode
	MERGE_READS(
		ch_barcodes
	)

	// Stage 2: Classify reads by primer sequence
	CLASSIFY_BY_PRIMER(
		MERGE_READS.out
			.combine( heavy_primers )
			.combine( light_primers )
	)

	// Combine heavy and light outputs, filter empty files
	ch_classified = CLASSIFY_BY_PRIMER.out.heavy
		.mix( CLASSIFY_BY_PRIMER.out.light )
		.filter { barcode_id, chain, fastq ->
			fastq.countFastq() >= params.min_reads
		}

	// Stage 3: Quality filter
	QUALITY_FILTER(
		ch_classified
	)

	FIND_ADAPTERS(
		QUALITY_FILTER.out
	)

	// Trim adapters only (no primer trimming)
	TRIM_PRIMERS(
		QUALITY_FILTER.out
			.join( FIND_ADAPTERS.out, by: [0, 1] )
	)

	// Stage 4: Clustering & Consensus
	CONVERT_TO_FASTA(
		TRIM_PRIMERS.out
	)

	CLUSTER_READS(
		CONVERT_TO_FASTA.out
	)

	// Stage 5: Extract majority consensus sequence per barcode/chain for Geneious import
	EXTRACT_MAJORITY_CONSENSUS(
		CLUSTER_READS.out
	)

	ch_majority = EXTRACT_MAJORITY_CONSENSUS.out
		.filter { barcode_id, chain, fasta -> fasta.size() > 0 }

	// Stage 6: Annotation (conditional on germline files being available)
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

		// Stage 6: Reporting
		COLLECT_STATS(
			PARSE_ANNOTATIONS.out.collect()
		)

		// Stage 7: VDJ annotation on majority consensus sequences
		BUILD_VDJ_DB(
			ch_germlines
		)

		ch_vdj_db = BUILD_VDJ_DB.out.first()

		ANNOTATE_MAJORITY_VDJ(
			ch_vdj_db,
			ch_majority
		)

		SUMMARIZE_VDJ(
			ANNOTATE_MAJORITY_VDJ.out.collect()
		)
	} else {
		// Skip annotation, just collect consensus stats
		// Extract only the path (third element) from the tuple before collecting
		COLLECT_CONSENSUS_STATS(
			CLUSTER_READS.out.map { barcode_id, chain, consensus_dir -> consensus_dir }.collect()
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
params.classified_reads = params.results + "/2_classified_reads"
params.filtered_reads = params.results + "/3_filtered_reads"
params.consensus_seqs = params.results + "/4_consensus_sequences"
params.majority_consensus = params.results + "/5_majority_consensus"
params.annotations = params.results + "/6_annotations"
params.vdj_annotations = params.results + "/6_vdj_annotations"
params.reports = params.results + "/7_reports"
// --------------------------------------------------------------- //



// PROCESS SPECIFICATION
// --------------------------------------------------------------- //

process MERGE_READS {

	tag "${barcode_id}"
	publishDir "${params.merged_reads}/${barcode_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), path(read_dir)

	output:
	tuple val(barcode_id), path("${barcode_id}_merged.fastq.gz")

	script:
	"""
	seqkit scat -j ${task.cpus} -f `realpath ${read_dir}` -o ${barcode_id}_merged.fastq.gz
	"""

}

process CLASSIFY_BY_PRIMER {

	tag "${barcode_id}"
	publishDir "${params.classified_reads}/${barcode_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), path(merged_reads), val(heavy_primer_list), val(light_primer_list)

	output:
	tuple val(barcode_id), val("heavy"), path("${barcode_id}_heavy.fastq.gz"), emit: heavy
	tuple val(barcode_id), val("light"), path("${barcode_id}_light.fastq.gz"), emit: light
	tuple val(barcode_id), val("unmatched"), path("${barcode_id}_unmatched.fastq.gz"), emit: unmatched

	script:
	def heavy_list = heavy_primer_list instanceof List ? heavy_primer_list : [heavy_primer_list]
	def light_list = light_primer_list instanceof List ? light_primer_list : [light_primer_list]
	def heavy_fasta = heavy_list.withIndex().collect { seq, i -> ">heavy_${i}\n${seq}" }.join('\n')
	def light_fasta = light_list.withIndex().collect { seq, i -> ">light_${i}\n${seq}" }.join('\n')
	"""
	# Create primer reference files in FASTA format
	cat << 'HEAVY_EOF' > heavy_primers.fasta
${heavy_fasta}
HEAVY_EOF

	cat << 'LIGHT_EOF' > light_primers.fasta
${light_fasta}
LIGHT_EOF

	# Step 1: Extract heavy chain reads (matching heavy primers)
	bbduk.sh in=`realpath ${merged_reads}` \
		outm=${barcode_id}_heavy.fastq.gz \
		outu=temp_not_heavy.fastq.gz \
		ref=heavy_primers.fasta \
		k=15 hdist=${params.primer_mismatch} \
		qin=33 threads=${task.cpus}

	# Step 2: From remaining reads, extract light chain reads
	bbduk.sh in=temp_not_heavy.fastq.gz \
		outm=${barcode_id}_light.fastq.gz \
		outu=${barcode_id}_unmatched.fastq.gz \
		ref=light_primers.fasta \
		k=15 hdist=${params.primer_mismatch} \
		qin=33 threads=${task.cpus}

	# Cleanup temp file
	rm -f temp_not_heavy.fastq.gz
	"""

}

process QUALITY_FILTER {

	tag "${barcode_id}_${chain}"
	publishDir "${params.filtered_reads}/${barcode_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), val(chain), path(reads)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_filtered.fastq.gz")

	script:
	"""
	seqkit seq \
	--min-len ${params.min_len} \
	--max-len ${params.max_len} \
	--min-qual ${params.min_qual} \
	--validate-seq \
	--threads ${task.cpus} \
	${reads} \
	-o ${barcode_id}_${chain}_filtered.fastq.gz
	"""

}

process FIND_ADAPTERS {

	tag "${barcode_id}_${chain}"

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	tuple val(barcode_id), val(chain), path(reads)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_adapters.fasta")

	script:
	"""
	bbmerge.sh in=`realpath ${reads}` outa="${barcode_id}_${chain}_adapters.fasta" ow qin=33
	"""

}

process TRIM_PRIMERS {

	tag "${barcode_id}_${chain}"
	publishDir "${params.filtered_reads}/${barcode_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), val(chain), path(reads), path(adapters)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_trimmed.fastq.gz")

	script:
	"""
	# Check if adapters file has real sequences (not just N's)
	if grep -v '^>' ${adapters} | grep -q '[ACGT]'; then
		bbduk.sh in=`realpath ${reads}` out=${barcode_id}_${chain}_trimmed.fastq.gz \
		ref=`realpath ${adapters}` \
		ktrim=r k=19 mink=11 hdist=2 \
		minlength=${params.min_len} maxlength=${params.max_len} \
		qin=33 threads=${task.cpus}
	else
		# No real adapters found, just copy input to output
		cp `realpath ${reads}` ${barcode_id}_${chain}_trimmed.fastq.gz
	fi
	"""

}

process CONVERT_TO_FASTA {

	tag "${barcode_id}_${chain}"

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	tuple val(barcode_id), val(chain), path(reads)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}.fasta")

	script:
	"""
	seqkit fq2fa `realpath ${reads}` > ${barcode_id}_${chain}.fasta
	"""

}

process CLUSTER_READS {

	tag "${barcode_id}_${chain}"
	publishDir "${params.consensus_seqs}/${barcode_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(barcode_id), val(chain), path(fasta)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_consensus")

	script:
	def amb_flag = params.use_ambiguous ? '-amb' : ''
	"""
	amplicon_sorter.py \
	-i ${fasta} \
	-o ${barcode_id}_${chain}_consensus \
	-min ${params.min_len} -max ${params.max_len} \
	-sg ${params.similar_genes} \
	-ss ${params.similar_species} \
	-sc ${params.similar_consensus} \
	-ldc ${params.length_diff_consensus} \
	${amb_flag} \
	-ar -maxr ${params.max_reads} -ra -np ${task.cpus}
	"""

}

process EXTRACT_MAJORITY_CONSENSUS {

	tag "${barcode_id}_${chain}"
	publishDir "${params.majority_consensus}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	tuple val(barcode_id), val(chain), path(consensus_dir)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_majority.fasta")

	script:
	"""
	# Find the top-level *_consensussequences.fasta (exclude group-specific *_N_consensussequences.fasta)
	consensus_file=\$(find ${consensus_dir} -name "*_consensussequences.fasta" \
		| grep -vE '_[0-9]+_consensussequences\\.fasta\$' | head -1)

	if [ -z "\$consensus_file" ]; then
		touch ${barcode_id}_${chain}_majority.fasta
		exit 0
	fi

	# Read counts are encoded in the header as (N) e.g. >consensus_barcode01_heavy_0_0(581)
	# Find the header with the highest read count
	grep "^>" "\$consensus_file" > headers.txt
	best_header=""
	best_count=0
	while IFS= read -r header; do
		count=\$(echo "\$header" | grep -oE '\\([0-9]+\\)' | tr -d '()')
		count=\${count:-0}
		if [ "\$count" -gt "\$best_count" ]; then
			best_count="\$count"
			best_header="\$header"
		fi
	done < headers.txt

	if [ -z "\$best_header" ]; then
		touch ${barcode_id}_${chain}_majority.fasta
		exit 0
	fi

	# Extract the sequence for that header and write with an informative name
	awk -v target="\$best_header" \
		-v new_header=">${barcode_id} chain=${chain} majority_consensus reads=\$best_count" '
		\$0 == target { found=1; next }
		/^>/ { if (found) exit }
		found { seq = seq \$0 }
		END { if (seq != "") { print new_header; print seq } }
	' "\$consensus_file" > ${barcode_id}_${chain}_majority.fasta
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

	tag "${barcode_id}_${chain}"
	publishDir "${params.annotations}/${barcode_id}", mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	input:
	path db_files
	tuple val(barcode_id), val(chain), path(consensus_dir)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_igblast.tsv"), path(consensus_dir)

	script:
	"""
	# Find consensus fasta files
	consensus_fasta=\$(find ${consensus_dir} -name "*.fasta" -o -name "*consensus*.fa" | head -1)

	if [ -z "\$consensus_fasta" ]; then
		# If no fasta found, create empty output
		echo "No consensus sequences found" > ${barcode_id}_${chain}_igblast.tsv
	else
		# Run IgBLAST with custom bovine database
		igblastn \
		-germline_db_V bovine_ig_db \
		-germline_db_J bovine_ig_db \
		-germline_db_D bovine_ig_db \
		-auxiliary_data optional_file/human_gl.aux \
		-query "\$consensus_fasta" \
		-outfmt "7 std qseq sseq" \
		-out ${barcode_id}_${chain}_igblast.tsv \
		|| echo "IgBLAST completed with warnings" > ${barcode_id}_${chain}_igblast.tsv
	fi
	"""

}

process PARSE_ANNOTATIONS {

	tag "${barcode_id}_${chain}"
	publishDir "${params.annotations}/${barcode_id}", mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	input:
	tuple val(barcode_id), val(chain), path(igblast_out), path(consensus_dir)

	output:
	tuple val(barcode_id), val(chain), path("${barcode_id}_${chain}_annotations.tsv"), path("${barcode_id}_${chain}_cdr3.fasta")

	script:
	"""
	parse_igblast.py \
	--input ${igblast_out} \
	--consensus_dir ${consensus_dir} \
	--chain ${chain} \
	--output_tsv ${barcode_id}_${chain}_annotations.tsv \
	--output_cdr3 ${barcode_id}_${chain}_cdr3.fasta
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
	echo -e "barcode\tchain\tnum_sequences\tnum_unique_v\tnum_unique_j\tavg_cdr3_len" > summary_stats.tsv

	for f in *_annotations.tsv; do
		# Extract barcode and chain from filename (format: barcode_chain_annotations.tsv)
		basename=\$(basename "\$f" _annotations.tsv)
		barcode=\$(echo "\$basename" | rev | cut -d'_' -f2- | rev)
		chain=\$(echo "\$basename" | rev | cut -d'_' -f1 | rev)
		if [ -s "\$f" ]; then
			num_seqs=\$(tail -n +2 "\$f" | wc -l)
			echo -e "\${barcode}\t\${chain}\t\${num_seqs}\tNA\tNA\tNA" >> summary_stats.tsv
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
	echo -e "barcode\tchain\tnum_consensus_sequences\ttotal_reads" > summary_stats.tsv

	for dir in */; do
		# Extract barcode and chain from directory name (format: barcode_chain_consensus)
		dirname=\$(basename "\$dir" _consensus)
		barcode=\$(echo "\$dirname" | rev | cut -d'_' -f2- | rev)
		chain=\$(echo "\$dirname" | rev | cut -d'_' -f1 | rev)
		if [ -d "\$dir" ]; then
			# Count consensus sequences
			num_seqs=\$(find "\$dir" -name "*.fasta" -exec grep -c "^>" {} + 2>/dev/null | awk -F: '{sum+=\$2} END {print sum}' || echo 0)
			echo -e "\${barcode}\t\${chain}\t\${num_seqs}\tNA" >> summary_stats.tsv
		fi
	done
	"""

}

process BUILD_VDJ_DB {

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2

	input:
	path fastas

	output:
	path "vdj_db_*"

	script:
	"""
	# Build separate V, D, J BLAST databases from IMGT-format germline files.
	# Files are identified by standard IMGT naming (e.g. IgHV, IgHD, IGHJ, IgKV, IgLJ, etc.).
	# Dots in sequence lines are IMGT alignment gaps and are stripped before DB construction.
	heavy_v=\$(ls *.fasta 2>/dev/null | grep -iE 'IGHV' | tr '\\n' ' ')
	heavy_d=\$(ls *.fasta 2>/dev/null | grep -iE 'IGHD' | tr '\\n' ' ')
	heavy_j=\$(ls *.fasta 2>/dev/null | grep -iE 'IGHJ' | tr '\\n' ' ')
	light_v=\$(ls *.fasta 2>/dev/null | grep -iE 'IG[KL]V' | tr '\\n' ' ')
	light_j=\$(ls *.fasta 2>/dev/null | grep -iE 'IG[KL]J' | tr '\\n' ' ')

	if [ -n "\$heavy_v" ]; then
		cat \$heavy_v | sed '/^>/!s/\\.//g' > heavy_V_nogaps.fasta
		makeblastdb -parse_seqids -dbtype nucl -in heavy_V_nogaps.fasta -out vdj_db_heavy_V
	fi
	if [ -n "\$heavy_d" ]; then
		cat \$heavy_d | sed '/^>/!s/\\.//g' > heavy_D_nogaps.fasta
		makeblastdb -parse_seqids -dbtype nucl -in heavy_D_nogaps.fasta -out vdj_db_heavy_D
	fi
	if [ -n "\$heavy_j" ]; then
		cat \$heavy_j | sed '/^>/!s/\\.//g' > heavy_J_nogaps.fasta
		makeblastdb -parse_seqids -dbtype nucl -in heavy_J_nogaps.fasta -out vdj_db_heavy_J
	fi
	if [ -n "\$light_v" ]; then
		cat \$light_v | sed '/^>/!s/\\.//g' > light_V_nogaps.fasta
		makeblastdb -parse_seqids -dbtype nucl -in light_V_nogaps.fasta -out vdj_db_light_V
	fi
	if [ -n "\$light_j" ]; then
		cat \$light_j | sed '/^>/!s/\\.//g' > light_J_nogaps.fasta
		makeblastdb -parse_seqids -dbtype nucl -in light_J_nogaps.fasta -out vdj_db_light_J
	fi

	ls vdj_db_* 2>/dev/null || touch vdj_db_placeholder
	"""

}

process ANNOTATE_MAJORITY_VDJ {

	tag "${barcode_id}_${chain}"
	publishDir { "${params.vdj_annotations}/${barcode_id}" }, mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	cpus 4

	input:
	path db_files
	tuple val(barcode_id), val(chain), path(majority_fasta)

	output:
	path("${barcode_id}_${chain}_vdj.airr.tsv")

	script:
	"""
	export IGDATA=/opt/ncbi-igblast-1.22.0

	if [ "${chain}" == "heavy" ]; then
		d_flag=""
		if ls vdj_db_heavy_D.* 2>/dev/null | head -1 > /dev/null 2>&1; then
			d_flag="-germline_db_D vdj_db_heavy_D"
		fi
		igblastn \\
			-germline_db_V vdj_db_heavy_V \\
			\$d_flag \\
			-germline_db_J vdj_db_heavy_J \\
			-query ${majority_fasta} \\
			-ig_seqtype Ig \\
			-organism human \\
			-num_threads ${task.cpus} \\
			-outfmt 19 \\
			-out ${barcode_id}_${chain}_vdj.airr.tsv \\
			|| echo -e "sequence_id\\tv_call\\td_call\\tj_call\\tproductive" > ${barcode_id}_${chain}_vdj.airr.tsv
	else
		igblastn \\
			-germline_db_V vdj_db_light_V \\
			-germline_db_J vdj_db_light_J \\
			-query ${majority_fasta} \\
			-ig_seqtype Ig \\
			-organism human \\
			-num_threads ${task.cpus} \\
			-outfmt 19 \\
			-out ${barcode_id}_${chain}_vdj.airr.tsv \\
			|| echo -e "sequence_id\\tv_call\\td_call\\tj_call\\tproductive" > ${barcode_id}_${chain}_vdj.airr.tsv
	fi
	"""

}

process SUMMARIZE_VDJ {

	publishDir params.reports, mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	input:
	path airr_files

	output:
	path "vdj_summary.tsv"

	script:
	"""
	echo -e "barcode\\tchain\\tv_call\\td_call\\tj_call\\tproductive\\tv_identity\\tjunction_aa_length" > vdj_summary.tsv

	for f in *_vdj.airr.tsv; do
		base=\$(basename "\$f" _vdj.airr.tsv)
		chain=\$(echo "\$base" | rev | cut -d'_' -f1 | rev)
		barcode=\$(echo "\$base" | rev | cut -d'_' -f2- | rev)

		if [ -s "\$f" ]; then
			awk -v barcode="\$barcode" -v chain="\$chain" '
			BEGIN { FS="\\t"; OFS="\\t" }
			NR==1 { for(i=1;i<=NF;i++) col[\$i]=i; next }
			NF > 1 {
				v_call    = (col["v_call"]     && \$col["v_call"]     != "") ? \$col["v_call"]     : "NA"
				d_call    = (col["d_call"]     && \$col["d_call"]     != "") ? \$col["d_call"]     : "NA"
				j_call    = (col["j_call"]     && \$col["j_call"]     != "") ? \$col["j_call"]     : "NA"
				prod      = (col["productive"] && \$col["productive"] != "") ? \$col["productive"] : "NA"
				v_id      = (col["v_identity"] && \$col["v_identity"] != "") ? \$col["v_identity"] : "NA"
				jaa_len   = (col["junction_aa"] && \$col["junction_aa"] != "") ? length(\$col["junction_aa"]) : "NA"
				print barcode, chain, v_call, d_call, j_call, prod, v_id, jaa_len
			}' "\$f" >> vdj_summary.tsv
		fi
	done
	"""

}

// --------------------------------------------------------------- //

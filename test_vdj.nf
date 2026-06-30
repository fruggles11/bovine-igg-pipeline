#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Minimal test pipeline for the VDJ annotation step only.
// Feed existing majority consensus FASTAs directly — skips all read-processing steps.
//
// Usage:
//   nextflow run test_vdj.nf -profile docker -work-dir /tmp/nf_work \
//     --majority_dir "/path/to/dir/with/*_majority.fasta files" \
//     --germline_dir "/path/to/bovine_germline" \
//     --results test_vdj_results

params.majority_dir = ""
params.germline_dir = "$projectDir/resources/germlines"
params.results      = "$launchDir/test_vdj_results"
params.vdj_annotations = params.results + "/vdj_annotations"
params.reports         = params.results + "/reports"


workflow {

	ch_germlines = Channel
		.fromPath( "${params.germline_dir}/*.fasta", checkIfExists: false )
		.collect()
		.ifEmpty( [] )

	// Parse barcode and chain from filenames: {barcode}_{chain}_majority.fasta
	ch_majority = Channel
		.fromPath( "${params.majority_dir}/*_majority.fasta" )
		.filter { it.size() > 0 }
		.map { fasta ->
			def stem  = fasta.getBaseName()                        // e.g. barcode03_heavy_majority
			def parts = stem.tokenize('_')
			def chain   = parts[-2]                                // heavy | light
			def barcode = parts[0..-3].join('_')                   // everything before _chain_majority
			tuple( barcode, chain, fasta )
		}

	BUILD_VDJ_DB( ch_germlines )

	ch_vdj_db = BUILD_VDJ_DB.out.first()

	ANNOTATE_MAJORITY_VDJ(
		ch_vdj_db,
		ch_majority
	)

	SUMMARIZE_VDJ(
		ANNOTATE_MAJORITY_VDJ.out.collect()
	)

}


process BUILD_VDJ_DB {

	errorStrategy 'retry'
	maxRetries 2

	input:
	path fastas

	output:
	path "vdj_db_*"

	script:
	"""
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
				v_call  = (col["v_call"]     && \$col["v_call"]     != "") ? \$col["v_call"]     : "NA"
				d_call  = (col["d_call"]     && \$col["d_call"]     != "") ? \$col["d_call"]     : "NA"
				j_call  = (col["j_call"]     && \$col["j_call"]     != "") ? \$col["j_call"]     : "NA"
				prod    = (col["productive"] && \$col["productive"] != "") ? \$col["productive"] : "NA"
				v_id    = (col["v_identity"] && \$col["v_identity"] != "") ? \$col["v_identity"] : "NA"
				jaa_len = (col["junction_aa"] && \$col["junction_aa"] != "") ? length(\$col["junction_aa"]) : "NA"
				print barcode, chain, v_call, d_call, j_call, prod, v_id, jaa_len
			}' "\$f" >> vdj_summary.tsv
		fi
	done
	"""

}

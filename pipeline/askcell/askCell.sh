#!/bin/bash
"""
# HSAPIEN
## DNA
HUMAN_RESOURCE=/mnt/BioPipeline/resource/Hsapien
GENOME=${HUMAN_RESOURCE}/GRCh38/genome/GRCh38.p13.genome.fa
GENOME_INDEX=${HUMAN_RESOURCE}/GRCh38/index/bwa/GRCh38.p13.genome.fa
SWIFT_SID_SNP=${HUMAN_RESOURCE}/GRCh38/panels/swift_sid/swift_sid.snp.grch38.bed
SWIFT_SID_AMPLICON=${HUMAN_RESOURCE}/GRCh38/panels/swift_sid/swift_sid.amplicons.grch38.bed
SWIFT_LUNG_AMPLICON=${HUMAN_RESOURCE}/GRCh38/panels/swift_lung/swift_lung.amplicons.grch38.bed
SWIFT_LUNG_TARGET=${HUMAN_RESOURCE}/GRCh38/panels/swift_lung/swift_lung.targets.grch38.bed
SWIFT_72G_AMPLICON=${HUMAN_RESOURCE}/GRCh38/panels/swift_72g/swift_72g.amplicons.grch38.bed
SWIFT_72G_TARGET=${HUMAN_RESOURCE}/GRCh38/panels/swift_72g/swift_72g.targets.grch38.bed
SWIFT_BRCA_AMPLICON=${HUMAN_RESOURCE}/GRCh38/panels/swift_brca/swift_brca.amplicons.grch38.bed
SWIFT_BRCA_TARGET=${HUMAN_RESOURCE}/GRCh38/panels/swift_brca/swift_brca.targets.grch38.bed
BULK_PURE_SID_SNP_DB=${HUMAN_RESOURCE}/GRCh38/pon/bulk_pure_cell_line_swift_sid/bulk_pure_cell_line_swift_sid.snp.rds
STAR_GENOME_INDEX=${HUMAN_RESOURCE}/GRCh38/index/star_151bp
rRNA_FA=${HUMAN_RESOURCE}/GRCh38/index/bwa_rRNA/rRNA.fasta
GNOMAD=${HUMAN_RESOURCE}/GRCh38/gnomad_v3/af-only-gnomad.hg38.vcf.gz
## RNAseq
STAR_MM_GENOME_INDEX=/mnt/BioPipeline/resource/genome_index/mm_star_151bp
GENCODE_GENE_BED=${HUMAN_RESOURCE}/GRCh38/annotation/gencode.v36.annotation.genes.bed6

# MM
## RNAseq
MM_RESOURCE=/mnt/BioPipeline/resource/Mouse
STAR_MM_GENOME_INDEX=${MM_RESOURCE}/GRCm39/index/mm_star_151bp
MM_rRNA_FA=${MM_RESOURCE}/GRCm39/index/bwa_rRNA/mm_rRNA.fasta

# Universal primer
P5="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
P7="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
TRIM_PARAM="-m 50 -q 20 --trim-n --max-n 5 --op-order GAWCQ"
"""


function error() {
		echo "[$(basename $1):$2:$3] ERROR: $4" 1>&2
	}

function info() {
		echo "[$1] INFO: $2"
	}

# check resource file
function check_resources () {

	echo "funciton to check resource files"
}

# check environment
function check_env () {
	echo "funciton to check resource files"
}

function check_file_exits () {
	local files=("$@")
	for f in "${files[@]}"; do
		if [[ ! -f ${f} ]]; then
			error $0 ${LINENO} "${FUNCNAME[0]}" "Cannot find the given file ${f}"
			exit 1
		fi
	done
}

function create_dirs () {
	local dirs=("$@")
	for d in "${dirs[@]}"; do
		if [[ ! -d ${d} ]]; then
			mkdir -p ${d}
		fi
	done
}

function trim_adapter () {
	local sample=$1
	local out_trimmed_dir=$2
	local out_log_dir=$3
	local molecule=$4
	local r1=$5
	local r2=$6

	local o1="${out_trimmed_dir}/${sample}.R1.trimmed.fastq.gz"
	local o2="${out_trimmed_dir}/${sample}.R2.trimmed.fastq.gz"
	local oinfo="${out_trimmed_dir}/${sample}.trimmed.info.txt"
	local ostats="${out_trimmed_dir}/${sample}.trimmed.stats.txt"
	local olog="${out_log_dir}/${sample}.adapter_trim.log"
	local adapter_trim_param=${TRIM_PARAM}
	case ${platform} in
		novaseq|nextseq|mn) adapter_trim_param="${adapter_trim_param} --nextseq-trim 20";;
		*) adapter_trim_param=${adapter_trim_param};;
	esac
	case ${workflow} in
		rna) adapter_trim_param="${adapter_trim_param}";;
		*) adapter_trim_param=${adapter_trim_param};;
	esac

	if [[ ${molecule} == "DNA" ]]; then
		local cmd="atropos trim --aligner insert ${adapter_trim_param} \
				 --no-cache-adapters \
				 -a ${P5} -A ${P7} \
				 -pe1 ${r1} -pe2 ${r2} \
				 -o ${o1} -p ${o2} \
				 --info-file ${oinfo} \
				 -T ${n_thread} \
				 > ${ostats} 2>${olog}"
	elif [[ ${molecule} == "RNA" ]]; then
		local cmd="atropos trim ${adapter_trim_param} \
				 -a ${P5} \
				 -a A{5}G{5} -a A{10}G{10} -a A{5}G{10} -a A{10}G{5} -a A{15}G{15} \
				 -a A{15} -a A{20} -a A{30} -a A{40} -a A{50} -a A{60} \
				 -a G{20} -a G{30} -a G{40} \
				 -se ${r1} \
				 -o ${o1}\
				 --info-file ${oinfo} \
				 -T ${n_thread} \
				 > ${ostats} 2>${olog}"
	else
		error $0 ${LINENO} ${FUNCNAME[0]} "Not supported molecules for adapter trimming"
		return 1
	fi
	eval ${cmd}
	if [[ $? -ne 0 ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Atropos returns non-zero code"
		error $0 ${LINENO} ${FUNCNAME[0]} "Please check log file ${olog}"
		return 1
	fi
	return 0
}

function map_reads () {
	local sample=$1
	local out_aln_dir=$2
	local out_log_dir=$3
	local molecules=$4
	local type=$5
	local r1=$6
	local r2=$7

	local olog="${out_log_dir}/${sample}.alignments.log"
	local rg="@RG\tID:${sample}\tSM:${sample}\tLB:${panel}\tPL:${platform}\tPU:run"

	if [[ ${molecule} == "DNA" ]]; then
		local bam="${out_aln_dir}/${sample}.so.bam"
		local otmp="${out_aln_dir}/${sample}.tmp."
		local cmd="bwa mem -M -t ${n_thread} -R \"${rg}\" ${GENOME_INDEX} ${r1} ${r2} 2>${olog} | \
				samtools view -bh | \
				samtools sort -T${otmp} -o ${bam} 2>>${olog}"
	elif [[ ${molecule} == "RNA" ]]; then
		local outp="${out_aln_dir}/${sample}/"
		local genome_index=${STAR_GENOME_INDEX}
		if [[ ${type} == "Ctrl" ]]; then
			genome_index=${STAR_MM_GENOME_INDEX}
		fi

		# because NAS does not support FIFO,
		# dumping tmp results to /tmp
		local cmd="STAR --runMode alignReads \
				   --genomeDir ${genome_index} \
				   --readFilesIn ${r1} \
				   --outSAMtype BAM SortedByCoordinate \
				   --quantMode GeneCounts \
				   --readFilesCommand \"zcat <\"  \
				   --runThreadN ${n_thread} \
				   --runDirPerm All_RWX \
				   --outFileNamePrefix ${outp} \
				   --outTmpDir /tmp/${sample} \
				   --outReadsUnmapped Fastx --outSAMunmapped Within --outSAMattributes All \
				   >${olog} 2>&1"
		local bam="${outp}/Aligned.sortedByCoord.out.bam"
	fi
	eval ${cmd}
	if [[ $? -ne 0 ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Alignment returns non-zero code"
		error $0 ${LINENO} ${FUNCNAME[0]} "Please check log file ${olog}"
		return 1
	fi

	if [[ ${workflow} == "wga" ]]; then
		#********************
		local mdup_bam="${out_aln_dir}/${sample}.so.mdup.bam"
		local mdup_metric="${out_aln_dir}/${sample}.mdup.metric.txt"
		cmd="picard "
		eval ${cmd}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "PICARD MarkDuplicates returns non-zero code"
			return 1
		fi
		bam=${mdup_bam}
		#********************
	fi
	cmd="samtools index ${bam}"
	eval ${cmd}
	if [[ $? -ne 0 ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Samtools index returns non-zero code"
		return 1
	fi

	return 0
}

function call_germlines () {
	local bam=$1
	local sample=$2
	local out_var_dir=$3
	local out_log_dir=$4

	local repo_dir=$(dirname "$0")
	local naive_call_py="${repo_dir}/naive_caller.py"
	if [[ ! -f ${naive_call_py} ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the naive_caller.py script in the current folder"
		return 1
	fi

	local out_cnt="${out_var_dir}/${sample}.swift_sid.cnt.csv"
	local olog="${out_log_dir}/${sample}.swift_sid.snp.log"
	local cmd="python3 ${naive_call_py} --nproc ${n_thread} --bam ${bam} --ref ${GENOME} --target ${SWIFT_SID_SNP} --sample ${sample} --out ${out_cnt} >${olog} 2>&1"
	eval ${cmd}
	if [[ $? -ne 0 ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Naive caller returns non-zero code"
		error $0 ${LINENO} ${FUNCNAME[0]} "Please check log file ${olog}"
		return 1
	fi
	return 0
}

function collect_dna_metrics () {
	local bam=$1
	local sample=$2
	local out_metrics_dir=$3
	local out_log_dir=$4
	local panel=$5

	local repo_dir=$(dirname "$0")
	local collect_metrics_py="${repo_dir}/collect_metrics.py"
	if [[ ! -f ${collect_metrics_py} ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the collect_metrics.py script in the current folder"
		return 1
	fi

	local target_bed=""
	if [ ${panel} = "swift_sid" ]; then
		target_bed=${SWIFT_SID_AMPLICON}
	elif [ ${panel} = "swift_lung" ]; then
		target_bed=${SWIFT_LUNG_TARGET}
	elif [ ${panel} = "swift_72g" ]; then
		target_bed=${SWIFT_72G_TARGET}
	elif [ ${panel} = "swift_brca" ]; then
		target_bed=${SWIFT_BRCA_TARGET}
	else
		error $0 ${LINENO} ${FUNCNAME[0]} "Unsupported panel"
		error $0 ${LINENO} ${FUNCNAME[0]} "Current supported panels: swift_sid, swift_lung, swift_72g"
		return 1
	fi
	local olog="${out_log_dir}/${sample}.collect_dna_metrics.log"
	local cmd="python3 ${collect_metrics_py} dna \
				--bam ${bam} \
				--sample ${sample} \
				--out ${out_metrics_dir} \
				--panel ${target_bed} \
				--genome ${GENOME} >${olog} 2>&1"
	eval ${cmd}
	if [[ $? -ne 0 ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Collect metrics returns non-zero code"
		error $0 ${LINENO} ${FUNCNAME[0]} "Please check log file ${olog}"
		return 1
	fi
	return 0
}

function collect_rna_metrics() {
	local sample=$1
	local out_trimmed_dir=$2
	local out_aln_dir=$3
	local out_metrics_dir=$4
	local out_log_dir=$5

	# do the adapter summary
	info ${FUNCNAME[0]} "Getting trimming stats for sample ${sample}"
	local trim_info_file=$(find ${out_trim_dir} -name "${sample}*info.txt")
	if [[ -z ${trim_info_file} ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the trimmed info file"
		return 1
	fi
	# 2nd column: -1=no adapter other value=adapter
	# 3rd column: when $2==-1, 3rd column is the read
	# 5th column: when $2!=-1, 5th column is the trimmed read
	# 6th adn 7th columns: when $2!=-1, 6th column is the adapter seq and 7th gives the sequence to the right of adapter seq
	local cmd="shuf -n 100000 ${trim_info_file} | awk -v FS=\"\t\" '{if(\$2 == -1) {n+=1; len_remain+=length(\$3); if(length(\$3) < 50) {nlt50+=1}} else{n_adapter+=1; if(length(\$5 < 50)) {nlt50+=1} len_remain+=length(\$5); adapter=\$6\"\"\$7; adapter_len+=length(adapter); nA+=gsub(\"A\", \"\", adapter); nG+=gsub(\"G\", \"\", adapter);}} END{print n_adapter/(n+n_adapter)*100, nlt50/(n+n_adapter)*100, len_remain/(n+n_adapter), adapter_len/n_adapter, nA/adapter_len*100, nG/adapter_len*100}'"
	read -r adapter_pct lt50_pct ave_read_len ave_adapter_len a_pct g_pct <<<$(eval ${cmd});

	# summarize alignment stats
	info ${FUNCNAME[0]} "Getting alignment stats for sample ${sample}"
	local final_aln_stats=$(find ${out_aln_dir}/ -name "Log.final.out")
	if [[ -z ${final_aln_stats} ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the ReadsPerGene.out.tab"
		return 1
	fi
	# please refer to Log.final.out
	cmd="awk '{if(\$0 ~ /Uniquely mapped reads number/) {n_unique=\$NF} else if(\$0 ~ /Number of input reads/) {n_tot=\$NF} else if(\$0 ~ /Number of reads mapped to multiple loci/) {n_multi=\$NF} else if(\$0 ~ /Number of reads mapped to too many loci/) {n_multi_many=\$NF} else if(\$0 ~ /Number of reads unmapped: too short/) {n_too_short=\$NF} else if(\$0 ~ /Number of reads unmapped: other/) {n_unmap_other=\$NF} else if(\$0 ~ /Mismatch rate per bas/) {mm_rate=\$NF}} END{print n_tot, n_unique, n_multi, n_multi_many, n_too_short, n_unmap_other, mm_rate}' ${final_aln_stats}"
	read -r n_tot n_uniq n_multi n_multi_many n_short n_unmap_other mm_rate <<<$(eval ${cmd});

	# get rRNA read counts
	info ${FUNCNAME[0]} "Mapping reads to rRNA sequence for sample ${sample}"
	unmapped_fq=$(find ${out_aln_dir}/ -name "Unmapped.out.mate1")
	if [[ -z ${unmapped_fq} ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the unmapped reads file"
		return 1
	fi
	local out_rRNA_metrics_dir="${out_metrics_dir}/rRNA_contam_alignments"
	create_dirs "${out_rRNA_metrics_dir}"
	local unmapped_to_rRNA_bam="${out_rRNA_metrics_dir}/${sample}.unmapped_to_rRNA.bam"
	local olog="${out_log_dir}/${sample}.rRNA.alignments.log"
	local rRNA_index=${rRNA_FA}
	if [[ ${type} == "Ctrl" ]]; then
		rRNA_index=${MM_rRNA_FA}
	fi
	local r1trim=$(find ${out_trim_dir} -name "${sample}*R1*fastq.gz")
	if [[ -z ${r1trim} ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the R1 read to do rRNA mapping"
		return 1
	fi
	cmd="bwa mem -M -t${n_thread} ${rRNA_index} ${r1trim} 2>${olog} | samtools view -bh - >${unmapped_to_rRNA_bam}"
	eval ${cmd}
	if [[ $? -ne 0 ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Alignments of unmapped reads to rRNA for sample ${sample} FAILED"
		error $0 ${LINENO} ${FUNCNAME[0]} "Please check the log file ${olog}"
		return 1
	fi
	info ${FUNCNAME[0]} "Getting # reads mapped to rRNA sequence for sample ${sample}"
	cmd="samtools view -F256 ${unmapped_to_rRNA_bam} | awk 'BEGIN{n_mapped=0} {if(\$2 != 4) {n_mapped+=1}} END{print n_mapped}'"
	read -r n_rRNA <<<$(eval ${cmd});

	# get stranded read counts
	info ${FUNCNAME[0]} "Getting stranded counts for sample ${sample}"
	local gene_count_file=$(find ${out_aln_dir}/ -name "ReadsPerGene.out.tab")
	if [[ -z ${gene_count_file} ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the ReadsPerGene.out.tab"
		return 1
	fi

	cmd="awk 'BEGIN{n_fwd=0; n_rev=0} (\$0 !~ /^N_/) {n_fwd+=\$3; n_rev+=\$4} END{print n_fwd, n_rev}' ${gene_count_file}"
	read -r n_fr_second n_fr_first <<<$(eval ${cmd});

	# get percentages
	local out_metrics="${out_metrics_dir}/${sample}.metrics.csv"
	local pct_uniq=$(awk -v n=${n_tot} -v m=${n_uniq} 'BEGIN {printf "%.2f", m/n*100}')
	local pct_multi=$(awk -v n=${n_tot} -v m=${n_multi} 'BEGIN {printf "%.2f", m/n*100}')
	local n_unmap=$(bc <<< ${n_multi_many}+${n_short}+${n_unmap_other})
	local pct_unmap=$(awk -v n=${n_tot} -v m=${n_unmap} 'BEGIN {printf "%.2f", m/n*100}')
	local pct_fwd=$(awk -v n=${n_tot} -v m=${n_fr_second} 'BEGIN {printf "%.2f", m/n*100}')
	local pct_rev=$(awk -v n=${n_tot} -v m=${n_fr_first} 'BEGIN {printf "%.2f", m/n*100}')
	local pct_rRNA=$(awk -v n=${n_tot} -v m=${n_rRNA} 'BEGIN {printf "%.2f", m/n*100}')
	local pct_mtRNA=$(awk -v n=${n_tot} -v m=${n_mt} 'BEGIN {printf "%.2f", m/n*100}')
	echo "SampleName,NumFrags,PCT_UniqMapped,PCT_MultiMapped,PCT_Unmapped,PCT_FwdStranded,PCT_RevStranded,MismatchRate,MeanFragLen,MeanAdapterLen,PCTAdapter,PCT_A,PCT_G,PCT_rRNA" > ${out_metrics}
	echo "${sample},${n_tot},${pct_uniq},${pct_multi},${pct_unmap},${pct_fwd},${pct_rev},${mm_rate},${ave_read_len},${ave_adapter_len},${adapter_pct},${a_pct},${g_pct},${pct_rRNA}" >> ${out_metrics}
	return 0
}

function process_one_sample () {
	local sample=$1
	local run_dir=$2
	local jira=$3
	local molecule=$4
	local type=$5
	local panel=$6

	# set up output directory
	out_jira_dir="${outdir}/${jira}"
	create_dirs ${out_jira_dir}

	local out_log_dir="${out_jira_dir}/log"
	local out_trim_dir="${out_jira_dir}/trimmed"
	sample_dirs=(${out_trim_dir} ${out_log_dir})
	create_dirs "${sample_dirs[@]}"

	# trim adapters
	sample=$(echo ${sample} | sed 's/_/-/g')
	info ${FUNCNAME[0]} "Start to process sample ${sample}"
	local odone="${out_log_dir}/${sample}.adapter_trim.done"
	if [[ ! -f ${odone} ]]; then
		# find sequencing data in the given run folder
		local seq_run_folder="${seq_run_folder}/${run_dir}"
		#local seq_run_folder="/Users/simo/Desktop/work/BulkRNA/Lexogen/${run_dir}"
		local r1=$(find ${seq_run_folder} -name "${sample}*R1*.fastq.gz")
		local r2=$(find ${seq_run_folder} -name "${sample}*R2*.fastq.gz")

		if [[ -z ${r1} ]]; then
			r1=$(find ${seq_run_folder} -name "${sample}.1*.fastq.gz")
			if [[ -z ${r1} ]]; then
				error $0 ${LINENO} ${FUNCNAME[0]} "Read1 is empty"
				return 1
			fi
		fi
		if [[ ${molecule} == "DNA" ]]; then
			if [[ -z ${r2} ]]; then
				error $0 ${LINENO} ${FUNCNAME[0]} "Read2 is empty"
				return 1
			fi
		fi

		local reads=(${r1} ${r2})
		check_file_exits "${reads[@]}"
		info ${FUNCNAME[0]} "Trimming adapters for sample ${sample}"
		#trim_adapter ${r1} ${r2} ${sample} ${out_trim_dir} ${out_log_dir} ${molecule}
		trim_adapter ${sample} ${out_trim_dir} ${out_log_dir} ${molecule} ${r1} ${r2}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Triming adapters for sample ${sample} FAILED"
			return 1
		fi
		touch ${odone}
	else
		info ${FUNCNAME[0]} "Previous results exist. Skip adapter trimming..."
	fi

	# reads alignment
	local out_aln_dir="${out_jira_dir}/alignments"
	create_dirs "${out_aln_dir}"
	odone="${out_log_dir}/${sample}.alignment.done"
	if [[ ! -f ${odone} ]]; then
		local r1trim=$(find ${out_trim_dir} -name "${sample}*R1*fastq.gz")
		local r2trim=$(find ${out_trim_dir} -name "${sample}*R2*fastq.gz")
		if [[ -z ${r1trim} ]] && [[ -f ${odone} ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Read1 trimmed is empty while trimming done file exists"
			return 1
		fi
		if [[ -z ${r2trim} ]] && [[ -f ${odone} ]] && [[ ${molecule} == "DNA" ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Read2 trimmed is empty while trimming done file exists"
			return 1
		fi
		info ${FUNCNAME[0]} "Aligning reads for sample ${sample}"
		#map_reads ${r1trim} ${r2trim} ${sample} ${out_aln_dir} ${out_log_dir} ${molecule}
		map_reads ${sample} ${out_aln_dir} ${out_log_dir} ${molecule} ${type} ${r1trim} ${r2trim}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Aligning reads for sample ${sample} FAILED"
			return 1
		fi
		touch ${odone}
	else
		info ${FUNCNAME[0]} "Previous results exist. Skip alignments..."
	fi

	# post alignment
	if [[ ${molecule} == "DNA" ]]; then
		local bam=$(find ${out_aln_dir} -name "${sample}.so.bam")
		# when panel is swift_sid, just call germlines
		if [[ ${panel} == "swift_sid" ]]; then
			local out_var_dir="${out_jira_dir}/germlines"
			create_dirs "${out_var_dir}"
			odone="${out_log_dir}/${sample}.call_germlines.done"
			if [[ ! -f ${odone} ]]; then
				info ${FUNCNAME[0]} "Calling germline variants for sample ${sample}"
				call_germlines ${bam} ${sample} ${out_var_dir} ${out_log_dir}
				if [[ $? -ne 0 ]]; then
					error $0 ${LINENO} ${FUNCNAME[0]} "Calling germline variants for sample ${sample} FAILED"
					return 1
				fi
				touch ${odone}
			else
				info ${FUNCNAME[0]} "Previous results exist. Skip calling germline variants..."
			fi
		fi

		local out_metrics_dir="${out_jira_dir}/metrics"
		create_dirs "${out_metrics_dir}"
		odone="${out_log_dir}/${sample}.collect_metrics.done"
		if [[ ! -f ${odone} ]]; then
			info ${FUNCNAME[0]} "Collecting metrics for sample ${sample}"
			collect_dna_metrics ${bam} ${sample} ${out_metrics_dir} ${out_log_dir} ${panel}
			if [[ $? -ne 0 ]]; then
				error $0 ${LINENO} ${FUNCNAME[0]} "Collecting metrics for sample ${sample} FAILED"
				return 1
			fi
			touch ${odone}
		else
			info ${FUNCNAME[0]} "Previous results exist. Skip collecting metrics..."
		fi
	elif [[ ${molecule} == "RNA" ]]; then
		local out_metrics_dir="${out_jira_dir}/metrics"
		odone="${out_log_dir}/${sample}.collect_rna_metrics.done"
		if [[ ! -f ${odone} ]]; then
			collect_rna_metrics ${sample} ${out_trim_dir} "${out_aln_dir}/${sample}" ${out_metrics_dir} ${out_log_dir}
			if [[ $? -ne 0 ]]; then
				error $0 ${LINENO} ${FUNCNAME[0]} "Collecting RNA metrics for sample ${sample} FAILED"
				return 1
			fi
			touch ${odone}
		fi
	else
		error $0 ${LINENO} ${FUNCNAME[0]} "Unsupported molecule type. Either DNA or RNA"
		return 1
	fi
	info ${FUNCNAME[0]} "Process sample ${sample} [DONE]"
	return 0
}

function check_sample_swap () {
	local t_sample=$1
	local n_sample=$2
	local jira=$3
	local target=$4

	local t_sample=$(echo ${t_sample} | sed 's/_/-/g')
	local n_sample=$(echo ${n_sample} | sed 's/_/-/g')
	local out_jira_dir="${outdir}/${jira}"
	local olog="${out_jira_dir}/log/${t_sample}.vs.${n_sample}.sample_swap_check.log"
	local odone="${out_jira_dir}/log/${t_sample}.vs.${n_sample}.sample_swap_check.done"
	if [ ! -f ${odone} ]; then
		local t_bam=$(find ${out_jira_dir} -name "${t_sample}*.bam")
		local n_bam=$(find ${out_jira_dir} -name "${n_sample}*.bam")
		if [ -z ${t_bam} ]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Failed to find BAM for tumor sample ${t_sample}"
			return 1
		fi
		if [ -z ${n_bam} ]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Failed to find BAM for normal sample ${n_sample}"
			return 1
		fi
		info ${FUNCNAME[0]} "Start checking potential swap between ${t_sample} and ${n_sample}"
		info ${FUNCNAME[0]} "Calling germline variants for ${t_sample} and ${n_sample}"
		local out_check_dir="${out_jira_dir}/swap_checker/${t_sample}.vs.${n_sample}"
		create_dirs ${out_check_dir}
		local ovcf="${out_check_dir}/${t_sample}.vs.${n_sample}.vcf.gz"
		local cmd="bcftools mpileup -f ${GENOME} -R ${target} \
					--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
					-Ou -q 30 -Q 25 -d 999999 --threads 8 --rf 2 \
					${t_bam} ${n_bam} 2>${olog} | bcftools call -m -Oz -v -o ${ovcf} 2>>${olog}"
		eval ${cmd}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "bcftools calling returns non-zero code"
			error $0 ${LINENO} ${FUNCNAME[0]} "Please check log file ${olog}"
			return 1
		fi
		cmd="tabix ${ovcf}"
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Indexing VCF returns non-zero code"
			return 1
		fi
		local repo_dir=$(dirname "$0")
		local check_swap_py="${repo_dir}/check_swap.py"
		if [[ ! -f ${check_swap_py} ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the check_swap.py script in the current folder"
			return 1
		fi
		info ${FUNCNAME[0]} "Calculating relatedness score"
		local out_swap_metric="${out_check_dir}/${t_sample}.vs.${n_sample}.sample_swap_check.metrics.tsv 1>>${olog} 2>&1"
		cmd="python3 ${check_swap_py} --vcf ${ovcf} --out ${out_swap_metric}"
		eval ${cmd}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Check sample swap returns non-zero code"
			error $0 ${LINENO} ${FUNCNAME[0]} "Please check log file ${olog}"
			return 1
		fi
		touch ${odone}
		info ${FUNCNAME[0]} "Checking potential swap between ${t_sample} and ${n_sample} [DONE]"
	else
		info ${FUNCNAME[0]} "Previous swap check results exist. Skip..."
	fi
}

function call_somatic () {
	local t_sample=$1
	local n_sample=$2
	local jira=$3
	local target=$4

	local t_sample=$(echo ${t_sample} | sed 's/_/-/g')
	local n_sample=$(echo ${n_sample} | sed 's/_/-/g')
	local out_jira_dir="${outdir}/${jira}"
	local olog="${out_jira_dir}/log/${t_sample}.vs.${n_sample}.somatics.log"
	local odone="${out_jira_dir}/log/${t_sample}.vs.${n_sample}.somatics.done"
	if [ ! -f ${odone} ]; then
		# find tumor and normal bams
		local t_bam=$(find ${out_jira_dir} -name "${t_sample}*.bam")
		local n_bam=$(find ${out_jira_dir} -name "${n_sample}*.bam")
		if [ -z ${t_bam} ]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Failed to find BAM for tumor sample ${t_sample}"
			return 1
		fi
		if [ -z ${n_bam} ]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Failed to find BAM for normal sample ${n_sample}"
			return 1
		fi
		info ${FUNCNAME[0]} "Calling somatic mutations for ${t_sample}"
		local out_soma_dir="${out_jira_dir}/somatics"
		local out_sample_dir="${out_soma_dir}/${t_sample}"
		create_dirs "${out_soma_dir}"
		create_dirs "${out_sample_dir}"
		# using m2
		local ovcf="${out_sample_dir}/${t_sample}.somatics.raw.vcf"
		local obam="${out_sample_dir}/${t_sample}.soamtics.bamout"
		local cmd="gatk Mutect2 --minimum-mapping-quality 30 --min-base-quality-score 25 \
					--read-filter FragmentLengthReadFilter --max-fragment-length 1000 \
					--dont-use-soft-clipped-bases \
					--read-filter OverclippedReadFilter --filter-too-short 50 \
					-R ${GENOME} \
					--max-reads-per-alignment-start 0 \
					--linked-de-bruijn-graph \
					-L ${target} \
					-germline-resource ${GNOMAD} \
					-I ${t_bam} \
					-I ${n_bam} --normal-sample ${n_sample} \
					-O ${ovcf} --bamout ${obam} \
					>${olog} 2>&1"
		eval ${cmd}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Somatic calling returns non-zero code"
			error $0 ${LINENO} ${FUNCNAME[0]} "Please check log file ${olog}"
			return 1
		fi
		local ofiltvcf="${out_sample_dir}/${t_sample}.somatics.filtered.vcf"
		# specific to mutect2
		local stats="${out_sample_dir}/${t_sample}.somatics.raw.vcf.stats"
		if [ ! -f ${stats} ]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Failed to find the somatic stats file"
			return 1
		fi
		info ${FUNCNAME[0]} "Filtering somatic mutations for ${t_sample}"
		cmd="gatk FilterMutectCalls \
				-R ${GENOME} \
				--stats ${stats} \
				--unique-alt-read-count 20 \
				--min-reads-per-strand 2 \
				--max-events-in-region 2 \
				--min-allele-fraction 0.01 \
				-V ${ovcf} -O ${ofiltvcf} \
				>> ${olog} 2>&1"
		eval ${cmd}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Somatic filtering returns non-zero code"
			error $0 ${LINENO} ${FUNCNAME[0]} "Please check log file ${olog}"
			return 1
		fi
		# bgzip and tabix
		touch ${odone}
	else
		info ${FUNCNAME[0]} "Previous somatic results exist. Skip..."
	fi

	return 0
}

function aggregate_and_summarize () {
	local jiras=$1

	for j in ${jiras[@]}; do
		info ${FUNCNAME[0]} "Start summarizing ${j}"
		local out_jira_dir="${outdir}/${j}"
		if [[ ! -d ${out_jira_dir} ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the output JIRA directory"
			exit 1
		fi
		local out_log_dir="${out_jira_dir}/log"
		local olog="${out_log_dir}/aggregate_and_summary.log"
		local odone="${out_log_dir}/aggregate_and_summary.done"
		if [[ ! -f ${odone} ]]; then
			local out_log_dir="${out_jira_dir}/log"
			local repo_dir=$(dirname "$0")
			local summary_R="${repo_dir}/summarize.R"
			if [[ ! -f ${summary_R} ]]; then
				error $0 ${LINENO} ${FUNCNAME[0]} "Cannot find the collect_metrics.py script in the current folder"
				return 1
			fi
			local cmd=""
			if [[ ${workflow} == "rna" ]]; then
				cmd="Rscript ${summary_R} \
					 --sheet ${samplesheet} \
					 --wdir ${out_jira_dir} \
					 --workflow ${workflow} \
					 --genes ${GENCODE_GENE_BED} "#\
					 #>${olog} 2>&1"
			elif [[ ${workflow} == "sid" ]] || [[ ${workflow} == "tgs" ]]; then
				cmd="Rscript ${summary_R} \
					 --sheet ${samplesheet} \
					 --wdir ${out_jira_dir} \
					 --workflow ${workflow} \
					 --bulk_db ${BULK_PURE_SID_SNP_DB}\
					 >${olog} 2>&1"
			fi
			eval ${cmd}
			if [[ $? -ne 0 ]]; then
				error $0 ${LINENO} ${FUNCNAME[0]} "Summarizing results FAILED"
				return 1
			fi
			touch ${odone}
			info ${FUNCNAME[0]} "Summarizing ${jira} [DONE]"
		else
			info ${FUNCNAME[0]} "Previous results for ${jira} exist. Skip summarizing..."
		fi
	done
	info ${FUNCNAME[0]} "Aggregating processed data and plot [DONE]"
	return 0

}

function run_rna () {
	local tmp_sheet=$(echo ${samplesheet} | sed 's/\.tsv/.tmp.tsv/g')
	tail -n +2 ${samplesheet} > ${tmp_sheet}
	echo "" >> ${tmp_sheet}
	local jiras=()
	while IFS=$'\t' read -r sample run_dir type well spikeinprop live input_cell input_dna instrument spikein_sample background jira revp fwdp il lb direct note; do
		if [ -z ${sample} ]; then
			break
		fi
		if [ ${sample} = "SampleName" ]; then
			continue
		fi
		echo ${sample} ${run_dir} ${type}
		echo ${panel} ${instrument} ${spikein_sample} ${background} ${jira}
		echo ${revp} ${fwdp} ${il} ${lb} ${direct}

		process_one_sample ${sample} ${run_dir} ${jira} "RNA" ${type} "wts"
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Processing sample ${sample} TERMINATED"
			exit 1
		fi
		jiras+=( ${jira} )
	done <"${tmp_sheet}"
	rm ${tmp_sheet}
	info ${FUNCNAME[0]} "Process all samples defined in the metasheet [DONE]"

	# the summary is done per JIRA cohort
	# all samples with the same JIRA have the same experimental purpose
	aggregate_and_summarize ${jiras}
	if [[ $? -ne 0 ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Summarize results TERMINATED"
		exit 1
	fi
	info ${FUNCNAME[0]} "Workflow ${workflow} [DONE]"
	return 0
}

function run_swift_sid () {
	local tmp_sheet=$(echo ${samplesheet} | sed 's/\.tsv/.tmp.tsv/g')
	tail -n +2 ${samplesheet} > ${tmp_sheet}
	echo "" >> ${tmp_sheet}
	local jiras=()
	while IFS=$'\t' read -r sample run_dir type well spikeinprop live input_cell input_dna panel instrument spikein_sample background jira note; do
		if [ -z ${sample} ]; then
			break
		fi
		if [ ${sample} = "SampleName" ]; then
			continue
		fi

		echo ${sample} "|" ${run_dir} ${type} ${well} ${spikeinprop} ${live} ${input_cell} ${input_dna}
		echo ${panel} ${instrument} ${spikein_sample} ${background} ${jira}

		process_one_sample ${sample} ${run_dir} ${jira} "DNA" ${type} ${panel}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Processing sample ${sample} TERMINATED"
			exit 1
		fi
		jiras+=( ${jira} )
	done <${tmp_sheet}
	info ${FUNCNAME[0]} "Process all samples defined in the metasheet [DONE]"
	rm ${tmp_sheet}

	aggregate_and_summarize ${jiras}
	if [[ $? -ne 0 ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Summarize results TERMINATED"
		exit 1
	fi
	info ${FUNCNAME[0]} "Workflow ${workflow} [DONE]"
	return 0
}

function run_tgs () {
	local tmp_sheet=$(echo ${samplesheet} | sed 's/\.tsv/.tmp.tsv/g')
	tail -n +2 ${samplesheet} > ${tmp_sheet}
	echo "" >> ${tmp_sheet}
	local jiras=()
	while IFS=$'\t' read -r t_sample t_run_dir type cell_type disease spikeinprop live t_input t_input n_sample n_run_dir n_input n_input jira lb wga panel instrument note; do
		if [ -z ${t_sample} ]; then
			break
		fi
		if [ ${t_sample} = "SampleName" ]; then
			continue
		fi

		# process tumor
		info ${FUNCNAME[0]} "----------------------------------------------------"
		process_one_sample ${t_sample} ${t_run_dir} ${jira} "DNA" ${type} ${panel}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Processing sample ${sample} TERMINATED"
			exit 1
		fi
		# process normal
		process_one_sample ${n_sample} ${n_run_dir} ${jira} "DNA" ${type} ${panel}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Processing sample ${sample} TERMINATED"
			exit 1
		fi

		target=""
		if [ ${panel} = "swift_72g" ]; then
			target=${SWIFT_72G_TARGET}
		elif [ ${panel} = "swift_lung" ]; then
			target=${SWIFT_LUNG_TARGET}
		elif [ ${panel} = "swift_brca" ]; then
			target=${SWIFT_BRCA_TARGET}
		else
			error $0 ${LINENO} ${FUNCNAME[0]} "Unrecognized panel name"
			error $0 ${LINENO} ${FUNCNAME[0]} "Current supported panel name for somatic calling: swift_72g, swift_lung, swift_brca"
			exit 1
		fi

		# check sample swap using something like somalier
		check_sample_swap ${t_sample} ${n_sample} ${jira} ${target}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Checking sample swap TERMINATED"
			exit 1
		fi

		# do somatic variant calling
		call_somatic ${t_sample} ${n_sample} ${jira} ${target}
		if [[ $? -ne 0 ]]; then
			error $0 ${LINENO} ${FUNCNAME[0]} "Calling somatic mutation TERMINATED"
			exit 1
		fi
		info ${FUNCNAME[0]} "----------------------------------------------------"

		jiras+=( ${jira} )
	done <${tmp_sheet}
	info ${FUNCNAME[0]} "Process all samples defined in the metasheet [done]"
	rm ${tmp_sheet}

	# aggregate results
	echo ${jiras}
	aggregate_and_summarize ${jiras}
	if [[ $? -ne 0 ]]; then
		error $0 ${LINENO} ${FUNCNAME[0]} "Summarize results TERMINATED"
		exit 1
	fi
	return 0
}

function get_abs_path () {
	local abs_path=$(cd -P -- "$(dirname -- "$1")" && printf '%s\n' "$(pwd -P)/$(basename -- "$1")")
	echo "${abs_path}"
}

function usage(){
	local program=$(basename $0)
	cat << EO
Usage: $program [options]
Options:
EO
    cat << EO | column -s\& -t
	-s or --samplesheet    & Path to the metasheet file to run the workflow [REQUIRED]
	-o or --outdir    & Path to the output base directory [REQUIRED]
	-w or --workflow    & Specify the workflow you want to run [REQUIRED] (sid,tgs,rna)
	-S or --seq_run    & Path to sequencing run folder where raw Fastq files can be found
	-p or --panel  & Specify the panel (not used)
	-P or --plt  & Specify the sequencing platform [mn]
	-t or --thread    & number of process [OPTIONAL] [8]
	--overwrite    & Force to overwrite previous results (not implemented yet)
EO
}

# main function below

n_thread=8
samplesheet=
workflow=
panel=
outdir=
platform="mn"
overwrite=false
seq_run_folder="/mnt/Data/MiniSeq"

while [ $# -gt 0 ]; do
	case $1 in
		-h|--help)
			usage
			exit 0;;
		-s|--samplesheet)
			shift; samplesheet=$(get_abs_path $1);;
		-w|--workflow)
			shift; workflow=$1;;
		-p|--panel)
			shift; panel=$1;;
		-S|--seq_run)
			shift; seq_run_folder=$1;;
		-o|-outdir)
			shift; outdir=$(get_abs_path $1);;
		-P|-plt)
			shift; platform=$1;;
		-t|--thread)
			shift; n_thread=$1;;
		--overwrite)
			overwrite=true;;
		--)
			shift; break;;
		*)
			echo "Invalid option: $1" 1>&2
			usage; exit 1;
			;;
	esac
	shift
done

if [[ -z ${samplesheet} ]]; then
	error $0 ${LINENO} "askCell_main" "askCell requires a samplesheet to work with"
	exit 1
fi

if [[ -z ${workflow} ]]; then
	error $0 ${LINENO} "askCell_main" "askCell needs to know what workflow you ask for"
	exit 1
fi

case ${workflow} in
	sid|tgs|wga|rna) workflow=${workflow};;
	*) error $0 ${LINENO} "askCell_main" "askCell yet to implement the provided workflow \"${workflow}\"";
		error $0 ${LINENO} "askCell_main" "askCell supports 3 workflows: sid, tgs, and wga";
		exit 1;;
esac

case ${platform} in
	novaseq|nextseq|mn) platform=${platform};;
	*) error $0 ${LINENO} "askCell_main" "askCell yet to support sequencing data from the given platform \"${platform}\"";
		error $0 ${LINENO} "askCell_main" "askCell supports sequencing data from 3 platforms: mn, nextseq, novaseq";
		exit 1;;
esac

#case ${panel} in
#	swift_sid|swift_lung|swift_72g|swift_56g|wts|wgs) panel=${panel};;
#	*) error $0 ${LINENO} "askCell_main" "askCell yet to support sequencing data from the given panel \"${panel}\"";
#		error $0 ${LINENO} "askCell_main" "askCell supports sequencing data from 4 panels: swift_sid, swift_lung, swift_72g, swift_56g";
#		exit 1;;
#esac

if [[ -z ${outdir} ]]; then
	error $0 ${LINENO} "askCell_main" "askCell requires a output directory to dump results"
	exit 1
fi

if [[ ! -d ${outdir} ]]; then
	mkdir -p ${outdir}
fi

check_file_exits "${samplesheet}"

#panel=""
#case ${workflow} in
#	sid) panel="swift_sid";;
#	tgs) panel="swift_lung";;
#	*) panel="";;
#esac

case ${workflow} in
	sid) run_swift_sid;;
	tgs) run_tgs;;
	wga) run_wga;;
	rna) run_rna;;
	*) error $0 ${LINENO} "askCell_main" "askCell yet to implement the provided workflow \"${workflow}\"";
	   exit 1;;
esac

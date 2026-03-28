rieux_short_reads_mapping() {
    local MAPQ=0
    local HOST_EXCLUDE=""
    local CONTROL_REF=""
    local INPUT_DIR=""
    local REF_MAP=""
    local GFF_FILE=""
    local OUT_DIR="$PWD/resultados_pipeline"
    local PARTITION=""
    local THREADS=""
    local MEM=""

    local OPTIND opt
    while getopts "i:r:g:h:c:m:o:p:t:M:" opt; do
        case $opt in
            i) INPUT_DIR=$(realpath "$OPTARG") ;;
            r) REF_MAP=$(realpath "$OPTARG") ;;
            g) GFF_FILE=$(realpath "$OPTARG") ;;
            h) HOST_EXCLUDE=$(realpath "$OPTARG") ;;
            c) CONTROL_REF=$(realpath "$OPTARG") ;;
            m) MAPQ="$OPTARG" ;;
            o) OUT_DIR=$(realpath "$OPTARG") ;;
            p) PARTITION="$OPTARG" ;;
            t) THREADS="$OPTARG" ;;
            M) MEM="$OPTARG" ;;
            *) echo "Uso inválido."; return 1 ;;
        esac
    done

    if [ -z "$INPUT_DIR" ] || [ -z "$REF_MAP" ] || [ -z "$GFF_FILE" ]; then
        echo -e "\n❌ ERRO: Faltam argumentos obrigatórios."
        echo "Uso: rieux_short_reads_mapping -i <pasta_fastq> -r <alvo.fasta> -g <alvo.gff3> [opcionais]"
        return 1
    fi

    local LOG_DIR="$OUT_DIR/logs"
    mkdir -p "$LOG_DIR" "$OUT_DIR"/{01_qc_pre,02_trimmed,03_qc_post,04_host_depletion,05_mapping_target,07_multiqc_report,ref_index,final_reports}

    echo -e "\n🔍 Lendo diretório de FASTQs..."
    find "$INPUT_DIR" -type f \( -name "*_R1_*.fastq.gz" -o -name "*_1.fastq.gz" \) | sort > "$LOG_DIR/lista_r1.txt"
    local NUM_SAMPLES=$(wc -l < "$LOG_DIR/lista_r1.txt")

    if [ "$NUM_SAMPLES" -eq 0 ]; then
        echo "❌ ERRO: Nenhum arquivo R1 detectado em $INPUT_DIR"
        return 1
    fi

    # ---------------------------------------------------------
    # 2. IA de Recursos com Trava de Segurança
    # ---------------------------------------------------------
    echo -e "🧠 Calculando recursos (Array Architecture)..."
    
    local HOST_MB=$([ -n "$HOST_EXCLUDE" ] && du -mL "$HOST_EXCLUDE" 2>/dev/null | cut -f1 || echo 0)
    local TARGET_MB=$(du -mL "$REF_MAP" 2>/dev/null | cut -f1)
    local CONTROL_MB=$([ -n "$CONTROL_REF" ] && du -mL "$CONTROL_REF" 2>/dev/null | cut -f1 || echo 0)
    
    local PEAK_FASTA_MB=$HOST_MB
    local COMPETITIVE_MB=$(( TARGET_MB + CONTROL_MB ))
    [ "$COMPETITIVE_MB" -gt "$PEAK_FASTA_MB" ] && PEAK_FASTA_MB=$COMPETITIVE_MB

    local CALC_MEM_GB
    local CALC_THREADS
    local AUTO_PARTITION="cpu"

    if [ -z "$PEAK_FASTA_MB" ] || [ "$PEAK_FASTA_MB" -eq 0 ]; then
        echo "⚠️  AVISO: Falha ao ler tamanho dos genomas. Ativando TRAVA DE SEGURANÇA (Modo Array)!"
        CALC_MEM_GB=40  
        CALC_THREADS=8  
    else
        CALC_MEM_GB=$(( (PEAK_FASTA_MB * 15 / 1024) + 8 ))
        [ $CALC_MEM_GB -lt 16 ] && CALC_MEM_GB=16
        CALC_THREADS=8 
        
        if [ $CALC_MEM_GB -gt 250 ]; then
            AUTO_PARTITION="fat"
            [ $CALC_MEM_GB -gt 4000 ] && CALC_MEM_GB=4000
        fi
    fi

    PARTITION=${PARTITION:-$AUTO_PARTITION}
    THREADS=${THREADS:-$CALC_THREADS}
    MEM=${MEM:-${CALC_MEM_GB}G}

    echo "=========================================================="
    echo "🚀 RIEUX SHORT READS MAPPING (SLURM ARRAYS - 3 STAGES)"
    echo "=========================================================="
    echo "⚙️  Estratégia    : Index -> Array ($NUM_SAMPLES tasks) -> Merge"
    echo "⚙️  Recursos/Task : Partição=[$PARTITION] | Threads=[$THREADS] | Mem=[$MEM]"
    [ -n "$CONTROL_REF" ] && echo "⚖️  Ref Controle  : $CONTROL_REF (Mapeamento Competitivo)"
    echo "📂 Diretório Out : $OUT_DIR"

    # ---------------------------------------------------------
    # JOB 0: Indexação Segura
    # ---------------------------------------------------------
    local INDEX_SCRIPT="${LOG_DIR}/00_index_ref.sbatch"
    local TARGET_INDEX="$REF_MAP"
    if [ -n "$CONTROL_REF" ]; then TARGET_INDEX="$OUT_DIR/ref_index/combined_target.fasta"; fi

    cat << 'EOF' > "$INDEX_SCRIPT"
#!/bin/bash
#SBATCH --job-name=rx_index
#SBATCH --output=logs/index_%j.out
#SBATCH --error=logs/index_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

REF_MAP=$1
CONTROL_REF=$2
OUT_DIR=$3

cd "$OUT_DIR"
if [ -n "$CONTROL_REF" ]; then
    module load bwa-mem2/2.3
    cat "$REF_MAP" "$CONTROL_REF" > ref_index/combined_target.fasta
    bwa-mem2 index ref_index/combined_target.fasta
    echo "✅ Índice competitivo gerado com sucesso."
fi
EOF

    # ---------------------------------------------------------
    # JOB 1: Array Processing
    # ---------------------------------------------------------
    local ARRAY_SCRIPT="${LOG_DIR}/01_array_process.sbatch"
    cat << 'EOF' > "$ARRAY_SCRIPT"
#!/bin/bash
#SBATCH --job-name=rx_array
#SBATCH --output=logs/array_%A_%a.out
#SBATCH --error=logs/array_%A_%a.err
#SBATCH --time=12:00:00

INPUT_DIR=$1
TARGET_INDEX=$2
HOST_EXCLUDE=$3
MAPQ=$4
OUT_DIR=$5
REF_MAP=$6

module load fastqc/0.12.1
module load fastp/0.23.4
module load bwa-mem2/2.3
module load samtools

R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$OUT_DIR/logs/lista_r1.txt")

if [[ "$R1" == *"_R1_"* ]]; then
    R2="${R1/_R1_/_R2_}"
    SAMPLE=$(basename "$R1" | sed -E 's/_S[0-9]+_L[0-9]+_R1_001.fastq.gz//' | sed 's/_R1_001.fastq.gz//')
elif [[ "$R1" == *"_1.fastq.gz" ]]; then
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    SAMPLE=$(basename "$R1" | sed 's/_1.fastq.gz//')
fi

echo "⚙ Iniciando Amostra: $SAMPLE"
TASK_SCRATCH="/scratch/$USER/rx_array_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TASK_SCRATCH"
cd "$TASK_SCRATCH"

fastp -w $SLURM_CPUS_PER_TASK -i "$R1" -I "$R2" -o ${SAMPLE}_R1_trim.fq.gz -O ${SAMPLE}_R2_trim.fq.gz -j ${SAMPLE}_fastp.json 2> fastp.log
fastqc -t $SLURM_CPUS_PER_TASK -q ${SAMPLE}_R1_trim.fq.gz ${SAMPLE}_R2_trim.fq.gz 2> fastqc.log

if [ -n "$HOST_EXCLUDE" ]; then
    bwa-mem2 mem -t $SLURM_CPUS_PER_TASK "$HOST_EXCLUDE" ${SAMPLE}_R1_trim.fq.gz ${SAMPLE}_R2_trim.fq.gz 2> bwa_host.log | samtools view -@ 4 -b -f 12 - > unmapped.bam
    samtools fastq -@ 4 -1 ${SAMPLE}_R1_clean.fq.gz -2 ${SAMPLE}_R2_clean.fq.gz -0 /dev/null -s /dev/null -n unmapped.bam 2> /dev/null
else
    mv ${SAMPLE}_R1_trim.fq.gz ${SAMPLE}_R1_clean.fq.gz
    mv ${SAMPLE}_R2_trim.fq.gz ${SAMPLE}_R2_clean.fq.gz
fi

# FILTRO BIOLÓGICO: Extrai apenas os nomes dos cromossomos do genoma alvo original
grep "^>" "$REF_MAP" | sed 's/>//' | awk '{print $1}' > target_contigs.txt

# MAPEAMENTO + ISOLAMENTO DO ALVO
bwa-mem2 mem -t $SLURM_CPUS_PER_TASK "$TARGET_INDEX" ${SAMPLE}_R1_clean.fq.gz ${SAMPLE}_R2_clean.fq.gz 2> bwa_target.log | \
    samtools view -@ 4 -h -q $MAPQ - | \
    awk 'NR==FNR{a[$1]; next} /^@/{print; next} $3 in a' target_contigs.txt - | \
    samtools view -@ 4 -b - | \
    samtools sort -@ 4 -o ${SAMPLE}_target_mapped_sorted.bam -

samtools index ${SAMPLE}_target_mapped_sorted.bam

cp -a ${SAMPLE}_fastp.json "$OUT_DIR/02_trimmed/"
cp -a *_fastqc.zip "$OUT_DIR/03_qc_post/" 2>/dev/null
cp -a ${SAMPLE}_R1_clean.fq.gz "$OUT_DIR/04_host_depletion/"
cp -a ${SAMPLE}_R2_clean.fq.gz "$OUT_DIR/04_host_depletion/"
cp -a ${SAMPLE}_target_mapped_sorted.bam "$OUT_DIR/05_mapping_target/"
cp -a ${SAMPLE}_target_mapped_sorted.bam.bai "$OUT_DIR/05_mapping_target/"

rm -rf "$TASK_SCRATCH"
echo "✅ Amostra $SAMPLE finalizada e purificada."
EOF

    # ---------------------------------------------------------
    # JOB 2: Consolidação (Merge)
    # ---------------------------------------------------------
    local MERGE_SCRIPT="${LOG_DIR}/02_merge_tables.sbatch"
    cat << 'EOF' > "$MERGE_SCRIPT"
#!/bin/bash
#SBATCH --job-name=rx_merge
#SBATCH --output=logs/merge_%j.out
#SBATCH --error=logs/merge_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

REF_MAP=$1
GFF_FILE=$2
OUT_DIR=$3

cd "$OUT_DIR"
module load samtools
module load bedtools
module load multiqc

multiqc 01_qc_pre/ 02_trimmed/ 03_qc_post/ -o 07_multiqc_report/ -q
TSV_FILE="07_multiqc_report/rastreabilidade_reads.tsv"
echo -e "Amostra\t# input reads\t# reads after trimming\t# reads not host\t# reads target\t% reads target\tBases_Cobertas\tProf_Efetiva(X)\tMAPQ_Medio" > $TSV_FILE

mapfile -t BAMS < <(ls 05_mapping_target/*_target_mapped_sorted.bam 2>/dev/null)
declare -a PROCESSED_SAMPLES

for bam in "${BAMS[@]}"; do
    SAMPLE=$(basename "$bam" | sed 's/_target_mapped_sorted.bam//')
    PROCESSED_SAMPLES+=("$SAMPLE")
    
    RAW_READS=$(grep -A 5 '"before_filtering":' 02_trimmed/${SAMPLE}_fastp.json 2>/dev/null | grep '"total_reads":' | tr -dc '0-9')
    TRIMMED_READS=$(grep -A 5 '"after_filtering":' 02_trimmed/${SAMPLE}_fastp.json 2>/dev/null | grep '"total_reads":' | tr -dc '0-9')
    [ -z "$RAW_READS" ] && RAW_READS=0
    [ -z "$TRIMMED_READS" ] && TRIMMED_READS=0
    
    if [ -s "04_host_depletion/${SAMPLE}_R1_clean.fq.gz" ]; then
        NON_HOST_R1=$(zcat "04_host_depletion/${SAMPLE}_R1_clean.fq.gz" | wc -l | awk '{print $1/4}')
        NON_HOST_READS=$((NON_HOST_R1 * 2))
    else
        NON_HOST_READS=$TRIMMED_READS
    fi

    # Como o BAM foi purificado, isso conta exatamente o alvo
    TARGET_READS=$(samtools view -c -F 260 "$bam")
    COV_STATS=$(samtools coverage "$bam" | grep -v "^#rname" | head -n 1)
    
    COV_BASES=0; EFF_DEPTH="0.00"; MEAN_MAPQ="0.0"
    if [ -n "$COV_STATS" ]; then
        GENOME_LEN=$(echo "$COV_STATS" | awk '{print $3}')
        COV_BASES=$(echo "$COV_STATS" | awk '{print $5}')
        MEAN_DEPTH=$(echo "$COV_STATS" | awk '{print $7}')
        MEAN_MAPQ=$(echo "$COV_STATS" | awk '{print $9}')
        if [ "$COV_BASES" -gt 0 ]; then EFF_DEPTH=$(awk -v md="$MEAN_DEPTH" -v gl="$GENOME_LEN" -v cb="$COV_BASES" 'BEGIN { printf "%.2f", (md * gl) / cb }'); fi
    fi
    
    PCT_TARGET=$(awk -v t=$NON_HOST_READS -v tr=$TARGET_READS 'BEGIN { if(t>0) printf "%.4f", (tr/t)*100; else print "0.0000" }')
    echo -e "${SAMPLE}\t${RAW_READS}\t${TRIMMED_READS}\t${NON_HOST_READS}\t${TARGET_READS}\t${PCT_TARGET}\t${COV_BASES}\t${EFF_DEPTH}\t${MEAN_MAPQ}" >> $TSV_FILE
done

cd 05_mapping_target/
grep "^>NC_" "$REF_MAP" | sed 's/>//' | awk '{print $1"\t1139633"}' > ref_genome.txt 2>/dev/null || samtools view -H ${PROCESSED_SAMPLES[0]}_target_mapped_sorted.bam | grep "@SQ" | awk '{split($2,a,":"); split($3,b,":"); print a[2]"\t"b[2]}' > ref_genome.txt

bedtools makewindows -g ref_genome.txt -w 500 > windows_500bp.bed
sort -k1,1 -k2,2n windows_500bp.bed > windows_sorted.bed
awk '$3=="CDS" || $3=="rRNA" || $3=="tRNA"' "$GFF_FILE" | sort -k1,1 -k4,4n > gff_filtered.gff3
bedtools map -a windows_sorted.bed -b gff_filtered.gff3 -c 3,9 -o collapse,collapse > windows_mapped.txt

awk -F'\t' -v OFS='\t' '{ if ($4 == ".") { feat = "intergenic"; nome = "NA"; } else { split($4, f, ","); feat = f[1]; if (match($5, /Name=[^;,]+/)) { nome = substr($5, RSTART+5, RLENGTH-5); } else if (match($5, /locus_tag=[^;,]+/)) { nome = substr($5, RSTART+10, RLENGTH-10); } else { nome = "NA"; } } print $1, $2, $3, feat, nome; }' windows_mapped.txt > windows_annotated.bed

TEMP_FILES=""
for SAMPLE in "${PROCESSED_SAMPLES[@]}"; do
    bam="${SAMPLE}_target_mapped_sorted.bam"
    bedtools bamtobed -i $bam | awk -F'\t' -v OFS='\t' '{mid=int(($2+$3)/2); print $1, mid, mid+1}' | bedtools coverage -a windows_annotated.bed -b - -counts | awk '{print $NF}' > "temp_${SAMPLE}.txt"
    TEMP_FILES="$TEMP_FILES temp_${SAMPLE}.txt"
done

paste windows_annotated.bed $TEMP_FILES | awk -F'\t' -v OFS='\t' '{ soma = 0; for(i=6; i<=NF; i++) soma += $i; if(soma > 0) print $0; }' > matriz_final_elegante.tsv

HEADER_SAMPLES=$(printf "%s\t" "${PROCESSED_SAMPLES[@]}" | sed 's/\t$//')
echo -e "Chr\tStart\tEnd\tFeature\tName\t$HEADER_SAMPLES" > cabecalho.txt
cat cabecalho.txt matriz_final_elegante.tsv > matriz_final_com_cabecalho.tsv

# Geração da Tabela de Blocos Contíguos
awk -F'\t' -v OFS='\t' 'NR==1 { print $0; num_cols = NF; next; } { if (prev_chr == "") { prev_chr = $1; prev_start = $2; prev_end = $3; prev_feat = $4; prev_name = $5; for(i=6; i<=NF; i++) sums[i] = $i; } else if ($1 == prev_chr && $4 == prev_feat && $5 == prev_name && $2 == prev_end) { prev_end = $3; for(i=6; i<=NF; i++) sums[i] += $i; } else { printf "%s\t%s\t%s\t%s\t%s", prev_chr, prev_start, prev_end, prev_feat, prev_name; for(i=6; i<=num_cols; i++) printf "\t%s", sums[i]; printf "\n"; prev_chr = $1; prev_start = $2; prev_end = $3; prev_feat = $4; prev_name = $5; for(i=6; i<=NF; i++) sums[i] = $i; } } END { if (prev_chr != "") { printf "%s\t%s\t%s\t%s\t%s", prev_chr, prev_start, prev_end, prev_feat, prev_name; for(i=6; i<=num_cols; i++) printf "\t%s", sums[i]; printf "\n"; } }' matriz_final_com_cabecalho.tsv > matriz_blocos_contiguos.tsv

# Geração do Super Resumo
awk -F'\t' -v OFS='\t' 'NR==1 { num_samples = NF - 5; for(i=6; i<=NF; i++) sample_names[i] = $i; next; } { feat = $4; if (!(feat in feat_seen)) { feat_seen[feat] = 1; feats[++num_feats] = feat; } for(i=6; i<=NF; i++) matrix[i, feat] += $i; } END { printf "Amostra"; for(f=1; f<=num_feats; f++) printf "\t%s", feats[f]; printf "\n"; for(i=6; i<=num_samples+5; i++) { printf "%s", sample_names[i]; for(f=1; f<=num_feats; f++) { val = matrix[i, feats[f]] + 0; printf "\t%s", val; } printf "\n"; } }' matriz_final_com_cabecalho.tsv > resumo_super_features.tsv

# Extração de Ouro: Copia as tabelas finais e apaga o lixo cru
cp ../07_multiqc_report/rastreabilidade_reads.tsv "$OUT_DIR/final_reports/"
cp ../07_multiqc_report/multiqc_report.html "$OUT_DIR/final_reports/"
cp matriz_blocos_contiguos.tsv "$OUT_DIR/final_reports/"
cp resumo_super_features.tsv "$OUT_DIR/final_reports/"

rm temp_*.txt cabecalho.txt ref_genome.txt windows_500bp.bed windows_sorted.bed gff_filtered.gff3 windows_mapped.txt windows_annotated.bed matriz_final_elegante.tsv matriz_final_com_cabecalho.tsv
echo "🏆 PIPELINE FINALIZADO COM SUCESSO!"
EOF

    # ---------------------------------------------------------
    # 6. Orquestração e Submissão via SLURM Dependencies
    # ---------------------------------------------------------
    cd "$OUT_DIR"
    echo -e "\n⏳ Orquestrando Jobs no SLURM..."
    
    local INDEX_ID=$(sbatch --parsable "$INDEX_SCRIPT" "$REF_MAP" "$CONTROL_REF" "$OUT_DIR")
    echo "   [ID: $INDEX_ID] 1. Job de Indexação submetido."
    
    local MAX_CONCURRENT=15
    local ARRAY_ID=$(sbatch --parsable --dependency=afterok:$INDEX_ID --array=1-${NUM_SAMPLES}%${MAX_CONCURRENT} --partition="$PARTITION" --cpus-per-task="$THREADS" --mem="$MEM" "$ARRAY_SCRIPT" "$INPUT_DIR" "$TARGET_INDEX" "$HOST_EXCLUDE" "$MAPQ" "$OUT_DIR" "$REF_MAP")
    echo "   [ID: $ARRAY_ID] 2. Job Array de Mapeamento ($NUM_SAMPLES tasks, máx $MAX_CONCURRENT simultâneas) engatilhado."
    
    local MERGE_ID=$(sbatch --parsable --dependency=afterok:$ARRAY_ID "$MERGE_SCRIPT" "$REF_MAP" "$GFF_FILE" "$OUT_DIR")
    echo "   [ID: $MERGE_ID] 3. Job de Consolidação engatilhado."
    
    echo -e "\n✅ Cadeia de execução completa! Monitore com: squeue -u \$USER"
}

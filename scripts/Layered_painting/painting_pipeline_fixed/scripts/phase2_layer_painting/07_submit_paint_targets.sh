#!/bin/bash
#
# Step 7: Paint targets (interactive or SLURM array)
#

set -e

LAYER=$1
CHR=$2

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname $(dirname $SCRIPT_DIR))"
CONFIG_FILE="$BASE_DIR/config/master_config.yaml"

# Parse config
OUTPUT_BASE=$(grep "output_base:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
EXECUTION_MODE=$(grep "execution_mode:" $CONFIG_FILE | awk '{print $2}' | tr -d '"')
PARTITION=$(grep "partition:" $CONFIG_FILE | awk '{print $2}' | tr -d '"')
TIME_LIMIT=$(grep "time_limit:" $CONFIG_FILE | awk '{print $2}' | tr -d '"')
MEM=$(grep "mem_per_job:" $CONFIG_FILE | awk '{print $2}' | tr -d '"')

LAYER_DIR="$OUTPUT_BASE/data/layer${LAYER}"
BATCH_DIR="$LAYER_DIR/target/batches/chr${CHR}"

# Get number of batches
N_BATCHES=$(cat $BATCH_DIR/n_batches.txt)

echo "Step 7: Painting targets for chr$CHR"
echo "  Execution mode: $EXECUTION_MODE"
echo "  Number of batches: $N_BATCHES"

# Create the painting job script
JOB_SCRIPT="$BATCH_DIR/paint_batch.sh"

cat > $JOB_SCRIPT << 'EOFSCRIPT'
#!/bin/bash

# Load modules (module load will fail gracefully if already loaded)
module load perl/5.38.0 2>/dev/null || true
module load gsl/2.5 2>/dev/null || true
module load bcftools/1.20 2>/dev/null || true
module load python/3.9.5 2>/dev/null || true

# Variables
LAYER={LAYER}
CHR={CHR}
BATCH=$1  # Batch number passed as argument
OUTPUT_BASE={OUTPUT_BASE}
LAYER_DIR=$OUTPUT_BASE/data/layer$LAYER
BATCH_DIR=$LAYER_DIR/target/batches/chr$CHR
RESULTS_DIR=$OUTPUT_BASE/results/layer$LAYER

# Paths
SPARSEPAINTER={SPARSEPAINTER}
PBWT={PBWT}
REF_PHASE=$LAYER_DIR/ref/chr${CHR}.ref.phase
GENETIC_MAP=$LAYER_DIR/genetic_maps/chr${CHR}.map
REF_POPFILE={REF_POPFILE}
LAMBDA=$(cat $LAYER_DIR/lambda/chr${CHR}_lambda.txt)

# Batch-specific files
BATCH_SAMPLES=$BATCH_DIR/batch${BATCH}_samples.txt
BATCH_VCF=$LAYER_DIR/target/chr${CHR}.target.vcf.gz
BATCH_PHASE=$BATCH_DIR/batch${BATCH}_chr${CHR}.phase
BATCH_NAMEFILE=$BATCH_DIR/batch${BATCH}_namefile.txt

# Create batch phase file using PBWT
echo "Processing batch $BATCH for chr$CHR ($(date))"
echo "  Extracting batch samples..."
$PBWT -readVcfGT $BATCH_VCF -selectSamples $BATCH_SAMPLES -writePhase $BATCH_PHASE

# Create namefile
cp $BATCH_SAMPLES $BATCH_NAMEFILE

# Run SparsePainter
echo "  Running SparsePainter..."
cd $BATCH_DIR

$SPARSEPAINTER \
    -reffile $REF_PHASE \
    -targetfile $BATCH_PHASE \
    -mapfile $GENETIC_MAP \
    -popfile $REF_POPFILE \
    -namefile $BATCH_NAMEFILE \
    -haplambda $LAMBDA \
    -prob \
    -chunklength \
    -probstore raw \
    -out chr${CHR}_batch${BATCH}

# Gzip outputs
gzip -f chr${CHR}_batch${BATCH}_prob.txt 2>/dev/null || true
gzip -f chr${CHR}_batch${BATCH}_chunklength.txt 2>/dev/null || true

# Clean up phase file to save space
rm -f $BATCH_PHASE

echo "  âœ“ Batch $BATCH complete for chr$CHR ($(date))"
EOFSCRIPT

# Replace placeholders
sed -i "s|{CHR}|$CHR|g" $JOB_SCRIPT
sed -i "s|{LAYER}|$LAYER|g" $JOB_SCRIPT
sed -i "s|{OUTPUT_BASE}|$OUTPUT_BASE|g" $JOB_SCRIPT

# Replace tool paths
SPARSEPAINTER=$(grep "sparsepainter:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
PBWT=$(grep "pbwt:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
REF_POPFILE="$BASE_DIR/config/samples/layer1/ref_popfile.txt"

sed -i "s|{SPARSEPAINTER}|$SPARSEPAINTER|g" $JOB_SCRIPT
sed -i "s|{PBWT}|$PBWT|g" $JOB_SCRIPT
sed -i "s|{REF_POPFILE}|$REF_POPFILE|g" $JOB_SCRIPT

chmod +x $JOB_SCRIPT

if [ "$EXECUTION_MODE" == "interactive" ]; then
    echo "  Running batches interactively (sequentially)..."
    echo "  This will take a while - consider running in screen/tmux"
    echo ""
    
    for BATCH_NUM in $(seq 1 $N_BATCHES); do
        echo "  Starting batch $BATCH_NUM/$N_BATCHES..."
        bash $JOB_SCRIPT $BATCH_NUM 2>&1 | tee -a $OUTPUT_BASE/logs/paint_layer${LAYER}_chr${CHR}_batch${BATCH_NUM}.log
    done
    
    echo "  âœ“ All batches complete for chr$CHR"
    
elif [ "$EXECUTION_MODE" == "slurm" ]; then
    echo "  Submitting SLURM array job..."
    
    # Create SLURM wrapper
    SLURM_SCRIPT="$BATCH_DIR/paint_slurm.sh"
    cat > $SLURM_SCRIPT << EOFSLURM
#!/bin/bash
#SBATCH --job-name=paint_${CHR}_${LAYER}
#SBATCH --output=$OUTPUT_BASE/logs/paint_layer${LAYER}_chr${CHR}_batch%a.out
#SBATCH --error=$OUTPUT_BASE/logs/paint_layer${LAYER}_chr${CHR}_batch%a.err
#SBATCH --partition=$PARTITION
#SBATCH --time=$TIME_LIMIT
#SBATCH --mem=$MEM
#SBATCH --cpus-per-task=1

bash $JOB_SCRIPT \$SLURM_ARRAY_TASK_ID
EOFSLURM
    
    sbatch --array=1-$N_BATCHES $SLURM_SCRIPT
    echo "  âœ“ Submitted $N_BATCHES jobs for chr$CHR"
    echo "  Monitor with: squeue -u \$USER"
else
    echo "  Error: Unknown execution mode: $EXECUTION_MODE"
    exit 1
fi

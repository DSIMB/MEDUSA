#! /bin/bash

# Show this for help
usage() {
cat << EOF
Usage:
        -i    | --seq         (Required)     Path to input Fasta sequence file. The path is relative to the project folder.
        -d    | --database    (Required)     Name of the database for HHBlits.
                                             Ex: UniRef30_2020_06
                                             If your path is '/path/to/Uniclust/UniRef30_2020_03_a3m.ffdata'
					     please provide: -d UniRef30_2020_03
        -o    | --outdir      (Required)     Path to the output directory.
	-c    | --cpus        (Optionnal)    Number of CPUs to use (for HHblits). Default is all virtual cores available.
	-m    | --memory      (Optionnal)    Maximum RAM to use in Gb (for HHblits). Default is all available.
        -h    | --help        (Optionnal)    Brings up this help
EOF
}

# Detect hardware based on OS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    NB_CPUS=$(getconf _NPROCESSORS_ONLN)
    MAX_MEMORY=$(expr $(free -g | awk '/^Mem:/{print $2}') - 2)
elif [[ "$OSTYPE" == "darwin"* ]]; then
    NB_CPUS=$(sysctl -n hw.ncpu)
    MAX_MEMORY=$(expr $(sysctl -n hw.memsize) / $((1024**3)) - 2 )
elif [[ "$OSTYPE" == "cygwin" ]]; then
    NB_CPUS=$(echo %NUMBER_OF_PROCESSORS%)
    MAX_MEMORY=4
fi

# Parse command line arguments
while [ "$1" != "" ]; do
    case $1 in
        -i | --seq )
            shift
            SEQ=$1
        ;;
        -d | --database )
            shift
            DATABASE=$1
        ;;
        -o | --outdir )
            shift
            OUTDIR=$1
        ;;
        -c | --cpus )
            shift
            NB_CPUS=$1
        ;;
        -m | --memory )
            shift
            MAX_MEMORY=$1
        ;;
        -h | --help )
            usage
            exit
        ;;
        * ) usage
            exit 1
    esac
    shift
done

# Check for valid sequence input
if [ ! -f /project/$SEQ ]; then
    printf "\nA valid input sequence file is required, provide it with the flag: -i sequence.fasta\nMake sure the path to the file is relative to your project folder.\n\n"
    exit
fi

# Check if database is accessible
if [ ! -f /database/$DATABASE"_a3m.ffindex" ]; then
    printf "\n The database $DATABASE seem unusable, please link the right database."
    printf "\nIf your path is '/path/to/database/UniRef30_2020_03_a3m.ffdata' provide: -d UniRef30_2020_03\n\n"
    exit
fi

# Check if the output directory exists
if [ ! -d /project/$OUTDIR ]; then
    printf "\n$OUTDIR does not exist, please provide a right output directory.\n\n"
    exit
fi

if [[ ! $NB_CPUS =~ ^[0-9]+$ ]]; then
    printf "\nThe number of CPUs argument should be an integer.\n\n"
    exit
fi

if [[ ! $MAX_MEMORY =~ ^[0-9]+$ ]]; then
    printf "\nThe memory argument should be an integer.\n\n"
    exit
fi

# Create a directory for the job output
export JOBDIR=$(mktemp -d -t MEDUSA.XXXX --suffix=-$(date +%Y%m%d%H%M%S) -p /project/$OUTDIR)
export ORIGINAL_OUTDIR_NAME=$OUTDIR/`basename $JOBDIR`

# Set global paths

# The project is mounted as docker bind volume 
# into /project directory in the docker container
export PROJECT=/project
export HHSUITE=$PROJECT/hh-suite
export HHBLITS=$HHSUITE/bin/hhblits
export HHFILTER=$HHSUITE/bin/hhfilter
# path to database is a mounted volume to /database when docker is launched
# the name of the database is given in command line argument
export DBHHBLITS=/database/$DATABASE
export OUTDIR=$JOBDIR
export SEQ=$PROJECT/$SEQ
export SCRIPTS=$PROJECT/scripts
export DATA=$PROJECT/data
export MODELS=$PROJECT/models
export WINDOW=15



#############
### START ###
#############


### Step 1. Create pssm data.
### HHblits is time consuming, probably remove the -realign_old_hits option. Probably decrease -max_filt and -realign_max.
printf "Run HHblits ... "
$HHBLITS -cpu $NB_CPUS -maxmem $MAX_MEMORY -maxfilt 10000 -diff inf -B 10000 -Z 10000 -e 0.0001 -cov 75 -realign_old_hits -realign_max 10000 -n 3 -i $SEQ -d $DBHHBLITS -oa3m $OUTDIR/job.a3m -o $OUTDIR/job.hhr 1>/dev/null 2> $OUTDIR/log_hhblits
printf "done\n"

printf "Run HHfilter ... "
$HHFILTER -v 2 -id 99 -neff 20 -qsc -30 -cov 75 -i $OUTDIR/job.a3m -o $OUTDIR/job\_filter.a3m 2>/dev/null
printf "done\n"

#perl $HHSUITE/scripts/reformat.pl -M -uc first a3m fas job\_filter.a3m pdb_mfasta/job.mfasta
# without insertions in the first (query) sequence 
printf "Reformat ... "
perl $HHSUITE/scripts/reformat.pl -r -M -uc first a3m fas $OUTDIR/job\_filter.a3m $OUTDIR/job.mfasta_woi 1>/dev/null 2>&1
printf "done\n"

printf "Transform alignment to frequence matrix ... "
$SCRIPTS/ali2freq-py3.py -first -al $OUTDIR/job.mfasta_woi -m $DATA/homstradfreq.txt -gapaa 1> $OUTDIR/job.aamtx_gaps 2>/dev/null
printf "done\n"

printf "Create AAindex and one-hot encodings ... "
### Step 2. Create aaindex data.
$SCRIPTS/fasta2vector_wgap $SEQ $DATA/selected_aaindex1_reformated_58_Van_Westen > $OUTDIR/job.aaindex

### Step 3. Create one-hot data.
$SCRIPTS/seq2onehot.pl $SEQ AA >> $OUTDIR/job.onehot
printf "done\n"

### Step 4. Translate each encoding to a vector of the given window size.
printf "Translate encodings into vectors ... "
$SCRIPTS/prepare_data_for_learning.pl $OUTDIR/job.aamtx_gaps $WINDOW > $OUTDIR/job.vector_aamtx
$SCRIPTS/prepare_data_for_learning.pl $OUTDIR/job.aaindex $WINDOW > $OUTDIR/job.vector_aaindex
$SCRIPTS/prepare_data_for_learning.pl $OUTDIR/job.onehot $WINDOW > $OUTDIR/job.vector_onehot
printf "done\n"

### Step 5. Merge all features to one vector
printf "Merge vectors ... "
python3 $SCRIPTS/create_vector_features.py -i $OUTDIR/job -o $OUTDIR/job.merged 
printf "done\n"


### Step 6. Make prediction
# Activate Tensorflow conda environment to run the script
# Run all predictions by default: strict, non strict, 3 & 5 classes
printf "Run predictions ... "
$SCRIPTS/medusa.py -i $OUTDIR/job.merged -o $OUTDIR/prediction -f $SEQ -m $MODELS -d 15 100 1>/dev/null 2>/dev/null
printf "done\n\n"

# Create downloadable tar.gz
mkdir $OUTDIR/logs
cp $OUTDIR/job* $OUTDIR/log_hhblits $OUTDIR/logs/
tar -czf $OUTDIR/medusa_job\_results.tar.gz $OUTDIR/prediction/*.csv $OUTDIR/logs 1>/dev/null 2>&1
rm $OUTDIR/job* $OUTDIR/log_hhblits
printf "Results can be found in $ORIGINAL_OUTDIR_NAME\n\n"


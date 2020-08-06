#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "No arguments supplied. Usage: benchmark.sh algorith_name [--dataset id] [--params program_params] [--extrainfo string] [--sequential|--thread numthreads|--mpi numprocs] --time HH:MM:SS"
    exit
fi

declare -a datasets=(GCA_000516895.1_LocustGenomeV1_genomic.fna \
        GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna \
        GCA_900067645.1_Triticum_aestivum_CS42_TGAC_v1_genomic.fna \
        GCF_000001405.26_GRCh38_genomic.fna \
        GCF_000001635.25_GRCm38.p5_genomic.fna \
        GCF_000001895.5_Rnor_6.0_genomic.fna \
        GCF_000002035.5_GRCz10_genomic.fna \
        NA12878.fastq
)


algname=$1
shift

index=1
if [[ "$algname" == "list" ]]; then
        for i in "${datasets[@]}"; do
                size=`du -hs $LUSTRE/TFM/datasets/$i | awk '{ print $1 }'`
                echo "${index}.- $i  ($size)"
                index=$((index + 1))
        done
        exit 0
fi

POSITIONAL=()
PARAMS=""

while [[ $# -gt 0 ]]
do
        key="$1"

        case $key in
            --dataset)
            ID="$2"
            shift # past argument
            shift # past value
            ;;
            --params)
            PARAMS="${PARAMS} $2"
            shift
            shift
            ;;
            --extrainfo)
            EXTRAINFO=$2
            shift
            shift
            ;;
            --sequential)
            SEQUENTIAL="true"
            suffix_outputname="sequential"
            NUMTHREADS=1
            shift
            ;;
            --thread)
            NUMTHREADS=$2
            DOTTHREADS="-t $NUMTHREADS"
            suffix_outputname="threads-$NUMTHREADS"
            THREAD="true"
            shift
            shift
            ;;
            --mpi)
            NUMTHREADS=$2
            MPI="true"
            MPIRUN="--mpi"
            suffix_outputname="threads-$NUMTHREADS"
            shift
            shift
            ;;
            --time)
            TIME="--time $2"
            shift
            shift
            ;;
            *)    # unknown option
            POSITIONAL+=("$1") # save it in an array for later
            shift # past argument
            ;;
        esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters


baseDir="${LUSTRE}/TFM"
resultadosDir="$baseDir/resultados"
datasetsDir="$baseDir/datasets"
executable="$baseDir/src/$algname/dot"
configfile="$baseDir/conf/configDot.cfg"

if [ ! -z $ID  ]; then
        data=${datasets[$((ID-1))]}
else
        data=(${datasets[@]})
fi



if [ "$SEQUENTIAL" != "true" ] && [ "$THREAD" != "true" ] && [ "$MPI" != "true" ]; then
        echo "Must provide a option: --sequential | --thread | --mpi"
        exit 1
fi

for i in "${data[@]}"; do
        for j in `seq 1 3`; do
                OUTPUT="$resultadosDir/$i/$algname/$suffix_outputname/test-$j"
                OUTPUTFILE="$OUTPUT/output"

                mkdir -p $OUTPUT

                cd $OUTPUT
                ln -s $executable dot

                outname="${i}"
                executable_baseopts="-s $datasetsDir/$i -c $configfile"
                executable_opts="$PARAMS ${executable_baseopts} $DOTTHREADS"
                slanzarv --jobname dot-d$ID-$NUMTHREADS-$j $MPIRUN -c $NUMTHREADS $TIME ./dot $executable_opts -o $OUTPUTFILE
        done
done



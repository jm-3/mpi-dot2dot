#!/bin/bash


read -p  "You are going to download and unzip the datasets used in Dot2Dot paper. It will require about 300GB. Do you want to continue? [(a)ll | (c)hoose which ones | (N)o]: " -n 1 -r

if [[ "$REPLY" =~ ^[aAcC]$ ]]; then
        echo
        OPTION=$REPLY

        declare -a datasets=(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz \
                ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/447/015/GCA_001447015.2_Sugar_pine_JHU_assembly/GCA_001447015.2_Sugar_pine_JHU_assembly_genomic.fna.gz \
                ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/067/645/GCA_900067645.1_Triticum_aestivum_CS42_TGAC_v1/GCA_900067645.1_Triticum_aestivum_CS42_TGAC_v1_genomic.fna.gz \
                ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/516/895/GCA_000516895.1_LocustGenomeV1/GCA_000516895.1_LocustGenomeV1_genomic.fna.gz \
                ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.25_GRCm38.p5/GCF_000001635.25_GRCm38.p5_genomic.fna.gz \
                ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.fna.gz \
                ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.5_GRCz10/GCF_000002035.5_GRCz10_genomic.fna.gz \
                http://s3.amazonaws.com/nanopore-human-wgs/rel6/rel_6.fastq.gz
                )

        declare -a sizes=(3.1 27 13 5.6 2.7 2.8 1.3 250) 

        index=0
        for i in "${datasets[@]}"; do
                if [[ "$REPLY" =~ ^[cC]$ ]]; then
                        read -p "Download ${i##*/} (${sizes[index]}GB)? (y/N)" -n 1 -r download
                fi
                if [[ "$download" =~ ^[Yy]$ ]] || [[ "$OPTION" =~ ^[Aa]$ ]]; then
                        echo "Downloading $i ..."
                        wget $i 
                        wait
                        echo "Unzipping ${i##*/*/} ..."
                        gunzip ${i##*/*/} &
                fi
                index=$((index + 1))
        done
        wait
        ls rel_6.fastq &> /dev/null
        if [[ "$?" == 0 ]]; then
                mv rel_6.fastq NA12878.fastq
        fi
fi
echo

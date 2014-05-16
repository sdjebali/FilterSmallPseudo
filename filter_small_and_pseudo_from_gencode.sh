#!/bin/bash

if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: filter_small_and_pseudo_from_gencode.sh annot.gtf >&2
    echo "       This script takes as input a Gencode annotation file in gtf format and outputs this same annotation" >&2
    echo "       without the small transcripts and the pseudogenes, in gtf format, in the directory where it is run." >&2
    echo "       This script assumes that it is located in a directory where there are two sub-directories: Data and Awk." >&2
    echo "       The Data sub-directory must contain a 1 column file with the biotypes of the small transcripts we want to exclude;" >&2
    echo "       The Awk sub-directory must contain two awk scripts: gff2gff.awk and compute_boundaries.awk." >&2
    echo "       This script requires awk to run." >&2
    exit 1
fi


# Initialize variables
######################
echo I am initializing the variables >&2
dir=`dirname $(readlink -m $0)`
annot=$1
basetmp=`basename $annot`
base=${basetmp%.gtf}
echo done >&2


# Filter the annotation
#######################
echo I am removing the objects related to small transcripts and to pseudogenes >&2
awk -v fileRef=$dir/Data/small_tx_types.txt 'BEGIN{while(getline < fileRef >0){small[$1]=1}} (($0!~/pseudogene/)&&($1!~/#/)){found=0; k=9; while((found==0)&&(k<=(NF-1))){if($k=="transcript_type"){found=1; split($(k+1),a,"\"")} k+=2} if(small[a[2]]!=1){print}}' $annot | awk -f $dir/Awk/gff2gff.awk > $base\_filtered_tmp.gtf
echo done >&2

# Compute the new gene boundaries
##################################
echo I am computing the new gene boundaries >&2
awk '$3=="exon"' $base\_filtered_tmp.gtf | sort -k12,12 -k4,4n -k5,5n | awk -v toadd=gene -v fldno=10 -f $dir/Awk/compute_boundaries.awk | awk '{print $10, $1"_"$4"_"$5"_"$7}' > $base\_gnid_coord.txt
echo done >&2

# Adjust gene boundaries
########################
echo I am reporting the new gene boundaries \in the filtered file >&2
awk -v fileRef=$base\_gnid_coord.txt 'BEGIN{while(getline < fileRef >0){coord[$1]=$2; split($2,a,"_"); chr[$1]=a[1]; beg[$1]=a[2]; end[$1]=a[3]; str[$1]=a[4]}} (($3!="gene")&&($1!~/#/)){print} ($3=="gene"){if(coord[$10]==$1"_"$4"_"$5"_"$7){print}else{$1=chr[$10]; $4=beg[$10]; $5=end[$10]; $7=str[$10]; print}}' $base\_filtered_tmp.gtf | awk -f $dir/Awk/gff2gff.awk > $base\_filtered.gtf
echo done >&2

# Clean
#######
echo I am cleaning >&2
rm $base\_filtered_tmp.gtf $base\_gnid_coord.txt 
echo done >&2

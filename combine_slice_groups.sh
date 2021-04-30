#!/bin/bash

if [ $# = 0 ]; then
    echo "Must pass the directories to convert as arguments (* for all, . for current)"
    exit
fi

cur=`pwd`

for filename in "$@"
do
    if [ ! -f "$filename" ]; then
        echo "$filename does not exist"
        continue
    fi

    echo "Re-interleaving $filename"

    nz=`${FSLDIR}/bin/fslval $1 dim3 | awk '{print $1}'`;
    dz=`${FSLDIR}/bin/fslval $1 pixdim3 | awk '{print $1}'`;

    echo "Number of slices is $nz"
    echo "Slice Resolution is $dz"

    mkdir temp~
    fslroi $1 temp~/il1 0 -1 0 -1 0 `expr $nz / 2`
    fslroi $1 temp~/il2 0 -1 0 -1 `expr $nz / 2` `expr $nz / 2`
    fslinterleave temp~/il1 temp~/il2 $1_reinterleaved 
    fslcpgeom temp~/il2 $1_reinterleaved -d
    
    sform=(`${FSLDIR}/bin/fslorient -getsform $1_reinterleaved`)
    echo "Before: ${sform[@]}"
    # echo ${sform[2]}
    sform[2]=`echo "${sform[2]}/2" | bc -l`
    sform[3]=`echo "${sform[3]}-${sform[2]}" | bc -l`
    echo "After:  ${sform[@]}"
    # echo ${sform[2]}
    # echo ${FSLDIR}/bin/fslorient -setsform ${sform[@]} $1_reinterleaved
    ${FSLDIR}/bin/fslorient -setsform ${sform[@]} $1_reinterleaved

    rm -r temp~

done
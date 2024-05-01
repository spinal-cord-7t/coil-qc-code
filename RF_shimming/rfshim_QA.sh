#!/usr/bin/env bash

################################################################################
# This script will generate the RF shim weights text file for the coilQA protocol

# Usage:
# Copy over the coil combined GRE scan and the tfl_rfmap scan to a local drive into DATADIR
# Copy over SarDataUser.mat from C/Medcom/MriProduct/PhysConfig into DATADIR
# Copy over dcm2bids_qa.json into DATADIR

# Input 1: Folder containing the dicoms and other: DATADIR
# Input 2: Output folder. (Should be outside of input 1)

# Example: ./rfshim_QA.sh ~/path/to/dicom/folder ~/path/to/output/folder

# Hard requirement (Shim calculation): Shimming Toolbox: https://shimming-toolbox.org/en/latest
# Hard requirement (Spinal cord masking): SCT: https://spinalcordtoolbox.com/
# Hard requirements (Viewer): FSLeyes: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes
################################################################################

# Create an absolute path of the dicom folder
INPUT_PATH="$(cd "$(dirname "$1")" || exit; pwd)/$(basename "$1")"
# Check if the directory exists
if [ ! -d "${INPUT_PATH}" ]
then
    echo "Input path does not exist"
    exit
fi

# Create the absolute path of the output variable provided
OUTPUT_PATH="$(cd "$(dirname "$2")" || exit; pwd)/$(basename "$2")"
# Check if the output directory exists, if not, create it
if [ ! -d "${OUTPUT_PATH}" ]
then
    echo "Creating output folder"
    mkdir "${OUTPUT_PATH}"
    echo "Creating folders for shim conditions"
    CVRED_PATH="${OUTPUT_PATH}/CVred"
    TARGET20_PATH="${OUTPUT_PATH}/Target20"
    PHASEONLY_PATH="${OUTPUT_PATH}/Phaseonly"
    SAROPT_PATH="${OUTPUT_PATH}/SAROPT"
    mkdir "${CVRED_PATH}"
    mkdir "${TARGET20_PATH}"
    mkdir "${PHASEONLY_PATH}"
    mkdir "${SAROPT_PATH}"
fi


echo "Running DCM2BIDS"
# Call dicom to nifti conversion and create a BIDS structure
st_dicom_to_nifti -i "${INPUT_PATH}" -o "${OUTPUT_PATH}" --subject "CoilQA" --config "${INPUT_PATH}/dcm2bids_QA.json" || exit

# Grab the centerline from the GRE scan
sct_get_centerline -i "${OUTPUT_PATH}/sub-CoilQA/anat/sub-CoilQA_gre.nii.gz" -c "t2s" -o "${OUTPUT_PATH}/sub-CoilQA/anat/centerline.nii.gz"

# Create a mask around this centerline
sct_create_mask -i "${OUTPUT_PATH}/sub-CoilQA/anat/sub-CoilQA_gre.nii.gz" -p centerline,"${OUTPUT_PATH}/sub-CoilQA/anat/centerline.nii.gz" -f cylinder -size 28mm -o "${OUTPUT_PATH}/sub-CoilQA/anat/mask.nii.gz"

#Visualize
fsleyes "${OUTPUT_PATH}/sub-CoilQA/anat/sub-CoilQA_gre.nii.gz" -cm greyscale "${OUTPUT_PATH}/sub-CoilQA/anat/mask.nii.gz" -cm red -a 70.0 &


#Visualize
#fsleyes "${OUTPUT_PATH}/sub-RFshim/anat/sub-RFshim_MP2RAGE_RFSHIM.nii.gz" -cm greyscale "${OUTPUT_PATH}/sub-RFshim/anat/mask.nii.gz" -cm red -a 70.0 &

#Cutting the top and bottom if necessary, and visualising again
echo "ENTER THE BOTTOM AND TOP SLICES OF THE MASK, SUCH THAT IT COVERS FROM C2/C3 TO T2/T3 (this is the Y direction in FSLEYES)!!!"
echo "ENTER BOTTOM SICE (Y direction in FSLEYES)"
read Bottom_sli
echo "BOTTOM SLICE IS: $Bottom_sli"
echo "ENTER TOP SICE"
read Top_sli
echo "TOP SLICE IS: $Top_sli"
sct_crop_image -i "${OUTPUT_PATH}/sub-CoilQA/anat/mask.nii.gz" -ymin $Bottom_sli -ymax $Top_sli -b 0 -o "${OUTPUT_PATH}/sub-CoilQA/anat/mask_crop.nii.gz"
fsleyes "${OUTPUT_PATH}/sub-CoilQA/anat/sub-CoilQA_gre.nii.gz" -cm greyscale "${OUTPUT_PATH}/sub-CoilQA/anat/mask_crop.nii.gz" -cm yellow -a 70.0 &


# TODO: add propt to user if they are happy with the mask
# TODO: figure out how to close FSLEyes instances
#Now that we have the mask, finally we can do RF shimming
st_b1shim --b1 "${OUTPUT_PATH}/sub-CoilQA/rfmap/sub-CoilQA_TB1map_uncombined.nii.gz" --mask "${OUTPUT_PATH}/sub-CoilQA/anat/mask_crop.nii.gz" --vop "${INPUT_PATH}/SarDataUser.mat" --target 20 --sar_factor 1.3 --output "${OUTPUT_PATH}/CVred" --algo 1
st_b1shim --b1 "${OUTPUT_PATH}/sub-CoilQA/rfmap/sub-CoilQA_TB1map_uncombined.nii.gz" --mask "${OUTPUT_PATH}/sub-CoilQA/anat/mask_crop.nii.gz" --vop "${INPUT_PATH}/SarDataUser.mat" --target 20 --sar_factor 1.3 --output "${OUTPUT_PATH}/Target20" --algo 2
st_b1shim --b1 "${OUTPUT_PATH}/sub-CoilQA/rfmap/sub-CoilQA_TB1map_uncombined.nii.gz" --mask "${OUTPUT_PATH}/sub-CoilQA/anat/mask_crop.nii.gz" --vop "${INPUT_PATH}/SarDataUser.mat" --target 20 --sar_factor 1.3 --output "${OUTPUT_PATH}/SAROPT" --algo 3
st_b1shim --b1 "${OUTPUT_PATH}/sub-CoilQA/rfmap/sub-CoilQA_TB1map_uncombined.nii.gz" --mask "${OUTPUT_PATH}/sub-CoilQA/anat/mask_crop.nii.gz" --vop "${INPUT_PATH}/SarDataUser.mat" --target 20 --sar_factor 1.3 --output "${OUTPUT_PATH}/Phaseonly" --algo 4



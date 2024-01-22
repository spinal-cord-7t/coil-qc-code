import os
import shutil
import glob
import json


# This script should be run from the directory containing the input sites directories
project_root = './'

sites = ["MGH", "MNI", "NYU"]

def copy_scan(input_path, output_path):
    nii_path = input_path
    json_path = input_path.replace(".nii.gz", ".json")
    shutil.copy2(nii_path, output_path + ".nii.gz")
    shutil.copy2(json_path, output_path + ".json")

def series_description(json_file_path):
    with open(json_file_path) as f:
        return json.load(f)["SeriesDescription"]


# Output directories are all generated in the outputs folder
output_path_root = os.path.join(project_root, "outputs")
if os.path.exists(output_path_root): shutil.rmtree(output_path_root)
os.makedirs(output_path_root)

for site in sites:
    # Input directories are named SITE-original, for example "MGH-original"
    input_site_path = os.path.join(project_root, site + "-original")

    for subject in [subject_dir for subject_dir in os.listdir(input_site_path) if os.path.isdir(os.path.join(input_site_path, subject_dir))]:
        input_path = os.path.join(input_site_path, subject)
        assert(os.path.exists(input_path))

        subject_prefix = "sub-" + site
        subject = subject.lower()
        if "sub" in subject: subject = subject_prefix + subject[3:]
        else: subject = subject_prefix + subject
        output_path = os.path.join(output_path_root, subject)

        output_anat_path = os.path.join(output_path, "anat")
        output_fmap_path = os.path.join(output_path, "fmap")
        os.makedirs(output_anat_path)
        os.makedirs(output_fmap_path)

        for dir_path in [dir[0] for dir in os.walk(input_path)]:
            dir_basename = os.path.basename(dir_path)

            if dir_basename == "COILQA_SAG_LARGE":
                sag_large_files = sorted(glob.glob(os.path.join(dir_path, "*.nii.gz")))
                snr_input_path = sag_large_files[1 if len(sag_large_files) > 1 else 0]
                copy_scan(snr_input_path, os.path.join(output_fmap_path, subject + "_coilqa-sag-large-SNR"))

            elif dir_basename == "COILQA_SAG_SMALL":
                gfactor_input_path = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))[0]
                copy_scan(gfactor_input_path, os.path.join(output_fmap_path, subject + "_coilqa-sag-small-GFactor"))

            elif dir_basename == "COILQA_TRA":
                gfactor_input_path = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))[0]
                copy_scan(gfactor_input_path, os.path.join(output_fmap_path, subject + "_coilqa-tra-GFactor"))

            elif dir_basename in ["DREAM_MEDIUM", "DREAM_MEDIUM_066", "DREAM_MEDIUM_HWLIMIT"]:
                dream_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                for dream_file_path in dream_file_paths:
                    dream_filename_tokens = os.path.basename(dream_file_path).split("_")
                    eVer = dream_filename_tokens[-1].split(".")[0]
                    series_number = dream_filename_tokens[0]
                    scan_type = dir_basename.lower().replace("_", "-")
                    copy_scan(dream_file_path, os.path.join(output_fmap_path, subject + "_" + scan_type + "-" + eVer + "-" + series_number))

            elif dir_basename == "GRE":
                gre_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                for gre_file_path in gre_file_paths:
                    gre_filename = os.path.basename(gre_file_path)
                    gre_filename_tokens = gre_filename.split("_")
                    if "RX" in gre_filename and not "ph" in gre_filename:
                        channel_name = gre_filename_tokens[5].split(".")[0]
                        copy_scan(gre_file_path, os.path.join(output_anat_path, subject + "_gre-uncombined-" + channel_name + "_T2starw"))

            elif dir_basename == "MP2RAGE":
                mp2rage_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                mp2rage_types = ["INV1", "INV2", "UNI"]
                for i in range(len(mp2rage_file_paths)):
                    series_desc = series_description(mp2rage_file_paths[i].replace(".nii.gz", ".json"))
                    for mp2rage_type_index in range(len(mp2rage_types)):
                        mp2rage_type = mp2rage_types[mp2rage_type_index]
                        if mp2rage_type in series_desc:
                            copy_scan(mp2rage_file_paths[i], os.path.join(output_anat_path, subject + "_mp2rage-" + mp2rage_type))

            elif dir_basename == "TFL_B1_C3C4":
                tfl_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                copy_scan(tfl_file_paths[0], os.path.join(output_fmap_path, subject + "_tfl-mag"))
                copy_scan(tfl_file_paths[2], os.path.join(output_fmap_path, subject + "_tfl-FAmap"))



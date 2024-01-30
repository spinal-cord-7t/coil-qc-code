import os
import shutil
import glob
import json
import re


# This script should be run from the directory containing the input sites directories
project_root = './'

sites = ["MGH", "MNI", "NYU"]
subject_input_names = ["SubD", "SubL", "SubR", "Spinoza6"]

def copy_scan(input_path, output_path, scan_additional=""):
    nii_path = input_path
    json_path = input_path.replace(".nii.gz", ".json")
    output_scan_path = output_path + ".nii.gz"
    if scan_additional != "":
        output_scan_path_tokens = output_scan_path.split("_")
        output_scan_path_tokens.insert(len(output_scan_path_tokens) - 1, json_additional)
        output_scan_path = "_".join(output_scan_path_tokens)
    shutil.copy2(nii_path, output_scan_path)
    shutil.copy2(json_path, output_path + ".json")

def json_attribute(attribute_name, json_file_path):
    with open(json_file_path) as f:
        return json.load(f)[attribute_name]


# Output directories are all generated in the outputs folder
output_path_root = os.path.join(project_root, "outputs")
if os.path.exists(output_path_root): shutil.rmtree(output_path_root)
os.makedirs(output_path_root)

participants_tsv_text = "participant_id\tspecies\tage\tsex\tpathology\tinstitution\tfield\n"

for site in sites:
    # Input directories are named SITE-original, for example "MGH-original"
    input_site_path = os.path.join(project_root, site + "-original")

    for subject in [subject_dir for subject_dir in os.listdir(input_site_path) if os.path.isdir(os.path.join(input_site_path, subject_dir))]:
        input_path = os.path.join(input_site_path, subject)
        assert(os.path.exists(input_path))

        subject_prefix = "sub-" + site
        subject = "sub-" + site + str(subject_input_names.index(subject) + 1)
        output_path = os.path.join(output_path_root, subject)

        participants_tsv_text += "\t".join([subject, "homo sapiens", "n/a", "n/a", "HC", site, "7T"]) + "\n"

        output_anat_path = os.path.join(output_path, "anat")
        output_fmap_path = os.path.join(output_path, "fmap")
        os.makedirs(output_anat_path)
        os.makedirs(output_fmap_path)

        for dir_path in [dir[0] for dir in os.walk(input_path)]:
            dir_basename = os.path.basename(dir_path)

            if dir_basename == "COILQA_SAG_LARGE":
                sag_large_files = sorted(glob.glob(os.path.join(dir_path, "*.nii.gz")))
                snr_input_path = sag_large_files[1 if len(sag_large_files) > 1 else 0]
                copy_scan(snr_input_path, os.path.join(output_fmap_path, subject + "_acq-coilQaSagLarge_SNR"))

            elif dir_basename == "COILQA_SAG_SMALL":
                gfactor_input_path = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))[0]
                copy_scan(gfactor_input_path, os.path.join(output_fmap_path, subject + "_acq-coilQaSagSmall_GFactor"))

            elif dir_basename == "COILQA_TRA":
                gfactor_input_path = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))[0]
                copy_scan(gfactor_input_path, os.path.join(output_fmap_path, subject + "_acq-coilQaTra_GFactor"))

            elif dir_basename in (dream_directory_names := ["DREAM_MEDIUM", "DREAM_MEDIUM_066", "DREAM_MEDIUM_HWLIMIT"]):
                dream_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                for dream_file_path in dream_file_paths:
                    dream_filename_tokens = os.path.basename(dream_file_path).split("_")
                    eVer = dream_filename_tokens[-1].split(".")[0]
                    series_number = dream_filename_tokens[0]
                    scan_type = ["dreamMediumRefV0.66", "dreamMediumRefV1", "dreamMediumRefV1.5"][dream_directory_names.index(dir_basename)]
                    copy_scan(dream_file_path, os.path.join(output_fmap_path, subject + "_acq-" + scan_type + "-" + eVer + "-" + series_number))

            elif dir_basename == "GRE":
                gre_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                for gre_file_path in gre_file_paths:
                    gre_filename = os.path.basename(gre_file_path)
                    gre_filename_tokens = gre_filename.split("_")
                    channel_name = gre_filename_tokens[-1].split(".")[0]
                    if not "ph" in gre_filename:
                        if ("RX" in channel_name or channel_name.startswith("i")):
                            channel_number = re.findall(r"\d+", channel_name)[0]
                            copy_scan(gre_file_path, os.path.join(output_anat_path, subject + "_rec-uncombined" + channel_number + "_T2starw"))
                        elif json_attribute("NonlinearGradientCorrection", gre_file_path.replace(".nii.gz", ".json")) == True:
                            copy_scan(gre_file_path, os.path.join(output_anat_path, subject + "_T2starw"))

            elif dir_basename == "MP2RAGE":
                mp2rage_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                mp2rage_types = ["INV1", "INV2", "UNI"]
                mp2rage_type_names = ["inv-1", "inv-2", "UNIT1"]
                for i in range(len(mp2rage_file_paths)):
                    series_desc = json_attribute("SeriesDescription", mp2rage_file_paths[i].replace(".nii.gz", ".json"))
                    for mp2rage_type_index in range(len(mp2rage_types)):
                        mp2rage_type = mp2rage_types[mp2rage_type_index]
                        if mp2rage_type in series_desc:
                            json_additional = "part-mag"
                            suffix = "_MP2RAGE"
                            if mp2rage_type == "UNI":
                                json_additional = ""
                                suffix = ""
                            mp2rage_type = mp2rage_type_names[mp2rage_type_index]
                            copy_scan(mp2rage_file_paths[i], os.path.join(output_anat_path, subject + "_" + mp2rage_type + suffix), json_additional)

            elif dir_basename == "TFL_B1_C3C4":
                tfl_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                copy_scan(tfl_file_paths[0], os.path.join(output_fmap_path, subject + "_acq-anat_TB1TFL"))
                copy_scan(tfl_file_paths[2], os.path.join(output_fmap_path, subject + "_acq-famp_TB1TFL"))

with open(os.path.join(output_path_root, "participants.tsv"), "w") as f:
    f.write(participants_tsv_text)

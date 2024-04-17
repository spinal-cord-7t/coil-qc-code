import os
import shutil
import glob
import json
import re


# This script should be run from the directory containing the input sites directories
project_root = './'

sites = ["MGH", "MNI", "NYU", "NTNU", "UCL", "Marseille", "MPI"]
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
        json_data = json.load(f)
        if not attribute_name in json_data: return None
        return json_data[attribute_name]


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

        subject_index = subject_input_names.index(subject)
        subject_prefix = "sub-" + site
        subject = "sub-" + site + str(subject_index + 1)
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

            elif dir_basename in (dream_directory_names := ["DREAM_MEDIUM", "DREAM", "DREAM_MEDIUM_066", "DREAM_MEDIUM_HWLIMIT"]):
                dream_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                dream_allowed_scan_types = {
                    "Reference Voltage Map": "refv",
                    "Flipangle Map": "famp",
                    "REFVOLTMAP": "refv",
                    "B1MAP": "famp",
                    "Transmitter Reference Map (Volt)": "refv",
                    "flip angle map": "famp",
                }
                dream_acq_voltage_token = {
                    "DREAM_MEDIUM_066": "0.66",
                    "DREAM_MEDIUM_HWLIMIT": "1.5"
                }
                for dream_file_path in dream_file_paths:
                    dream_json_path = dream_file_path.replace(".nii.gz", ".json")
                    if json_attribute("NonlinearGradientCorrection", dream_json_path): continue
                    scan_type = json_attribute("ImageComments", dream_json_path)
                    if scan_type is None:
                        scan_type = json_attribute("ImageType", dream_json_path)[-1]
                    else:
                        scan_type = scan_type.split(";")[0]
                    if scan_type not in dream_allowed_scan_types and scan_type: continue
                    scan_type = dream_allowed_scan_types[scan_type]
                    voltage_token = ""
                    if dir_basename in dream_acq_voltage_token:
                        voltage_token = "-" + dream_acq_voltage_token[dir_basename]
                    destination_path = os.path.join(output_fmap_path, subject + "_acq-" + scan_type + voltage_token + "_TB1DREAM")
                    
                    def get_fov(json_path):
                        desc_tokens = json_attribute("ProtocolName", json_path).split("_")
                        for token in desc_tokens:
                            if token.startswith("FOV"):
                                return int(token[3:])
                        return None
                    
                    if os.path.exists(destination_path + ".nii.gz"):
                        previous_fov = get_fov(destination_path + ".json")
                        new_fov = get_fov(dream_json_path)
                        if previous_fov == None:
                            print("Could not find FOV token", dream_json_path)
                            continue
                        if new_fov > previous_fov: continue

                    copy_scan(dream_file_path, destination_path)

            elif dir_basename == "GRE":
                gre_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                for gre_file_path in gre_file_paths:
                    gre_filename = os.path.basename(gre_file_path)
                    gre_filename_tokens = gre_filename.split("_")
                    channel_name = gre_filename_tokens[-1].split(".")[0]
                    if not "ph" in gre_filename:
                        if ("RX" in channel_name or channel_name.startswith("i") or channel_name.startswith("c")):
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

            elif dir_basename == "TFL_B1_C3C4" or dir_basename == "TFL" or dir_basename == "TFL_B1_OPT":
                tfl_file_paths = sorted(glob.glob(os.path.join(dir_path, "*nii.gz")))
                anat_found = False
                famp_found = False
                for tfl_file_path in tfl_file_paths:
                    tfl_json_path = tfl_file_path.replace(".nii.gz", ".json")
                    if json_attribute("NonlinearGradientCorrection", tfl_json_path): continue
                    image_comments = json_attribute("ImageComments", tfl_json_path)
                    if image_comments is None: continue
                    if "anatomical image" in image_comments:
                        copy_scan(tfl_file_path, os.path.join(output_fmap_path, subject + "_acq-anat_TB1TFL"))
                        anat_found = True
                    elif "flip angle map" in image_comments:
                        copy_scan(tfl_file_path, os.path.join(output_fmap_path, subject + "_acq-famp_TB1TFL"))
                        famp_found = True
                if not famp_found:
                    copy_scan(tfl_file_paths[[5, 0, 1][subject_index]], os.path.join(output_fmap_path, subject + "_acq-famp_TB1TFL"))

with open(os.path.join(output_path_root, "participants.tsv"), "w") as f:
    f.write(participants_tsv_text)

shutil.copy2(os.path.join(project_root, ".bidsignore"), os.path.join(output_path_root, ".bidsignore"))
shutil.copy2(os.path.join(project_root, "dataset_description.json"), os.path.join(output_path_root, "dataset_description.json"))

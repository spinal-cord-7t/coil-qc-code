#!/usr/bin/python3
# -*- coding: utf-8 -*-

# to run from the command line: 
# python3 dream_FAmap.py Iste.nii Ifid.nii

import numpy as np
import os
import nibabel as nib
import sys

Iste_fname = sys.argv[1]
Ifid_fname = sys.argv[2]

Iste_img = nib.load(Iste_fname)
Ifid_img = nib.load(Ifid_fname)

Iste_data = Iste_img.get_fdata()
Ifid_data = Ifid_img.get_fdata()


FA = np.array(np.zeros_like(Iste_data))
FA = np.arctan(np.sqrt(np.divide(2*Iste_data,Ifid_data, where=Ifid_data!=0)))

nii_FA = nib.Nifti1Image(FA, Iste_img.affine)
nib.save(nii_FA, 'dream_FA.nii.gz')

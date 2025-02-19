#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last updated: 27/1 2025

@author: Joachim Larsen
"""

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import pytesseract as pt
import numpy as np
import pydicom as pyd
from tkinter.filedialog import askopenfilename
import functools
import os
import pandas as pd
from glob import glob

save_path = '/home/jlar0426/Documents/csv/lfu.csv'

def read_dcm(path = None):
    if path == None:
        paths = askopenfilename()
    else:
        paths = path
    dataset = pyd.dcmread(paths)
    return dataset, paths

def get_paths(dicom, path):
    
    exam_paths = []
    exam_data = []
    dirloc = path.split('DICOMDIR')[0]
    for i, series in enumerate(dicom[0x00041220]):
        if i == 0:
            pid = series[0x00100020].value
        if (0x00080020) in series:
            date = series[0x00080020].value
            act_date = date[-2:] +'-' + date[-4:-2] + '-' + date[:4]
        if (0x00041510) in series:
            name = series[0x00041510].repval
            
            if name == 'Secondary Capture Image Storage':
                path = series[0x00041500]
                p = functools.reduce(lambda a,b : a + '/' + b, path)
                total_path = dirloc + p + '.JPG'
                total_path = total_path.replace('DICOM', 'IHE_PDI')
# Currently for testing only returns the first found
                idd = {'ID' : pid, 'Date' : act_date}
                exam_paths.append(total_path)
                exam_data.append(idd)
#                return total_path, idd
#                print(total_path)
    return exam_paths, exam_data


def read_img_data(test_path, idd):


#    test_path, idd = get_paths(dataset, path)


#test_path = '/home/jlar0426/Documents/LFU_data/LFU/0103460778/IHE_PDI/0000806C/AA919865/AA06D7C0/0000D659/EE7899CC.JPG'

    img = plt.imread(test_path)

#plt.imshow(img)
#plt.show()

# x range: 700 - 1300, y range: 200 - 470

    imgtest1 = pt.image_to_string(img[200:470, 1020:1300])
    data_info = imgtest1.split('\n')
    data_info = [s for s in data_info if s.strip()]

    data_info[1] = data_info[1].split(' ')[0]
    data_info[2] = data_info[2].split(' ')[0]
    data_info[3] = data_info[3].split(' ')[0]

    info_parts = ['gender', 'height', 'weight', 'age']

    info_dict = dict(zip(info_parts, data_info))

# CPR: y range : 215 - 270, x range: 360 - 630
#    imgtest2 = pt.image_to_string(img[215:270, 360:630])

#    cpr = imgtest2.replace(' ','').replace('\n','')


# Data: x range : 600 - 800, y range: 670 - 1350
# Final range: y 670 - 1390. Chosen very specifically since it is able to read
# properly both with and without the time stamp existing in a dataset.
    imgtest3 = pt.image_to_string(img[670:1390, 600:800], config = '--psm 6 -c tessedit_char_whitelist=0123456789.')

    data_vals = imgtest3.split('\n')
    data_vals = [s for s in data_vals if s.strip()]
    if len(data_vals) > 18:
        data_vals = data_vals[1:]
    elif len(data_vals) < 18:
        o_vals = [0] * (18 - len(data_vals))
        data_vals = data_vals + o_vals
    
# Mini tests
#imgtest4 = pt.image_to_string(img[680:890, 600:800], config = '--psm 6 -c tessedit_char_whitelist=0123456789.')
#imgtest5 = pt.image_to_string(img[910:1120, 600:800], config = '--psm 6 -c tessedit_char_whitelist=0123456789.')
#imgtest6 = pt.image_to_string(img[1140:1350, 600:800], config = '--psm 6 -c tessedit_char_whitelist=0123456789.')

# The whole data line is the same as creating 3 areas excluding the line between.

    parts = ['FEV 1', 'FVC', 'FEV 1 % FVC', 'FEV 1 % VC MAX', 'PEF', 'MEEF 75/25',
             'VC IN', 'ITGV', 'ERV', 'RV', 'TLC', 'RV % TLC',
             'DLCO SB', 'DLCO/VA', 'DLCOc SB', 'VIN', 'TLC-SB', 'Hb']

    data_dict = dict(zip(parts, data_vals))

#cpr_dict = {'ID' : cpr}

    info_dict.update(data_dict)
#info_dict.update(cpr_dict)
    info_dict.update(idd)

    for k in info_dict:
        info_dict[k] = [str(info_dict[k]).replace(',', '.')]

    df = pd.DataFrame.from_dict(info_dict)
    
# What needs to be given a new button
    if os.path.isfile(save_path):
        pid = idd['ID']
        date = idd['Date']
        ori_df = pd.read_csv(save_path, dtype = 'O')
        del_df = ori_df[(ori_df['ID'] == pid) & (ori_df['Date'] == date)]

        gen_check = (df['gender'][0] == 'female') or (df['gender'][0] == 'male')
        if (not gen_check) & (not del_df.empty):
            print('weird things in this image, skips')
            return
        if not del_df.empty:
            ori_df = ori_df.drop(del_df.index.values)

                
        new_df = pd.concat([ori_df, df])
        
        
        new_df.to_csv(save_path, mode='w', index=False)
    else:
        df.to_csv(save_path, mode='a', index=False)        
    return

def mapped(name):
    data, path = read_dcm(name)
    dicom, ex_data = get_paths(data, path)
    for i in range(len(dicom)):
        print(dicom[i], ex_data[i])
        read_img_data(dicom[i], ex_data[i])

if __name__ == '__main__':
    
    files = glob('/home/jlar0426/Documents/LFU_data/LFU/*')
    names = list(map(lambda x : os.path.basename(x) + '/DICOMDIR', files))
    directory = os.path.dirname(files[0])
    os.chdir(directory)

    list(map(lambda x : mapped(x), names))

    print('finished')
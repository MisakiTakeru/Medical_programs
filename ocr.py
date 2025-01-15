#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 09:11:04 2025

@author: jlar0426
"""

import matplotlib.pyplot as plt
import pytesseract as pt
import numpy as np
import pydicom as pyd
from tkinter.filedialog import askopenfilename
import functools
import os
import pandas as pd
import tkinter as tk
from tkinter import ttk
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
                return total_path, idd
                print(total_path)


def read_img_data(dataset, path):


    test_path, idd = get_paths(dataset, path)


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
    imgtest2 = pt.image_to_string(img[215:270, 360:630])

    cpr = imgtest2.replace(' ','').replace('\n','')


# Data: x range : 600 - 800, y range: 670 - 1350
# Final range: y 670 - 1390. Chosen very specifically since it is able to read
# properly both with and without the time stamp existing in a dataset.
    imgtest3 = pt.image_to_string(img[670:1390, 600:800], config = '--psm 6 -c tessedit_char_whitelist=0123456789.')

    data_vals = imgtest3.split('\n')
    data_vals = [s for s in data_vals if s.strip()]
    if len(data_vals) > 18:
        data_vals = data_vals[1:]

# Mini tests
#imgtest4 = pt.image_to_string(img[680:890, 600:800], config = '--psm 6 -c tessedit_char_whitelist=0123456789.')
#imgtest5 = pt.image_to_string(img[910:1120, 600:800], config = '--psm 6 -c tessedit_char_whitelist=0123456789.')
#imgtest6 = pt.image_to_string(img[1140:1350, 600:800], config = '--psm 6 -c tessedit_char_whitelist=0123456789.')

# The whole data line is the same as creating 3 areas excluding the line between.

    parts = ['FEV 1', 'FVC', 'FEV 1 % FVC', 'FEV 1 % VC MAX', 'PEF', 'MEEF 75/25',
             'VC IN', 'ITGV', 'ERV', 'RV', 'TLC', 'RV % TLC',
             'DLSO SB', 'DLCO/VA', 'DLCOc SB', 'VIN', 'TLC-SB', 'Hb']

    data_dict = dict(zip(parts, data_vals))

#cpr_dict = {'ID' : cpr}

    info_dict.update(data_dict)
#info_dict.update(cpr_dict)
    info_dict.update(idd)

    for k in info_dict:
        info_dict[k] = [info_dict[k]]

    df = pd.DataFrame.from_dict(info_dict)

    if os.path.isfile(save_path):
        pid = info_dict['ID'][0]
        date = info_dict['Date'][0]
        ori_df = pd.read_csv(save_path)
                                    
        new_df = ori_df[(ori_df['ID'] != pid) & (ori_df['Date'] != date)]
            
        new_df = pd.concat([new_df, df])
            
        new_df.to_csv(save_path, mode='w', index=False)
    else:
        df.to_csv(save_path, mode='a', index=False)

class GUIApp(tk.Frame):
    def __init__(self, parent, files, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.files = files
        
        self.option_value = tk.StringVar(self)

        self.drop_values = files

        self.dropbox = ttk.Combobox(
            parent,
            textvariable =self.option_value,
            values = self.drop_values,
            state = 'readonly')
                

        self.dropbox.place(relx = 0.1, rely = 0.05, relwidth= 0.15)

if __name__ == '__main__':
#    dataset, path = read_dcm()
#    read_img_data(dataset, path)
    window = tk.Tk()

    files = glob('/home/jlar0426/Documents/LFU_data/LFU/*')
    names = list(map(lambda x : os.path.basename(x), files))
    directory = os.path.dirname(files[0])
    os.chdir(directory)
    
    gui = GUIApp(window, names)

    window.minsize(1200,600)
    
# read_dcm(names[0] + '/DICOMDIR')

    window.mainloop()
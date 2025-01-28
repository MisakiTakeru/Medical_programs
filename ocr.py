#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 09:11:04 2025

@author: jlar0426
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
import tkinter as tk
from tkinter import ttk
from glob import glob
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

load_path = '/home/jlar0426/Documents/csv/lfu.csv'
save_path = '/home/jlar0426/Documents/csv/lfu_val.csv'

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

        # binding dropbox to the function update, when it's value is changed.
        self.dropbox.bind("<<ComboboxSelected>>", self.update)                

        self.dropbox.place(relx = 0.1, rely = 0.05, relwidth= 0.15)

        self.drop_values_2 = []
        
        self.option_value_2 = tk.StringVar(self)
        
        self.dropbox_2 = ttk.Combobox(
            parent,
            textvariable = self.option_value_2,
            values = self.drop_values_2,
            state = 'readonly')
        
        self.dropbox_2.place(relx = 0.3, rely = 0.05, relwidth = 0.15)

        self.button = tk.Button(
            self.parent,
            text = 'go through chosen examination',
            command = lambda : self.read_img())
        self.button.place(relx = 0.5, rely = 0.05, relwidth = 0.15)

        self.button_2 = tk.Button(
            self.parent,
            text = 'Save data',
            command = lambda : self.save_data())
        self.button_2.place(relx = 0.7, rely = 0.05, relwidth = 0.15)

# Long list of all labels and updatable entries to manually change.        
        self.gender_var = tk.StringVar()
        self.gender_label = tk.Label(self.parent, text = 'gender:')
        self.gender_label.place(relx = 0.52, rely = 0.2)
        self.gender = tk.Entry(self.parent, textvariable = self.gender_var)
        self.gender.place(relx = 0.60, rely = 0.2)

        self.height_var = tk.StringVar()
        self.height_label = tk.Label(self.parent, text = 'height:')
        self.height_label.place(relx = 0.52, rely = 0.25)
        self.height = tk.Entry(self.parent, textvariable = self.height_var)
        self.height.place(relx = 0.60, rely = 0.25)
        
        self.weight_var = tk.StringVar()
        self.weight_label = tk.Label(self.parent, text = 'weight:')
        self.weight_label.place(relx = 0.52, rely = 0.30)
        self.weight = tk.Entry(self.parent, textvariable = self.weight_var)
        self.weight.place(relx = 0.60, rely = 0.3)
        
        self.age_var = tk.StringVar()
        self.age_label = tk.Label(self.parent, text = 'age:')
        self.age_label.place(relx = 0.52, rely = 0.35)
        self.age = tk.Entry(self.parent, textvariable = self.age_var)
        self.age.place(relx = 0.60, rely = 0.35)

        self.fev1_var = tk.StringVar()
        self.fev1_label = tk.Label(self.parent, text = 'FEV 1:')
        self.fev1_label.place(relx = 0.52, rely = 0.40)
        self.fev1 = tk.Entry(self.parent, textvariable = self.fev1_var)
        self.fev1.place(relx = 0.60, rely = 0.40)

        self.fvc_var = tk.StringVar()
        self.fvc_label = tk.Label(self.parent, text = 'FVC:')
        self.fvc_label.place(relx = 0.52, rely = 0.45)
        self.fvc = tk.Entry(self.parent, textvariable = self.fvc_var)
        self.fvc.place(relx = 0.60, rely = 0.45)

        self.fevfvc_var = tk.StringVar()
        self.fevfvc_label = tk.Label(self.parent, text = 'FEV 1 % FVC:')
        self.fevfvc_label.place(relx = 0.52, rely = 0.50)
        self.fevfvc = tk.Entry(self.parent, textvariable = self.fevfvc_var)
        self.fevfvc.place(relx = 0.60, rely = 0.50)

        self.fevvcm_var = tk.StringVar()
        self.fevvcm_label = tk.Label(self.parent, text = 'FEV 1 % VC MAX:')
        self.fevvcm_label.place(relx = 0.52, rely = 0.55)
        self.fevvcm = tk.Entry(self.parent, textvariable = self.fevvcm_var)
        self.fevvcm.place(relx = 0.60, rely = 0.55)

        self.pef_var = tk.StringVar()
        self.pef_label = tk.Label(self.parent, text = 'PEF:')
        self.pef_label.place(relx = 0.52, rely = 0.60)
        self.pef = tk.Entry(self.parent, textvariable = self.pef_var)
        self.pef.place(relx = 0.60, rely = 0.60)

        self.meef_var = tk.StringVar()
        self.meef_label = tk.Label(self.parent, text = 'MEEF 75/25:')
        self.meef_label.place(relx = 0.52, rely = 0.65)
        self.meef = tk.Entry(self.parent, textvariable = self.meef_var)
        self.meef.place(relx = 0.60, rely = 0.65)

        self.vcin_var = tk.StringVar()
        self.vcin_label = tk.Label(self.parent, text = 'VC IN:')
        self.vcin_label.place(relx = 0.52, rely = 0.70)
        self.vcin = tk.Entry(self.parent, textvariable = self.vcin_var)
        self.vcin.place(relx = 0.60, rely = 0.70)
        
        self.itgv_var = tk.StringVar()
        self.itgv_label = tk.Label(self.parent, text = 'ITGV:')
        self.itgv_label.place(relx = 0.52, rely = 0.75)
        self.itgv = tk.Entry(self.parent, textvariable = self.itgv_var)
        self.itgv.place(relx = 0.60, rely = 0.75)
        
        self.erv_var = tk.StringVar()
        self.erv_label = tk.Label(self.parent, text = 'ERV:')
        self.erv_label.place(relx = 0.52, rely = 0.80)
        self.erv = tk.Entry(self.parent, textvariable = self.erv_var)
        self.erv.place(relx = 0.60, rely = 0.80)

        self.rv_var = tk.StringVar()
        self.rv_label = tk.Label(self.parent, text = 'RV:')
        self.rv_label.place(relx = 0.52, rely = 0.85)
        self.rv = tk.Entry(self.parent, textvariable = self.rv_var)
        self.rv.place(relx = 0.60, rely = 0.85)

        self.tlc_var = tk.StringVar()
        self.tlc_label = tk.Label(self.parent, text = 'TLC:')
        self.tlc_label.place(relx = 0.52, rely = 0.90)
        self.tlc = tk.Entry(self.parent, textvariable = self.tlc_var)
        self.tlc.place(relx = 0.60, rely = 0.90)

        self.rvtlc_var = tk.StringVar()
        self.rvtlc_label = tk.Label(self.parent, text = 'RV % TLC:')
        self.rvtlc_label.place(relx = 0.52, rely = 0.95)
        self.rvtlc = tk.Entry(self.parent, textvariable = self.rvtlc_var)
        self.rvtlc.place(relx = 0.60, rely = 0.95)

        self.dlco_var = tk.StringVar()
        self.dlco_label = tk.Label(self.parent, text = 'DLCO SB:')
        self.dlco_label.place(relx = 0.75, rely = 0.20)
        self.dlco = tk.Entry(self.parent, textvariable = self.dlco_var)
        self.dlco.place(relx = 0.83, rely = 0.20)

        self.dlcova_var = tk.StringVar()
        self.dlcova_label = tk.Label(self.parent, text = 'DLCO/VA:')
        self.dlcova_label.place(relx = 0.75, rely = 0.25)
        self.dlcova = tk.Entry(self.parent, textvariable = self.dlcova_var)
        self.dlcova.place(relx = 0.83, rely = 0.25)

        self.dlcosb_var = tk.StringVar()
        self.dlcosb_label = tk.Label(self.parent, text = 'DLCOc SB:')
        self.dlcosb_label.place(relx = 0.75, rely = 0.30)
        self.dlcosb = tk.Entry(self.parent, textvariable = self.dlcosb_var)
        self.dlcosb.place(relx = 0.83, rely = 0.30)

        self.vin_var = tk.StringVar()
        self.vin_label = tk.Label(self.parent, text = 'VIN:')
        self.vin_label.place(relx = 0.75, rely = 0.35)
        self.vin = tk.Entry(self.parent, textvariable = self.vin_var)
        self.vin.place(relx = 0.83, rely = 0.35)

        self.tlcsb_var = tk.StringVar()
        self.tlcsb_label = tk.Label(self.parent, text = 'TLC-SB:')
        self.tlcsb_label.place(relx = 0.75, rely = 0.40)
        self.tlcsb = tk.Entry(self.parent, textvariable = self.tlcsb_var)
        self.tlcsb.place(relx = 0.83, rely = 0.40)

        self.hb_var = tk.StringVar()
        self.hb_label = tk.Label(self.parent, text = 'Hb:')
        self.hb_label.place(relx = 0.75, rely = 0.45)
        self.hb = tk.Entry(self.parent, textvariable = self.hb_var)
        self.hb.place(relx = 0.83, rely = 0.45)

        self.f = plt.Figure()
        self.f.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        self.a = self.f.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.f, self.parent)
        self.canvas.get_tk_widget().place(relx = 0.1, rely = 0.10, relwidth = 0.4, relheight = 0.80)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self.parent)
        self.toolbar.update()

        
    def update(self, event):
        dataset, path = read_dcm(self.option_value.get() + '/DICOMDIR')
        exams, ex_data = get_paths(dataset, path)
        self.exams = list(zip(exams, ex_data))
        exam_dates = [x['Date'] for x in ex_data]
        self.drop_values_2 = exam_dates
        self.dropbox_2.config(values = self.drop_values_2)
        self.exams_dates = exam_dates
        self.dropbox_2.set('')

    def read_img(self):
        val = self.dropbox_2.get()
        print(val)
        if val == '':
            return
        else:
            index = [i for i,x in enumerate(self.exams_dates) if val == x][0]
            data = self.exams[index]

            data_generator = create_gens(load_path)
            
            ids = self.dropbox.get()
            
            df_gen = (
                d
                for d in data_generator
                if (val in d['Date']) & (ids in d['ID'])
                )
    

            for d in df_gen:
                df = d

            print(ids)
            print(val)
            print('df')
            print(df)
# set all of the different variables parameters (next time I should probably create them as a dictionary to loop through).
            self.gender_var.set(df['gender'])
            self.height_var.set(df['height'])
            self.weight_var.set(df['weight'])
            self.age_var.set(df['age'])
            self.fev1_var.set(df['FEV 1'])
            self.fvc_var.set(df['FVC'])
            self.fevfvc_var.set(df['FEV 1 % FVC'])
            self.fevvcm_var.set(df['FEV 1 % VC MAX'])
            self.pef_var.set(df['PEF'])
            self.meef_var.set(df['MEEF 75/25'])
            self.vcin_var.set(df['VC IN'])
            self.itgv_var.set(df['ITGV'])
            self.erv_var.set(df['ERV'])
            self.rv_var.set(df['RV'])
            self.tlc_var.set(df['TLC'])
            self.rvtlc_var.set(df['RV % TLC'])
            self.dlco_var.set(df['DLCO SB'])
            self.dlcova_var.set(df['DLCO/VA'])
            self.dlcosb_var.set(df['DLCOc SB'])
            self.vin_var.set(df['VIN'])
            self.tlcsb_var.set(df['TLC-SB'])
            self.hb_var.set(df['Hb'])
            self.id = df['ID']
            self.date = df['Date']
            
# read image, clear it and show the new image and draw it
            img = plt.imread(data[0])            
            self.a.cla()
            self.a.imshow(img)
            self.a.axis('off')

            self.canvas.draw()

    def save_data(self):

        keys = ['gender', 'height', 'weight', 'age', 'FEV 1', 'FVC',
                 'FEV 1 % FVC', 'FEV 1 % VC MAX', 'PEF', 'MEEF 75/25',
                 'VC IN', 'ITGV', 'ERV', 'RV', 'TLC', 'RV % TLC',
                 'DLCO SB', 'DLCO/VA', 'DLCOc SB', 'VIN', 'TLC-SB', 'Hb',
                 'ID', 'Date']

        vals = [self.gender.get(), self.height.get(), self.weight.get(),
                self.age.get(), self.fev1.get(), self.fvc.get(),
                self.fevfvc.get(), self.fevvcm.get(), self.pef.get(),
                self.meef.get(), self.vcin.get(), self.itgv.get(),
                self.erv.get(), self.rv.get(), self.tlc.get(),
                self.rvtlc.get(), self.dlco.get(), self.dlcova.get(),
                self.dlcosb.get(), self.vin.get(), self.tlcsb.get(),
                self.hb.get(), self.id, self.date]
        
        d = dict(map(lambda k,v : (k,[v]), keys, vals))
        print(d)
        df = pd.DataFrame.from_dict(d)

# What needs to be given a new button
        if os.path.isfile(save_path):
            pid = self.id
            date = self.date
            ori_df = pd.read_csv(save_path, dtype = 'O')
                                        
            new_df = ori_df[(ori_df['ID'] != pid) & (ori_df['Date'] != date)]
                
            new_df = pd.concat([new_df, df])
            
            new_df.to_csv(save_path, mode='w', index=False)
        else:
            df.to_csv(save_path, mode='a', index=False)        
        return


# Function that creates a small pipeline of generators for reading single lines
# and turn them into a dictionary og columns and the lines.
def create_gens(file_name):

# Generator for running through the lines one by one
    lines = (line for line in open(file_name, encoding='utf8'))
# Generator for tab splitting the lines
    list_lines = (line.strip().split(',') for line in lines)
# Gets the first line to use as dictionary keywords
    cols = next(list_lines)
# Generator that creates a dictionary of columns and the splittet data
    dicts = (dict(zip(cols, data)) for data in list_lines)
    return dicts



if __name__ == '__main__':
    window = tk.Tk()

    files = glob('/home/jlar0426/Documents/LFU_data/LFU/*')
    names = list(map(lambda x : os.path.basename(x), files))
    directory = os.path.dirname(files[0])
    os.chdir(directory)
    
    gui = GUIApp(window, names)

    window.minsize(1200,600)
    
# read_dcm(names[0] + '/DICOMDIR')

    window.mainloop()
    

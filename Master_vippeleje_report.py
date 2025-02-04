#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Last updated: 30/1 2025

@author: Joachim Larsen
"""


import tkinter as tk
import sin_fit
import pandas as pd
from tkinter import ttk
import choose_data
import os
import numpy as np

read_path = '/home/jlar0426/Documents/csv'

# Location the data is saved to.
save_path = '/home/jlar0426/Documents/csv/vip_fin.csv'

# Location from which database data is read from
read_from = read_path + '/vipexam.csv'

# The column order of which the data is saved in.
columns = ['resperation_hr_mean', 'resperation_hr_std', 'resperation_hr_maxs',
           'resperation_hr_mins', 'resperation_hr_chosen', 'resperation_hr_fit_mean',
           'resperation_hr_fit_std', 'resperation_hr_fit_maxs',
           'resperation_hr_fit_mins', 'resperation_sBP_mean',
           'resperation_sBP_std', 'resperation_sBP_maxs', 'resperation_sBP_mins',
           'resperation_sBP_chosen', 'resperation_sBP_fit_mean', 'resperation_sBP_fit_std',
           'resperation_sBP_fit_maxs', 'resperation_sBP_fit_mins',
           'resperation_dBP_mean', 'resperation_dBP_std', 'resperation_dBP_maxs',
           'resperation_dBP_mins', 'resperation_dBP_chosen', 'resperation_dBP_fit_mean',
           'resperation_dBP_fit_std', 'resperation_dBP_fit_maxs',
           'resperation_dBP_fit_mins', 'Active standing_hr_mean',
           'Active standing_hr_std', 'Active standing_sBP_mean',
           'Active standing_sBP_std', 'Active standing_dBP_mean',
           'Active standing_dBP_std', 'Active standing_intervals',
           'Active standing_sBP_min', 'Active standing_dBP_min',
           'Active standing_HR_max', 'Tilt_hr_mean', 'Tilt_hr_std',
           'Tilt_sBP_mean', 'Tilt_sBP_std', 'Tilt_dBP_mean', 'Tilt_dBP_std',
           'Tilt_intervals']

# resp [0 : 27]
# acitve [27 : 37]
# tilt [37 : 44]
# fit version is same but without chosen
# [mean_hr, std_hr, maxs_hr, mins_hr, chosen]
class Vip_data():
    def __init__(self):
        self.resp_data = {}
        self.active_data = {}
        self.tilt_data = {}
        self.carotis_data = {}
        self.id = None
        self.date = None

    def reset(self):
        self.resp_data = {}
        self.id = None
        self.date = None
        self.active_data = {}
        self.tilt_data = {}
        self.carotis_data = {}
    
    def update_choose_data(self,
        sBP_mean, sBP_std,
        dBP_mean, dBP_std,
        hr_mean, hr_std,
        areas, action, minmaxs = {}):
        
        dicts = {action + '_hr_mean' : hr_mean, action + '_hr_std' : hr_std,
                 action + '_sBP_mean' : sBP_mean, action + '_sBP_std' : sBP_std,
                 action + '_dBP_mean' : dBP_mean, action + '_dBP_std' : dBP_std,
                 action + '_intervals' : areas}
        
        if not minmaxs == dict():
            for k in minmaxs.keys():
                minmaxs[action + '_' + k] = minmaxs.pop(k)
            dicts.update(minmaxs)
        if action == 'Active standing':
            self.active_data = dicts
        elif action == 'Tilt':
            self.tilt_data = dicts
        elif action == 'Carotis':
            self.carotis_data = dicts
    
    def update_resp_data(self, hr, hr_fit, sBP, sBP_fit, dBP, dBP_fit):
        
# Calls helper function to create a dictionary of the hr data.
        hr_dict = self.resp_helper(hr, 'hr')
    
        hr_fit_dict = self.resp_helper(hr_fit, 'hr_fit', True)
# Merges the two dictionaries.
        hr_dict.update(hr_fit_dict)
        
        sBP_dict = self.resp_helper(sBP, 'sBP', False)

        hr_dict.update(sBP_dict)

        sBP_fit_dict = self.resp_helper(sBP_fit, 'sBP_fit', True)

        hr_dict.update(sBP_fit_dict)

        dBP_dict = self.resp_helper(dBP, 'dBP', False)

        hr_dict.update(dBP_dict)

        dBP_fit_dict = self.resp_helper(dBP_fit, 'dBP_fit', True)

        hr_dict.update(dBP_fit_dict)
        
        self.resp_data = hr_dict
        
    
    def resp_helper(self, data, name, fit = False):
        """
        Parameters
        ----------
        data : list of lists
            list of lists containing the data from sinudal fit.
        name : string
            The name of the type of data.
        fit : Boolean, optional
            Whether it is the fitted data, since the fit data does not have 
            the chosen parameter. The default is False.

        Returns
        -------
        res : dict
            dictionary with all of the data.

        """
        mean = []; std = []; maxs = []; mins = []
        if not fit:
            chosen = []
        for i in range(len(data)):
            
            mean.append(data[i][0])
            std.append(data[i][1])
            maxs.append(data[i][2])
            mins.append(data[i][3])
            
            if not fit:
                chosen.append(data[i][4])

        res = {'resperation_' + name + '_mean' : mean,
               'resperation_' + name + '_std' : std,
               'resperation_' + name + '_maxs' : maxs,
               'resperation_' + name + '_mins' : mins}

        if not fit:
            ch = {'resperation_' + name + '_chosen': chosen}
            res.update(ch)
        
        return res

    def to_dataframe(self):
        df = {}
        resp_df = pd.DataFrame(self.resp_data, columns = columns[:27])
        df.update(self.active_data)
        df.update(self.tilt_data)
        df.update(self.carotis_data)
        
        df_rest = pd.DataFrame([df], columns = columns[27:])
        return pd.concat([resp_df, df_rest], axis = 1)


class GUIApp(tk.Frame):
    def __init__(self, parent, database, df, *args, **kwargs):
        tk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.database = database
                
        # Initializing buttons and dropdown menu.
        self.button = tk.Button(self.parent, text = 'sinudal fit', command = lambda : self.run_fun())
        self.button.place(relx = 0.3, rely = 0.05)

        self.button_active = tk.Button(
            self.parent,
            text = 'Active Standing',
            command = lambda : self.run_choose('Active standing'))
        self.button_active.place(relx = 0.5, rely = 0.05)

        self.button_tilt = tk.Button(
            self.parent,
            text = 'Tilt',
            command = lambda : self.run_choose('Tilt'))
        self.button_tilt.place(relx = 0.7, rely = 0.05)
        
        self.button_save = tk.Button(
            self.parent,
            text = 'Save data',
            command = lambda : self.save())
        self.button_save.place(relx = 0.9, rely = 0.05)
        
        
        self.option_value = tk.StringVar(self)

        self.drop_values = list(np.unique(df.apply(lambda row : str(row.iloc[0]) + ',' + str(row.iloc[1]), axis = 1)))

        self.dropbox = ttk.Combobox(
            parent,
            textvariable =self.option_value,
            values = self.drop_values,
            state = 'readonly')
        
        # binding dropbox to the function update, when it's value is changed.
        self.dropbox.bind("<<ComboboxSelected>>", self.update)
        

        self.dropbox.place(relx = 0.1, rely = 0.05, relwidth= 0.15)
        
        self.drop_label = tk.Label(self.parent, text = 'Choose examination data found:')
        self.drop_label.place(relx = 0.1, rely = 0.025)
        
        self.sinudal_fit = None
    
        self.label_sinudal = tk.Label(self.parent, text = 'sinudal data: ')
        self.label_active = tk.Label(self.parent, text = 'Active standing data:')
        self.label_tilt = tk.Label(self.parent, text = 'Tilt test data:')
        
        self.label_sinudal.place(relx = 0.01, rely = 0.2)
        self.label_active.place(relx = 0.31, rely = 0.2)
        self.label_tilt.place(relx = 0.61, rely = 0.2)

# Sinudal fit labels
        self.label_hr_title = tk.Label(self.parent, text = 'HR:')
        self.label_hr_title.place(relx = 0.01, rely = 0.25)

        self.label_hr_mean = tk.Label(self.parent, text = 'Means: ')
        self.label_hr_mean.place(relx = 0.01, rely = 0.28)
        self.label_hr = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.label_hr.place(relx = 0.05, rely = 0.27)

        self.label_hr_std = tk.Label(self.parent, text = 'Stds: ')
        self.label_hr_std.place(relx = 0.01, rely = 0.33)
        self.label_hr1 = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.label_hr1.place(relx = 0.05, rely = 0.32)
        

        self.label_dBP_title = tk.Label(self.parent, text = 'dBP:')
        self.label_dBP_title.place(relx = 0.01, rely = 0.45)

        self.label_dBP_mean = tk.Label(self.parent, text = 'Means: ')
        self.label_dBP_mean.place(relx = 0.01, rely = 0.50)
        self.label_dBP = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.label_dBP.place(relx = 0.05, rely = 0.49)

        self.label_dBP_std = tk.Label(self.parent, text = 'Stds: ')
        self.label_dBP_std.place(relx = 0.01, rely = 0.55)
        self.label_dBP1 = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.label_dBP1.place(relx = 0.05, rely = 0.54)
        

        self.label_sBP_title = tk.Label(self.parent, text = 'sBP:')
        self.label_sBP_title.place(relx = 0.01, rely = 0.65)

        self.label_sBP_mean = tk.Label(self.parent, text = 'Means: ')
        self.label_sBP_mean.place(relx = 0.01, rely = 0.70)
        self.label_sBP = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.label_sBP.place(relx = 0.05, rely = 0.69)

        self.label_sBP_std = tk.Label(self.parent, text = 'Stds: ')
        self.label_sBP_std.place(relx = 0.01, rely = 0.75)
        self.label_sBP1 = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.label_sBP1.place(relx = 0.05, rely = 0.74)


# Active stand labels
        
        self.active_hr_title = tk.Label(self.parent, text = 'HR')
        self.active_hr_title.place(relx = 0.31, rely = 0.25)

        self.active_hr_mean = tk.Label(self.parent, text = 'Means: ')
        self.active_hr_mean.place(relx = 0.31, rely = 0.28)
        self.active_hr = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.active_hr.place(relx = 0.35, rely = 0.27)

        self.active_hr_std = tk.Label(self.parent, text = 'Stds: ')
        self.active_hr_std.place(relx = 0.31, rely = 0.33)
        self.active_hr1 = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.active_hr1.place(relx = 0.35, rely = 0.32)
        

        self.active_dBP_title = tk.Label(self.parent, text = 'dBP:')
        self.active_dBP_title.place(relx = 0.31, rely = 0.45)

        self.active_dBP_mean = tk.Label(self.parent, text = 'Means: ')
        self.active_dBP_mean.place(relx = 0.31, rely = 0.50)
        self.active_dBP = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.active_dBP.place(relx = 0.35, rely = 0.49)

        self.active_dBP_std = tk.Label(self.parent, text = 'Stds: ')
        self.active_dBP_std.place(relx = 0.31, rely = 0.55)
        self.active_dBP1 = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.active_dBP1.place(relx = 0.35, rely = 0.54)
        
        self.active_sBP_title = tk.Label(self.parent, text = 'sBP:')
        self.active_sBP_title.place(relx = 0.31, rely = 0.65)

        self.active_sBP_mean = tk.Label(self.parent, text = 'Means: ')
        self.active_sBP_mean.place(relx = 0.31, rely = 0.70)
        self.active_sBP = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.active_sBP.place(relx = 0.35, rely = 0.69)

        self.active_sBP_std = tk.Label(self.parent, text = 'Stds: ')
        self.active_sBP_std.place(relx = 0.31, rely = 0.75)
        self.active_sBP1 = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.active_sBP1.place(relx = 0.35, rely = 0.74)
        
# Tilt gui labels
        self.tilt_hr_title = tk.Label(self.parent, text = 'HR')
        self.tilt_hr_title.place(relx = 0.61, rely = 0.25)

        self.tilt_hr_mean = tk.Label(self.parent, text = 'Means: ')
        self.tilt_hr_mean.place(relx = 0.61, rely = 0.28)
        self.tilt_hr = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.tilt_hr.place(relx = 0.65, rely = 0.27)

        self.tilt_hr_std = tk.Label(self.parent, text = 'Stds: ')
        self.tilt_hr_std.place(relx = 0.61, rely = 0.33)
        self.tilt_hr1 = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.tilt_hr1.place(relx = 0.65, rely = 0.32)
        

        self.tilt_dBP_title = tk.Label(self.parent, text = 'dBP:')
        self.tilt_dBP_title.place(relx = 0.61, rely = 0.45)

        self.tilt_dBP_mean = tk.Label(self.parent, text = 'Means: ')
        self.tilt_dBP_mean.place(relx = 0.61, rely = 0.50)
        self.tilt_dBP = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.tilt_dBP.place(relx = 0.65, rely = 0.49)

        self.tilt_dBP_std = tk.Label(self.parent, text = 'Stds: ')
        self.tilt_dBP_std.place(relx = 0.61, rely = 0.55)
        self.tilt_dBP1 = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.tilt_dBP1.place(relx = 0.65, rely = 0.54)
        
        self.tilt_sBP_title = tk.Label(self.parent, text = 'sBP:')
        self.tilt_sBP_title.place(relx = 0.61, rely = 0.65)

        self.tilt_sBP_mean = tk.Label(self.parent, text = 'Means: ')
        self.tilt_sBP_mean.place(relx = 0.61, rely = 0.70)
        self.tilt_sBP = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.tilt_sBP.place(relx = 0.65, rely = 0.69)

        self.tilt_sBP_std = tk.Label(self.parent, text = 'Stds: ')
        self.tilt_sBP_std.place(relx = 0.61, rely = 0.75)
        self.tilt_sBP1 = tk.Label(self.parent, text = '', font = ('Arial', 15))
        self.tilt_sBP1.place(relx = 0.65, rely = 0.74)


    
    def save(self):
        if self.option_value.get() == '':
            print('please select a date + id in the dropdown menu')
            return
        df = self.database.to_dataframe()
        dateid = self.option_value.get().split(',')
        date = dateid[0]
        pid = dateid[1]
                
        
        df_date = df.apply(lambda _ : date, axis = 1)
        df_id = df.apply(lambda _ : pid, axis = 1)
        df_date = df_date.rename('date')
        df_id = df_id.rename('ID')
        df = pd.concat([df, df_date, df_id], axis = 1)

        
        if os.path.isfile(save_path):
            ori_df = pd.read_csv(save_path, dtype = 'O')

            del_df = ori_df[(ori_df['ID'] == pid) & (ori_df['date'] == date)]

            if not del_df.empty:
                ori_df = ori_df.drop(del_df.index.values)
            new_df = pd.concat([ori_df, df])
            
            new_df.to_csv(save_path, mode='w', index=False)
        else:
            df.to_csv(save_path, mode='a', index=False)

    
    def reset(self):
# Sinudal reset
        self.label_hr.config(text = '')
        self.label_hr1.config(text = '')
        self.label_dBP.config(text = '')
        self.label_dBP1.config(text = '')
        self.label_sBP.config(text = '')
        self.label_sBP1.config(text = '')

# Active reset
        self.active_hr.config(text = '')
        self.active_hr1.config(text = '')
        self.active_dBP.config(text = '')
        self.active_dBP1.config(text = '')
        self.active_sBP.config(text = '')
        self.active_sBP1.config(text = '')

# Tilt reset
        self.tilt_hr.config(text = '')
        self.tilt_hr1.config(text = '')
        self.tilt_dBP.config(text = '')
        self.tilt_dBP1.config(text = '')
        self.tilt_sBP.config(text = '')
        self.tilt_sBP1.config(text = '')
        
    def update(self,event):
        """
        Function to read database, check for valid existing data, display and 
        save to values if they exist. This way it allows one to finish some
        and come back later to finish all data, or the need to only correct a 
        little.
        """

        dateid = self.option_value.get().split(',')
        date = dateid[0]
        pid = dateid[1]

        self.database.reset()

        self.reset()

        if os.path.exists(save_path):
            valid = pd.read_csv(save_path)
            
            data = valid.loc[(valid['date'] == date) & (valid['ID'] == pid)]
            data = data.reset_index(drop = True)
            self.test = data
            if data.index.size == 0:
                return

            resp_data = data.iloc[:,:27]
#            self.database.resp_data = resp_data.apply(resp_data_order)
            l = len(resp_data)
            res = {}
            if l == 1:
                res = resp_data.to_dict()
            elif l == 2:
                res = {}
                for k in resp_data.keys():
                    res.update({k : [resp_data[k][0], resp_data[k][1]]})
            elif l == 3:
                res = {}
                for k in resp_data.keys():
                    res.update({k : [resp_data[k][0], resp_data[k][1], resp_data[k][2]]})
            
            self.database.resp_data = res
            self.database.active_data = data.iloc[0,27 : 37].to_dict()
            self.database.tilt_data = data.iloc[0, 37 : 44].to_dict()
            
            
            
# sindual gui update
            self.label_hr.config(text = round_lists(list(data['resperation_hr_mean'])))
            self.label_hr1.config(text = round_lists(list(data['resperation_hr_std'])))
            self.label_dBP.config(text = round_lists(list(data['resperation_dBP_mean'])))
            self.label_dBP1.config(text = round_lists(list(data['resperation_dBP_std'])))

            self.label_sBP.config(text = round_lists(list(data['resperation_sBP_mean'])))
            self.label_sBP1.config(text = round_lists(list(data['resperation_sBP_std'])))

# Active gui update
            if type(data['Active standing_hr_mean'].iloc[0]) == str:
                self.active_hr.config(text = round_lists(convert_string_list_to_list(
                                            data['Active standing_hr_mean'].iloc[0])))
                self.active_hr1.config(text = round_lists(convert_string_list_to_list(
                                            data['Active standing_hr_std'].iloc[0])))
                self.active_dBP.config(text = round_lists(convert_string_list_to_list(
                                            data['Active standing_dBP_mean'].iloc[0])))
                self.active_dBP1.config(text = round_lists(convert_string_list_to_list(
                                            data['Active standing_dBP_std'].iloc[0])))

                self.active_sBP.config(text = round_lists(convert_string_list_to_list(
                                            data['Active standing_sBP_mean'].iloc[0])))
                self.active_sBP1.config(text = round_lists(convert_string_list_to_list(
                                            data['Active standing_sBP_std'].iloc[0])))

# Tilt gui update
            if type(data['Tilt_hr_mean'].iloc[0]) == str:
                self.tilt_hr.config(text = round_lists(convert_string_list_to_list(
                                            data['Tilt_hr_mean'].iloc[0])))
                self.tilt_hr1.config(text = round_lists(convert_string_list_to_list(
                                            data['Tilt_hr_std'].iloc[0])))
                self.tilt_dBP.config(text = round_lists(convert_string_list_to_list(
                                            data['Tilt_dBP_mean'].iloc[0])))
                self.tilt_dBP1.config(text = round_lists(convert_string_list_to_list(
                                            data['Tilt_dBP_std'].iloc[0])))

                self.tilt_sBP.config(text = round_lists(convert_string_list_to_list(
                                            data['Tilt_sBP_mean'].iloc[0])))
                self.tilt_sBP1.config(text = round_lists(convert_string_list_to_list(
                                            data['Tilt_sBP_std'].iloc[0])))            
            
    
    def run_choose(self, choice):
        if self.option_value.get() == '':
            print('please select a date + id in the dropdown menu')
            return
        try:
            dateid = self.option_value.get().split(',')
            date = dateid[0]
            pid = dateid[1]
        except:
            print('information in the dropdown menu is either incorrect or has been corrupted')
            return
        if choice == 'Active standing':
            active_stand_rects = choose_data.create_dicts(
                ['Active standing', 'Active standing'],
                [30, 30], 
                [180, -30],
                ['lime','deepskyblue'],
                [10, 10])
    
            exam_data, mark_data = choose_data.read_examination(date, pid, 'Active standing', read_path)

            (sBP_mean_list, sBP_std_list, dBP_mean_list, dBP_std_list, 
                    hr_mean_list, hr_std_list, areas_list) = choose_data.plot_data(exam_data, mark_data, active_stand_rects)
            mark_start = mark_data.loc[mark_data['Mark'] == 'Active standing']['Time'].values[0]
            res = choose_data.find_min_max(exam_data,mark_start, 15)

            self.active_hr.config(text = round_lists(hr_mean_list))
            self.active_hr1.config(text = round_lists(hr_std_list))
            self.active_dBP.config(text = round_lists(dBP_mean_list))
            self.active_dBP1.config(text = round_lists(dBP_std_list))

            self.active_sBP.config(text = round_lists(sBP_mean_list))
            self.active_sBP1.config(text = round_lists(sBP_std_list))

        elif choice == 'Tilt':
            exam_data, mark_data = choose_data.read_examination(date, pid, 'Tilt', read_path, 240)
        
            one_fourth = (mark_data['Time'][1] - mark_data['Time'][0])/4

    # extra rectangle for after tilt down  x<  
            tilt_rects = choose_data.create_dicts(
                ['Tilt up', 'Tilt up', 'Tilt down', 'Tilt down'],
                [30, 30, 30, 30], 
                [-40, one_fourth, -one_fourth, 40], 
                ['deepskyblue', 'limegreen', 'c', 'darkorchid'], 
                [30, 30, 30, 30])

            (sBP_mean_list, sBP_std_list, dBP_mean_list, dBP_std_list, 
                    hr_mean_list, hr_std_list, areas_list) = choose_data.plot_data(exam_data, mark_data, tilt_rects)
            res = {}

            self.tilt_hr.config(text = round_lists(hr_mean_list))
            self.tilt_hr1.config(text = round_lists(hr_std_list))
            self.tilt_dBP.config(text = round_lists(dBP_mean_list))
            self.tilt_dBP1.config(text = round_lists(dBP_std_list))

            self.tilt_sBP.config(text = round_lists(sBP_mean_list))
            self.tilt_sBP1.config(text = round_lists(sBP_std_list))


        self.database.update_choose_data(
            sBP_mean_list, sBP_std_list,
            dBP_mean_list, dBP_std_list,
            hr_mean_list, hr_std_list,
            areas_list, choice, res)

    def run_fun(self):
        if self.option_value.get() == '':
            print('please select a date + id in the dropdown menu')
            return
        try:
            dateid = self.option_value.get().split(',')
            date = dateid[0]
            pid = dateid[1]
        except:
            print('information in the dropdown menu is either incorrect or has been corrupted')
            return
        
        resp_datas, resp_times = sin_fit.read_examination(date, pid, 0, read_path)

        if resp_datas == []:
            print('There have been found no proper resperation datas')
            return

        hr_list = []
        sBP_list = []
        dBP_list = []
        hr_fit_list = []
        dBP_fit_list = []
        sBP_fit_list = []        
        for i in range(len(resp_datas)):
            hr, sBP, dBP, hr_fit, sBP_fit, dBP_fit = sin_fit.ten_period_resps(resp_datas[i], resp_times[i], i)
            hr_list.append(hr)
            sBP_list.append(sBP)
            dBP_list.append(dBP)
            hr_fit_list.append(hr_fit)
            sBP_fit_list.append(sBP_fit)
            dBP_fit_list.append(dBP_fit)
        
        self.database.update_resp_data(
            hr_list, 
            hr_fit_list, 
            sBP_list, 
            sBP_fit_list, 
            dBP_list, 
            dBP_fit_list)
    
    
# Most likely inefficient update of data to GUI.
        
        t0 = []
        t1 = []
        t2 = []
        t3 = []
        
        d0 = []
        d1 = []
        s0 = []
        s1 = []
        for i in range(len(hr_list)):
            t0.append(round(hr_list[i][0], 2))
            t1.append(round(hr_list[i][1],2))
            t2.append(round_lists(hr_list[i][2]))
            t3.append(round_lists(hr_list[i][3]))
            
            d0.append(round(dBP_list[i][0],2))
            d1.append(round(dBP_list[i][1],2))
            
            s0.append(round(sBP_list[i][0],2))
            s1.append(round(sBP_list[i][1],2))
            
        self.label_hr.config(text = t0)
        self.label_hr1.config(text = t1)
        self.label_dBP.config(text = d0)
        self.label_dBP1.config(text = d1)

        self.label_sBP.config(text = s0)
        self.label_sBP1.config(text = s1)

def round_lists(l, decs = 2):
    return list(map(lambda x : round(x, decs), l))


def t(x):
    if x == ',' or x == ']':
        return True
    elif x == '[':
        return False
    else:
        return None

def convert_string_list_to_list(strs):
    strs = strs.replace(' ', '')
    l = []
    i = ''
    for s in strs:
        b = t(s)
        if b == True:
            l.append(float(i))
            i = ''
        elif b == None:
            i = i + s
    return l

def resp_data_order(x):
    l = len(x)
    if l == 1:
        return x
    elif l == 2:
        return [x[0], x[1]]
    elif l == 3:
        return [x[0], x[1], x[2]]

if __name__ == '__main__':
    vdata = Vip_data()
    window = tk.Tk()

# The database read from.    
    df = pd.read_csv(read_from)
    
    gui = GUIApp(window, vdata, df)

    window.minsize(1200,600)

    window.mainloop()

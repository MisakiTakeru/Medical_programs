#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on 13/1 2025

@author: Joachim Normann Larsen
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from stretch_rect import stretch_rect
from matplotlib.patches import Rectangle

def read_examination(date, pid, exam, ex_time_after = 0):
    """
    
    Parameters
    ----------
    date : string
        date to look up data in the csv files.
    pid : string
        patient id to look up data in the csv files.
    exam : string
        the type of data marking needed for extraction

    Returns
    -------
    resp_datas : list(pd.DataFrame)
        a list including pd.DataFrames of respiration periods.
    resps_times : list(tuples)
        tuple of the time markers for respiration [(resp_start, resp_slut)].

    """

# all of the marker names and their datatypes. Used to reduce the memory
# usage when loading with pandas.
    dtype_dict = {'Time': 'float64', 'Beat': 'float64', 'RRI': 'float64',
     'HR': 'float64', 'sBP': 'float64', 'dBP': 'float64', 'mBP': 'float64',
     'SV': 'float64', 'SI': 'float64', 'CO': 'float64', 'CI': 'float64',
     'TPR': 'float64', 'TPRI': 'float64', 'LFnu-dBP': 'float64',
     'HFnu-dBP': 'float64', 'VLF-dBP': 'float64', 'LF-dBP': 'float64',
     'HF-dBP': 'float64', 'PSD-dBP': 'float64', 'LF/HF-dBP': 'float64',
     'LF/HF': 'float64', 'LFnu-RRI': 'float64', 'HFnu-RRI': 'float64',
     'VLF-RRI': 'float64', 'LF-RRI': 'float64', 'HF-RRI': 'float64',
     'PSD-RRI': 'float64', 'LF/HF-RRI': 'float64', 'Type': 'O',
     'Idx': 'float64', 'Count': 'float64', 'RampUpCount': 'float64',
     'RampDownCount': 'float64', 'Slope': 'float64', 'SlopeMean': 'float64',
     'SlopeMeanUp': 'float64', 'SlopeMeanDown': 'float64', 'Offset': 'float64',
     'SBP_y': 'O', 'DBP_y': 'float64', 'HR_y': 'O', 'Duration': 'O',
     'date': 'O', 'ID': 'O', 'Label': 'O', 'Markers': 'O'}

    data = pd.read_csv('/home/jlar0426/Documents/csv/test.csv', dtype = dtype_dict)
    markers = pd.read_csv('/home/jlar0426/Documents/csv/marks.csv')
    
# Gets data from specific date and id 
    did = data.loc[(data['date'] == date) & (data['ID'] == pid)]
    mid = markers.loc[(markers['date'] == date) & (markers['ID'] == pid)]

    names, counts = np.unique(mid.Mark, return_counts = True)
    
    stop_loc = np.where(names == 'Stop Recording')

# Checks if part is including more than one experiment (only works for 2 as I 
# assume there is a maximum of 2 experiments.), changes the second experiments
# time to ensure that they do not overlap.
    if counts[stop_loc[0][0]] > 1:
        for i in range(len(did)-1):
            c_time = did.iloc[i].Time
            if c_time > did.iloc[i+1].Time:
                split = i+1
                time_add = c_time
                break
        start_add = np.where(mid.Mark == 'Stop Recording')[0][0] + 1
        
        for i in range(start_add, len(mid)):
            mid.at[i,'Time'] += time_add
        for i in range(split, len(did)):
            did.at[i,'Time'] += time_add
            
    try:
        if exam == 'Active standing':
            m1 = mid.loc[mid['Mark'] == 'Active standing']
            m2 = mid.loc[mid['Mark'] == 'Active standing - done']
            
            if m1.empty or m2.empty:
                raise Exception('Missing either Active standing or Active standing - done marker')

        elif exam == 'Tilt':
            m1 = mid.loc[mid['Mark'] == 'Tilt up']
            m2 = mid.loc[mid['Mark'] == 'Tilt down']            

            if m1.empty or m2.empty:
                raise Exception('Missing either Tilt up or Tilt down marker')
    except:
        raise Exception('Missing a crucial mark')

    mark_data = pd.concat([m1, m2])    
    mark_data = mark_data.reset_index(drop = True)
    
    if mark_data.empty:
        return pd.DataFrame([]), mark_data

    exam_data = did.loc[(did['Time'] < (mark_data['Time'][1] + ex_time_after)) & 
                         (did['Time'] > (mark_data['Time'][0] - 45))]

    return exam_data, mark_data


def calculate(data, x0, x1):
    """
    Parameters
    ----------
    data : pd.DataFrame
        all of the data from the given time interval.
    x0 : Float
        start time period to calculate data from.
    x1 : Float
        end time period to calculate data from.

    Returns
    -------
    sBP_mean : Float
        The mean of the systolic blood pressure.
    sBP_std : Float
        The standard deviation of the systolic blood pressure.
    dBP_mean : Float
        The mean of the diastolic blood pressure.
    dBP_std : Float
        The standard deviation of the diastolic blood pressure.
    hr_mean : Float
        The mean of the heart rate.
    hr_std : Float
        The standard deviation of the heart rate.
    """
    
    tdata = data.loc[(data['Time'] >= x0) & (data['Time'] <= x1)] 
    
# There can be nan numbers in our dataset and for mean and std to not return
# nan, we need to drop them.
    tdata = tdata.drop(tdata[np.isnan(tdata.HR)].index)
    tdata = tdata.drop(tdata[np.isnan(tdata.sBP)].index)
    tdata = tdata.drop(tdata[np.isnan(tdata.dBP)].index)
    
    sBP_mean = np.mean(tdata['sBP'].values)
    sBP_std = np.std(tdata['sBP'].values)
    dBP_mean = np.mean(tdata['dBP'].values)
    dBP_std = np.std(tdata['dBP'].values)
    hr_mean = np.mean(tdata['HR'].values)
    hr_std = np.std(tdata['HR'].values)
    
    
    return sBP_mean, sBP_std, dBP_mean, dBP_std, hr_mean, hr_std

def create_draggable_line(ax, time, length, y_pos, offset, 
                          colour = 'deepskyblue', dist_interact = 10):
    
    rect = Rectangle((time + offset, y_pos), length, 
                     0.2, fill = False, color = colour)
    ax.add_patch(rect)
    inter = stretch_rect(rect, dist_interact)
    return inter

def find_min_max(data, mark_start, interval, calc = ['sBP', 'dBP', 'HR'], 
                      minmax = ['min', 'min', 'max']):
    """
    Parameters
    ----------
    data : pd.DataFrame
        The data from which we find the minimum and maximum.
    mark_start : string
        The mark from which we start our interval.
    interval : int
        how long the interval will be
    calc : list, optional
        The data which it will look for. The default is ['sBP', 'dBP', 'HR'].
    minmax : list, optional
        Whether it's going to find the minimum or maximum. 
        The default is ['min', 'min', 'max'].

    Returns
    -------
    calcs : list
        returns the minimum or maximum depended on input of input types.

    """
    if len(calc) != len(minmax):
        print('length of types to use and length of which to gather are different')
        return
    
    # Active stand find within first 15 seconds
    # Carotis find within first 45 seconds.    
    tdata = data.loc[(data['Time'] >= mark_start) & 
                          (data['Time'] <= (mark_start + interval))] 
    
    res = []
    for i in range(len(calc)):
        mm = eval(minmax[i])(tdata[calc[i]].values)
        res.append(mm)
            
#    min_sBP = min(tdata['sBP'].values)
#    min_dBP = min(tdata['dBP'].values)
#    max_HR = max(tdata['HR'].values)
    calc_names = list(map(lambda x, y : x +'_' + y, calc, minmax))
    res = dict(zip(calc_names, res))
    return res

def plot_data(exam_data, mark_data, rectangles, min_max = False, dtypes= [], funcs = []):
    """
    Parameters
    ----------
    exam_data : pd.DataFrame
        All of the data from a given examination.
    mark_data : pd.DataFrame
        The markers made by the examinators during the examination.
    rectangles : list of dict
        List of rectangles needed to be used. One dictionary needs to have a
        mark, length, offset from mark, colour, and distance interaction.
    Returns
    -------
    None.
    """
    
    fig, ax = plt.subplots()
    plt.plot(exam_data['Time'], exam_data['sBP'], label = 'sBP')
    plt.plot(exam_data['Time'], exam_data['HR'], label = 'hr')
    plt.plot(exam_data['Time'], exam_data['dBP'], label = 'dBP')
    
    ax.legend()
    height_start = max(exam_data['sBP'])

    rects = []
    for d in rectangles:
        time = mark_data.loc[mark_data['Mark'] == d['mark']]['Time'].values[0]
        inter = create_draggable_line(ax, time, d['length'], height_start, 
                                      d['offset'], colour = d['colour'], 
                                      dist_interact=d['dist active'])
        rects.append(inter)


    ax.grid(True)
    plt.title('Drag and stretch the line along the x-axis for data to include. Left mouse click to stretch, right mouse click to drag')
    plt.show(block = True)
    

    

    sBP_mean_list = []
    sBP_std_list = []
    dBP_mean_list = []
    dBP_std_list = []
    hr_mean_list = []
    hr_std_list = []
    areas_list = []
# loops through the rectangles and prints out the results
    for i in range(len(rects)):
        x = list(rects[i].get_area())
        x.sort()
        x0 = x[0]
        x1 = x[1]

        sBP_mean, sBP_std, dBP_mean, dBP_std, hr_mean, hr_std = calculate(exam_data, x0, x1)    

        sBP_mean_list.append(sBP_mean)
        sBP_std_list.append(sBP_std)
        dBP_mean_list.append(dBP_mean)
        dBP_std_list.append(dBP_std)
        hr_mean_list.append(hr_mean)
        hr_std_list.append(hr_std)
        areas_list.append(x)

    return (sBP_mean_list, sBP_std_list, dBP_mean_list, dBP_std_list, 
            hr_mean_list, hr_std_list, areas_list)
#        print('Rect ' + str(i) + ' sBP :')
#        print('mean: ', sBP_mean, ', std: ', sBP_std)
#        print('dBP :')
#        print('mean: ', dBP_mean, ', std: ', dBP_std)
#        print('HR :')
#        print('mean: ', hr_mean, ', std: ', hr_std)
#        print('interval: ', x0, x1)
        
    

def create_dicts(marks, lengths, offsets, colours, dists):
    lm = len(marks)
    ll = len(lengths)
    lo = len(offsets)
    lc = len(colours)
    ld = len(dists)
    
    length_test = [ll, lo, lc, ld]
    length_test = [lm == x for x in length_test]
    
    if len(np.unique(length_test)) != 1:
        print('error not all inputs are of equal lengths')
        print('please check that all inputs are of equal lengths')
        return []

    res = []
    for i in range(lm):
        d = {'mark' : marks[i], 'length' : lengths[i], 'offset' : offsets[i],
             'colour' : colours[i], 'dist active' : dists[i]}
        res.append(d)
    
    return res
    
if __name__ == '__main__':

    date = '9/7/2021'
    pid = '021206-5242'

# Active standing test:    
    active_stand_rects = [{'mark' : 'Active standing', 'length' : 30, 
                           'offset' :180, 'colour' : 'lime', 'dist active' : 10},
                          {'mark' : 'Active standing', 'length' : 30, 
                           'offset' :-30, 'colour' : 'deepskyblue', 'dist active' : 10}]
    
    exam_data, mark_data = read_examination(date, pid, 'Active standing')

    plot_data(exam_data, mark_data, active_stand_rects)
    mark_start = mark_data.loc[mark_data['Mark'] == 'Active standing']['Time'].values[0]
    res = find_min_max(exam_data,mark_start, 15)
    print(res)
    
# Tilt test:
    exam_data, mark_data = read_examination(date, pid, 'Tilt', 240)
    
    one_fourth = (mark_data['Time'][1] - mark_data['Time'][0])/4

# extra rectangle for after tilt down  x<  
    tilt_rects = create_dicts(['Tilt up', 'Tilt up', 'Tilt down', 'Tilt down'],
                              [30, 30, 30, 30], 
                              [-40, one_fourth, -one_fourth, 40], 
                              ['deepskyblue', 'limegreen', 'c', 'darkorchid'], 
                              [30, 30, 30, 30])

    res = plot_data(exam_data, mark_data, tilt_rects)
    print(res)
# Carotis test:
#    exam_data, mark_data = read_examination(date, pid, 'Carotis')
    
#    print(exam_data, mark_data)

#    carotis_rects = create_dicts(['Carotis Sin'], [30], [-40], ['deepskyblue'],
#                                 [30])

#    print(carotis_rects)
#    mark_carrots =  mark_data.loc[mark_data['Mark'] == 'Carotis']['Time']
#    for time in mark_carrots:
#    mark_start =  mark_data.loc[mark_data['Mark'] == 'Carotis sin']['Time'].values[0]
#       res = find_min_max(exam_data,time, 45)
#       print(res)
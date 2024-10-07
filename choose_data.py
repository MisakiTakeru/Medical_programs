#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from stretch_rect import stretch_rect
from matplotlib.patches import Rectangle
from sin_fit import interact

def read_examination(date, pid, exam):
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
    data = pd.read_csv('/home/jlar0426/Documents/csv/test.csv')
    markers = pd.read_csv('/home/jlar0426/Documents/csv/marks.csv')

# Gets data from specific date and id 
    did = data.loc[(data['date'] == date) & (data['ID'] == pid)]
    mid = markers.loc[(markers['date'] == date) & (markers['ID'] == pid)]

    mark_data = mid[mid['Mark'].str.contains(exam)]
    mark_data = mark_data.reset_index(drop = True)

    print(mark_data['Time'][0])
    print(type(mark_data['Time'][0]))

    exam_data = did.loc[(did['Time'] < mark_data['Time'][1]) & 
                         (did['Time'] > (mark_data['Time'][0] - 45))]

    return exam_data, mark_data


def calculate(data, x0, x1):
    tdata = data.loc[(data['Time'] >= x0) & (data['Time'] <= x1)] 
    
    sBP_mean = np.mean(tdata['sBP'].values)
    sBP_std = np.std(tdata['sBP'].values)
    dBP_mean = np.mean(tdata['dBP'].values)
    dBP_std = np.std(tdata['dBP'].values)
    hr_mean = np.mean(tdata['HR'].values)
    hr_std = np.std(tdata['HR'].values)
    
    return sBP_mean, sBP_std, dBP_mean, dBP_std, hr_mean, hr_std

def plot_data(exam_data, mark_data):
    fig, ax = plt.subplots()
    plt.plot(exam_data['Time'], exam_data['sBP'], label = 'sBP')
    plt.plot(exam_data['Time'], exam_data['HR'], label = 'hr')
    plt.plot(exam_data['Time'], exam_data['dBP'], label = 'dBP')
    
    ax.legend()
    rect = Rectangle((mark_data['Time'][0] + 180,max(exam_data['sBP'])), 30, 0.2, fill = False, color = 'r')
    ax.add_patch(rect)
    inter = stretch_rect(rect, 10)
    
    
    rect2 = Rectangle((mark_data['Time'][0] - 30, min(exam_data['sBP'])), 30, 0.2, fill = False, color = 'c')
    ax.add_patch(rect2)
    inter2 = stretch_rect(rect2, 10)

    ax.grid(True)
    plt.title('Drag and stretch the line along the x-axis for data to include. Left mouse click to stretch, right mouse click to drag')
    plt.show(block = True)
    
    mark_start = mark_data.loc[mark_data['Mark'] == 'Active standing']['Time'].values[0]
    
    tdata = exam_data.loc[(exam_data['Time'] >= mark_start) & 
                          (exam_data['Time'] <= (mark_start + 20))] 
    
    min_sBP = min(tdata['sBP'].values)
    min_dBP = min(tdata['dBP'].values)
    max_HR = max(tdata['HR'].values)
    
    x = list(inter.get_area())
    
    x.sort()
    x0 = x[0]
    x1 = x[1]

    sBP_mean, sBP_std, dBP_mean, dBP_std, hr_mean, hr_std = calculate(exam_data, x0, x1)    

    print('Rect 1 sBP :')
    print('mean: ', sBP_mean, ', std: ', sBP_std)
    print('dBP :')
    print('mean: ', dBP_mean, ', std: ', dBP_std)
    print('HR :')
    print('mean: ', hr_mean, ', std: ', hr_std)
    print('interval: ', x0, x1)
    
    x_2 = list(inter2.get_area())
    x_2.sort()
    x2 = x_2[0]
    x3 = x_2[1]
    sBP_mean, sBP_std, dBP_mean, dBP_std, hr_mean, hr_std = calculate(exam_data, x2, x3)    

    print('Rect 2 sBP :')
    print('mean: ', sBP_mean, ', std: ', sBP_std)
    print('dBP :')
    print('mean: ', dBP_mean, ', std: ', dBP_std)
    print('HR :')
    print('mean: ', hr_mean, ', std: ', hr_std)
    print('interval: ', x2, x3)
    
    print('minima:')
    print('dBP: ', min_dBP, 'sBP', min_sBP)
    print('maxima')
    print('HR', max_HR)
if __name__ == '__main__':

    date = '2019-06-13'
    pid = '020698-0931'
    
    exam_data, mark_data = read_examination(date, pid, 'Active standing')

    plot_data(exam_data, mark_data)
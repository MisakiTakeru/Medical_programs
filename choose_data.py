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



if __name__ == '__main__':

    date = '2019-06-13'
    pid = '020698-0931'
    
    exam_data, mark_data= read_examination(date, pid, 'Active standing')

    fig, ax = plt.subplots()
    plt.plot(exam_data['Time'], exam_data['sBP'], label = 'sBP')
    plt.plot(exam_data['Time'], exam_data['HR'], label = 'hr')
    print(exam_data)
    
    ax.legend()
    rect = Rectangle((mark_data['Time'][0] + 180,max(exam_data['sBP'])), 30, 0.2, fill = False)
    ax.add_patch(rect)
    inter = stretch_rect(rect, 5)
    inter2 = interact(rect)

    plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on 3/2 2025

@author: Joachim Normann Larsen
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import matplotlib.patches as patches
from matplotlib.widgets import CheckButtons

def read_examination(date, pid, resp_nr, read_path):
    """
    
    Parameters
    ----------
    date : string
        date to look up data in the csv files.
    pid : string
        patient id to look up data in the csv files.

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

    dtype_marker = {'Time' : 'float64', 'Mark' : 'O', 'date' : 'O', 'ID' : 'O'}

    data = pd.read_csv(read_path + '/test.csv', dtype = dtype_dict)
    markers = pd.read_csv(read_path + '/marks.csv', dtype = dtype_marker)


# Gets data from specific date and id 
    did = data.loc[(data['date'] == date) & (data['ID'] == pid)].reset_index(drop = True)
    mid = markers.loc[(markers['date'] == date) & (markers['ID'] == pid)].reset_index(drop = True)
    
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

# gets all markers that have resperation in its name
    resps = mid[mid['Mark'].str.contains('resperation')]
    resps = resps.drop_duplicates()
    resps = resps.reset_index(drop = True)
    
    allowed_marks = ['_start_1', '_slut_1','_start_2', '_slut_2','_start_3', '_slut_3']
    df = pd.DataFrame()
    for am in allowed_marks:
        tmp_df = resps.loc[resps['Mark'] == 'resperation' + am]
    
        df = pd.concat([df,tmp_df], ignore_index = True)
    
    
# gets the data within the first resperation interval.
    resp_datas = []
    resp_times = []
    for r in take_two(df['Time'],2):
        resp_times.append(r)
        resp_data_i = resp_areas(r, did)
        resp_datas.append(resp_data_i)
#    resperation_1 = did.loc[(did['Time'] > resps['Time'][resp_nr]) & (did['Time'] < resps['Time'][resp_nr+1])]
    return resp_datas, resp_times

#    no_nan = resperation_1.drop(resperation_1[np.isnan(resperation_1.HR)].index)

def fit_resp(exam, data_type):
    """

    Parameters
    ----------
    exam : pd.DataFrame
        dataframe consisting of data for a respiration_start to respiration_slut period.

    Returns
    -------
    param which are the found optmized parameters for the sine function as well as,
    it's covariance variables the time data interval and the heart rate for that 
    time interval.

    """

    hr = exam[data_type].values

    time = exam['Time'].values

    nan_indices = np.argwhere(np.isnan(hr))
    
    if np.size(nan_indices) != 0:
        for ind in -np.sort(-nan_indices[:,0]):
            hr = np.delete(hr, ind)
            time = np.delete(time, ind)

    ff = fftfreq(len(time), (time[1] - time[0]))

    Fyy = abs(fft(hr))

# Finds frequency based on the peak frequency of the frequency space (fourier transformation)
#    B_guess_fft = 2*np.pi*abs(ff[np.argmax(Fyy[1:])+1])
# Finds frequency based on knwoledge that there will be 6 minima and maximas over the duration
#    B_guess = (2*np.pi)/((resps['Time'][resp_nr+1]-resps['Time'][resp_nr])/6)
    B_guess_fft = (2*np.pi)/(time[-1] - time[0])

# Finds the positions of the maximas and minimas
    maximas, _ = find_peaks(hr)
    minimas, _ = find_peaks(hr*(-1))
    
# Amplitude guess based on all of the maxima and minima
    A_guess = (np.mean(hr[maximas]) - np.mean(hr[minimas]))/2
# Set in for the 10 second intervals calculations as there is not always enough well defined points to find peaks.
    A_guess = np.std(hr) * 2.**0.5

# horisontal shift (Don't know what to guess, as it will be the direction the curve is starting at)
    h_guess = 0
# vertical shift uses mean of the dataset as it is the indicator of the mid point of the sinus curve when np.sin() == 0
    v_guess = np.mean(hr)

# setting maxfev (max function evaluations) up to 10000 since I have hit situations where it failed since it took more than 1000 steps.

    param, p_cov = curve_fit(sine, time, hr, p0 = [A_guess, B_guess_fft, h_guess, v_guess], maxfev = 10000)
    return param, p_cov, time, hr

# function to loop over to take two elements at a time
def take_two(l, n):
    return zip(*[iter(l)]*n)

# Function to get a resperation number and drops all Nans.
def resp_areas(resps_t, did):
    resp_nr = did.loc[(did['Time'] > resps_t[0]) & (did['Time'] < resps_t[1])]
    return resp_nr.drop(resp_nr[np.isnan(resp_nr.HR)].index)

def sine(t,A, B, H, V):
    return A * np.sin(B*t + H) + V

def ten_period_resps(resp_data, markers, res_num):
    """
    

    Parameters
    ----------
    resp_data : pd.DataFrame
        All data point for a given respiration period.
    markers : tuple
        A Tuple of when the period starts and end.
    res_num : int
        indication of which respiration is currently being worked on.

    Returns
    -------
    Returns the mean, std, maximas, minimas and which 10 sec intervals to be
    included of the HR, sBP, dBP and their fits.

    """
    
    t_range = np.arange(markers[0], markers[1], (markers[1] - markers[0]) / 6)
    t_range = np.append(t_range,markers[1])
    
    lambd_tuplify = lambda i : (t_range[i], t_range[i+1])
    
    time_defs = list(map(lambd_tuplify, range(len(t_range) -1)))

# Splits a repsiration period into 6 closely equal size periods.
    def split_data(interval):
        return resp_data.loc[(resp_data['Time'] > (interval[0] - 0.5)) & 
                             (resp_data['Time'] < (interval[1] + 0.5))]

    periods = list(map(split_data, time_defs))

    hr_data, hr_fit, chosen = image_handling(periods, res_num, 'HR')

    
    mean_hr, std_hr, maxs_hr, mins_hr = data_handling(hr_data, chosen)
    mean_hr_fit, std_hr_fit, maxs_hr_fit, mins_hr_fit = data_handling(hr_fit, chosen)

    bp_data, bp_fit, chosen1, chosen2 = image_handling(periods, res_num, ['sBP', 'dBP'])
    
    sBP_data = [x[0] for x in bp_data]
    dBP_data = [x[1] for x in bp_data]

    mean_sBP, std_sBP, maxs_sBP, mins_sBP = data_handling(sBP_data, chosen1)
    mean_dBP, std_dBP, maxs_dBP, mins_dBP = data_handling(dBP_data, chosen2)

    sBP_fit = [x[0] for x in bp_fit]
    dBP_fit = [x[1] for x in bp_fit]

    mean_sBP_fit, std_sBP_fit, maxs_sBP_fit, mins_sBP_fit = data_handling(sBP_fit, chosen1)
    mean_dBP_fit, std_dBP_fit, maxs_dBP_fit, mins_dBP_fit = data_handling(dBP_fit, chosen2)

    hr = [mean_hr, std_hr, maxs_hr, mins_hr, chosen]
    sBP = [mean_sBP, std_sBP, maxs_sBP, mins_sBP, chosen1]
    dBP = [mean_dBP, std_dBP, maxs_dBP, mins_dBP, chosen2]
    
    hr_stats_fit = [mean_hr_fit, std_hr_fit, maxs_hr_fit, mins_hr_fit]
    sBP_stats_fit = [mean_sBP_fit, std_sBP_fit, maxs_sBP_fit, mins_sBP_fit]
    dBP_stats_fit = [mean_dBP_fit, std_dBP_fit, maxs_dBP_fit, mins_dBP_fit]


    return hr, sBP, dBP, hr_stats_fit, sBP_stats_fit, dBP_stats_fit

def image_handling(periods, res_num, used_data):
    """
    

    Parameters
    ----------
    periods : list
        list of 10 second interval pd.DataFrame points.
    res_num : int
        Which respiration part is currently being worked on.
    used_data : str or list
        Which data is currently being used to plot data.

    Returns
    -------
    data : list
        list of original datapoints from the 10 seconds intervals.
    data_fit : list
        list of the fitted datapoints from the 10 seconds intervals
    chosen : list
        list of which 10 seconds intervals to be included for calculations.

    """
    fig, ax = plt.subplots(height_ratios = [0.1])
    ax.grid(True)

    tit = 'respiration %d' % (res_num+1)
    plt.title(tit)

    data = []
    data_fit = []
    for i in range(len(periods)):
        if type(used_data) == str:
            param, p_cov, time, d_type = fit_resp(periods[i], used_data)

            data.append(d_type)
            t2 = np.arange(time[0],time[-1], 0.1)
            t = sine(t2, param[0], param[1], param[2], param[3])
            data_fit.append(t)

            ax.plot(t2, t, label = 'line' + str(i))
            ax.scatter(time, d_type)
            fig.text(0.30, 0.22, 'hr')

        elif type(used_data) == list:
            param1, p_cov1, time1, d_type1 = fit_resp(periods[i], used_data[0])
            param2, p_cov2, time2, d_type2 = fit_resp(periods[i], used_data[1])

            data.append([d_type1, d_type2])
            time2 = np.arange(time1[0],time1[-1], 0.1)
            
            t1 = sine(time2, param1[0], param1[1], param1[2], param1[3])            
            t2 = sine(time2, param2[0], param2[1], param2[2], param2[3])
            data_fit.append([t1, t2])
            
            ax.plot(time2, t1, label = 'line' + str(i) + ', ' + used_data[0])
            ax.scatter(time1, d_type1)
            ax.plot(time2, t2, label = 'line' + str(i) + ', ' + used_data[1])
            ax.scatter(time1, d_type2)

            cb_pos2 = plt.axes([0.50, 0.05, 0.15, 0.15])
            check_but2 = CheckButtons(cb_pos2, list(range(6)), actives = [True, True, True, True, True, True])
            fig.text(0.50, 0.22, 'dBP')
            fig.text(0.30, 0.22, 'sBP')


    ax.legend()
    ax.set_position([0.1,0.30, 0.85, 0.65])

    cb_pos = plt.axes([0.30, 0.05, 0.15, 0.15])
    fig.text(0.15, 0.15, 'Which needs to be used?')
    check_but = CheckButtons(cb_pos, list(range(6)), actives = [True, True, True, True, True, True])

    plt.show(block = True)
    
    chosen = check_but.get_status()
    if type(used_data) == list:
        chosen2 = check_but2.get_status()
        return data, data_fit, chosen, chosen2
    return data, data_fit, chosen


def data_handling(data, chosen):
    """
    

    Parameters
    ----------
    data : list
        list of the six 10 seconds intervals datapoints.
    chosen : list
        list of which data intervals to be included in calculations.

    Returns
    -------
    mean : float.
        mean of the data maximas - minimas included.
    std : float
        std of the data maximas - minimas included.
    maxs : list
        list of the maxima points found for each 10 second interval.
    mins : list
        list of the minima points found for each 10 second interval.

    """
    used_data = np.array([])
    maxs = []
    mins = []
    for i in range(len(chosen)):
        d = data[i]
        maxima = max(d)
        minima = min(d)
        used_data = np.append(used_data, [maxima - minima])
        maxs.append(maxima)
        mins.append(minima)
    used_data = used_data * chosen
    used_data = [d for d in used_data if d != 0]
    mean = np.mean(used_data)
    std = np.std(used_data)
    return mean, std, maxs, mins
    
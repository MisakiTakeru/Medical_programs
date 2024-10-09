#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import matplotlib.patches as patches
import functools
from matplotlib.widgets import CheckButtons, TextBox

def read_examination(date, pid, resp_nr):
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
    data = pd.read_csv('/home/jlar0426/Documents/csv/test.csv')
    markers = pd.read_csv('/home/jlar0426/Documents/csv/marks.csv')

# Gets data from specific date and id 
    did = data.loc[(data['date'] == date) & (data['ID'] == pid)]
    mid = markers.loc[(markers['date'] == date) & (markers['ID'] == pid)]

# gets all markers that have resperation in its name
    resps = mid[mid['Mark'].str.contains('resperation')]
    resps = resps.reset_index(drop = True)
# gets the data within the first resperation interval.
    resp_datas = []
    resp_times = []
    for r in take_two(resps['Time'],2):
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


    ff = fftfreq(len(time), (time[1] - time[0]))

    Fyy = abs(fft(hr))

# Finds frequency based on the peak frequency of the frequency space (fourier transformation)
    B_guess_fft = 2*np.pi*abs(ff[np.argmax(Fyy[1:])+1])
# Finds frequency based on knwoledge that there will be 6 minima and maximas over the duration
#    B_guess = (2*np.pi)/((resps['Time'][resp_nr+1]-resps['Time'][resp_nr])/6)

 #   print(B_guess_fft, B_guess)

# Finds the positions of the maximas and minimas
    maximas, _ = find_peaks(hr)
    minimas, _ = find_peaks(hr*(-1))
    
#    print('extremas')
#    print(minimas, maximas)

# Amplitude guess based on all of the maxima and minima
    A_guess = (np.mean(hr[maximas]) - np.mean(hr[minimas]))/2
# Set in for the 10 second intervals calculations as there is not always enough well defined points to find peaks.
    A_guess = np.std(hr) * 2.**0.5

# horisontal shift (Don't know what to guess, as it will be the direction the curve is starting at)
    h_guess = 0
# vertical shift uses mean of the dataset as it is the indicator of the mid point of the sinus curve when np.sin() == 0
    v_guess = np.mean(hr)

#    print([A_guess, B_guess_fft, h_guess, v_guess])
#    try:
# setting maxfev (max function evaluations) up to 10000 since I have hit situations where it failed since it took more than 1000 steps.
    param, p_cov = curve_fit(sine, time, hr, p0 = [A_guess, B_guess_fft, h_guess, v_guess], maxfev = 10000)
#    except:
#        plt.scatter(time,hr)
#        t_new = np.arange(time[0],time[-1], 0.0001)
#        B_guess = (2*np.pi)/((time[-1]-time[0]))
#        print(B_guess, B_guess_fft)
#        param, p_cov = curve_fit(sine, time, hr, p0 = [A_guess, B_guess_fft, h_guess, v_guess])
#        plt.plot(t_new, t2)
#        plt.show()

#    param2, p_cov2 = curve_fit(sine, time, hr, p0 = [A_guess, B_guess, h_guess, v_guess])
#    print(param, param2)
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

def sine_fit(t, y, m_start, m_end):
    """
    Parameters
    ----------
    t : ndarray
        input of the x-axis value.
    y : ndarray
        input of the y-axis values.
    m_start : float
        time for marker start.
    m_end : float
        time for marker end.

    Returns
    -------
    p_opt : ndarray
        optimized parameters of the sine function.
    p_cov : ndarray
        covariance of the optimized parameters.

    """
    maximas, _ = find_peaks(y, distance = 7)
    minimas, _ = find_peaks(y*(-1), distance = 7)
    A_guess = (np.mean(y[maximas]) - np.mean(y[minimas]))/2
    B_guess = (2*np.pi)/((m_end - m_start)/6)
    H_guess = 0
    V_guess = np.mean(y)
    guess = [A_guess, B_guess, H_guess, V_guess]
    p_opt, p_cov = curve_fit(sine, t, y, p0 = guess)
    return p_opt, p_cov

# Function to calculate the 12 areas for one respiration to look for data for calculating mean and std
def set_rois(data, time):
    maximas, _ = find_peaks(data, distance = 7)
    minimas, _ = find_peaks(data*(-1), distance = 7)
#    print(maximas)
#    print(time[maximas])
#    print(minimas)
#    print(time[minimas])
    return
    
# testing plot interacbility
# A class for dragging patches in an axis. Input is a list of draggable patches.
class interact(object):
    
    def __init__(self, artists):
        self.artists = artists
        self.fig = self.artists[0].figure
        self.ax = self.artists[0].axes

        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

        self.dragging = False
        self.selecting = False
        self.selected_artists = []
        self.rect = patches.Rectangle((700,55),20,20, linestyle = '--', fill = False)
        self.ax.add_patch(self.rect)
        self.rect.set_visible(False)
        self.x0 = 0
        self.y0 = 0
        
    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        elif event.button == 1:
            is_on_artist = False
            self.x0, self.y0 = (event.xdata, event.ydata)
            for art in self.artists:
                if art.contains(event)[0]:
                    is_on_artist = art
            if is_on_artist:
# If there has been pressed on an patch outside the box we only want the new one.
                if is_on_artist not in self.selected_artists:
                    self.selected_artists = [is_on_artist]
# self.offset need to be calculated as an array since it is possible there are multiple
                self.offset = ([p.get_xy() for p in self.selected_artists] - 
                               np.array((event.xdata, event.ydata)))
                self.dragging = True
                return
            else:
                self.selecting = True
                self.selected_artists = []

    def on_motion(self, event):
        if event.inaxes != self.ax:
            return
        if self.selecting:
            self.rect.set_visible(True)
            xlim, ylim = self.sort_pos(event.xdata, event.ydata)
            self.rect.set_xy((xlim[0], ylim[0]))
            self.rect.set_width(np.diff(xlim)[0])
            self.rect.set_height(np.diff(ylim)[0])

        elif self.dragging:
            new_loc = np.array((event.xdata, event.ydata)) + self.offset
            for i,art in enumerate(self.selected_artists):
                art.set_xy(new_loc[i])
        self.fig.canvas.draw_idle()

    def on_release(self, event):
        if event.inaxes != self.ax:
            return
        if self.selecting:
            xlim, ylim = self.sort_pos(event.xdata, event.ydata)
            for art in self.artists:
                cx, cy = art.get_center()
                if (cx >= xlim[0] and cx < xlim[1] and cy >= ylim[0] 
                    and cy < ylim[1]):
                        self.selected_artists.append(art)
            self.selecting = False
            self.rect.set_visible(False)
            self.fig.canvas.draw_idle()
        elif self.dragging:
            self.dragging = False

        
    def sort_pos(self, x, y):
        self.x1, self.y1 = (x, y)
        xlim = np.sort([self.x0, self.x1])
        ylim = np.sort([self.y0, self.y1])
        return xlim, ylim        


def ten_period_resps(resp_data, markers, res_num):
#    print(resp_data)
#    print(markers)
    
    t_range = np.arange(markers[0], markers[1], (markers[1] - markers[0]) / 6)
    t_range = np.append(t_range,markers[1])
#    print(markers)
#    print(t_range)
    
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
    used_data = np.array([])
    maxs = []
    mins = []
    for i in range(len(chosen)):
#        d = data[i] * chosen[i]
        d = data[i]
        maxima = max(d)
        minima = min(d)
        used_data = np.append(used_data, [maxima - minima])
        maxs.append(maxima)
        mins.append(minima)
    used_data = used_data * chosen
    used_data = [d for d in used_data if d != 0]
    print(used_data)    
    mean = np.mean(used_data)
    std = np.std(used_data)
    return mean, std, maxs, mins
    
if __name__ == '__main__':

    date = '2019-06-13'
    pid = '020698-0931'
    
    resp_datas, resp_times = read_examination(date, pid, 0)
    for i in range(3):
#        param, p_cov, time, hr = fit_resp(resp_datas[i])
        hr, sBP, dBP, hr_fit, sBP_fit, dBP_fit = ten_period_resps(resp_datas[i], resp_times[i], i)
        print('respiration ' + str(i+1))
        print('hr')
        print(hr)
        print('hr fit')
        print(hr_fit)
        print('sBP')
        print(sBP)
        print('sBP fit')
        print(sBP_fit)
        print('dBP')
        print(dBP)
        print('dBP fit')
        print(dBP_fit)


    t = patches.Rectangle((700,60), 20, 20, fill = False)
    rect = patches.Rectangle((650,60), 20, 20, fill = False)

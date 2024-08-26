#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import matplotlib.patches as patches


def read_examination(date, pid):
    data = pd.read_csv('/home/jlar0426/Documents/csv/test.csv')
    markers = pd.read_csv('/home/jlar0426/Documents/csv/marks.csv')

# Gets data from specific date and id 
    did = data.loc[(data['date'] == date) & (data['ID'] == pid)]
    mid = markers.loc[(markers['date'] == date) & (markers['ID'] == pid)]

# gets all markers that have resperation in its name
    resps = mid[mid['Mark'].str.contains('resperation')]
    resps = resps.reset_index(drop = True)
# gets the data within the first resperation interval.
    resperation_1 = did.loc[(did['Time'] > resps['Time'][0]) & (did['Time'] < resps['Time'][1])]

    no_nan = resperation_1.drop(resperation_1[np.isnan(resperation_1.HR)].index)

    hr = no_nan['HR'].values

    time = no_nan['Time'].values

    ff = fftfreq(len(time), (time[1] - time[0]))

    Fyy = abs(fft(hr))

# Finds frequency based on the peak frequency of the frequency space (fourier transformation)
    B_guess_fft = 2*np.pi*abs(ff[np.argmax(Fyy[1:])+1])
# Finds frequency based on knwoledge that there will be 6 minima and maximas over the duration
    B_guess = (2*np.pi)/((resps['Time'][1]-resps['Time'][0])/6)

# Finds the positions of the maximas and minimas
    maximas, _ = find_peaks(hr)
    minimas, _ = find_peaks(hr*(-1))

# Amplitude guess based on all of the maxima and minima
    A_guess = (np.mean(hr[maximas]) - np.mean(hr[minimas]))/2

# horisontal shift (Don't know what to guess, as it will be the direction the curve is starting at)
    h_guess = 0
# vertical shift uses mean of the dataset as it is the indicator of the mid point of the sinus curve when np.sin() == 0
    v_guess = np.mean(hr)

    param, p_cov = curve_fit(sine, time, hr, p0 = [A_guess, B_guess_fft, h_guess, v_guess])
    param2, p_cov2 = curve_fit(sine, time, hr, p0 = [A_guess, B_guess, h_guess, v_guess])
    return param, p_cov, param2, p_cov2, time, hr

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
    maximas, _ = find_peaks(y)
    minimas, _ = find_peaks(y*(-1))
    A_guess = (np.mean(y[maximas]) - np.mean(y[minimas]))/2
    B_guess = (2*np.pi)/((m_end - m_start)/6)
    H_guess = 0
    V_guess = np.mean(y)
    guess = [A_guess, B_guess, H_guess, V_guess]
    p_opt, p_cov = curve_fit(sine, t, y, p0 = guess)
    return p_opt, p_cov

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
        self.rect = patches.Rectangle((0,0),20,20, linestyle = '--', fill = False)
        self.ax.add_patch(self.rect)
        self.rect.set_visible(False)
        self.x0 = 0
        self.y0 = 0
        
    def on_press(self, event):
        if event.inaxes == self.ax:
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


if __name__ == '__main__':

    date = '2019-06-13'
    pid = '020698-0931'
    
    param, p_cov, param2, p_cov2, time, hr = read_examination(date, pid)
    
    t = sine(time, param[0], param[1], param[2], param[3])
    #t = param[0] * np.sin(param[1]*time+param[2])+param[3]
    t2 = param2[0] * np.sin(param2[1]*time+param2[2])+param2[3]

    fig, ax = plt.subplots()
    ax.plot(time, t, label = 'sine fit')
    ax.scatter(time, hr, label = 'heart rate')
    ax.plot(time, t2, label = 'sine fit 2')
    ax.legend()
    ax.grid(True)


    t = patches.Rectangle((50,100), 20, 20, fill = True)
    rect = patches.Rectangle((100,100), 20, 20, fill = False)
    ax.add_patch(t)
    ax.add_patch(rect)
    inter = interact([t,rect])
# I will be able to create squares for indicating areas where my program will look for data using patches.Rectangle
    
    plt.ion()
    plt.show()
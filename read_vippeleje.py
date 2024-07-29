# -*- coding: utf-8 -*-
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, CheckButtons, TextBox, Button
from enum import Enum
from tkinter.filedialog import askopenfilename
import os
from pathlib import Path
#matplotlib.use('TkAgg')

# In case of this program not opening a new window, but instead tries to 
# show it in IDE run this line in console: %matplotlib qt


# Creates a generator to read a file line by line.
def reader(n):
    for row in open(n, 'r', encoding = 'cp1252'):
        yield row


# Function that splits the Unit lines such that it becomes easier to handle the data.
def read_units(line):
    return re.split(r'\s+', line)[1:-1]

 

"""
---------
read_data
---------
Inputs:
    line -> String
    units -> List

Function for ; separator read line, and removing leading and trailing spaces.
If line is a marker, then it returns the time and the marker as a list.
"""
def read_data(line, units):
    data = re.split(r';', line)
    if len(data) == 1:
        data = re.split(r'   ', line)
        res = [d.strip() for d in data]
    else:
        data = data[:-1]
# Everything is float except for Beat which is a int, and Type which is String.
        converter = lambda key, data : data if key == 'Type' else (
            int(data) if key == 'Beat' else float(data))
        data = [converter(k,d.strip()) for k,d in list(zip(units,data))]
        res =  dict(zip(units,data))
    return res



def assemble_data(name):
    """
    -------------
    assemble_data
    -------------
    Inputs:
        name -> String

    Function that reads the 5 csv files from a single Vippeleje measurement and
    creates a single dataframe containing all of the information.
    """
# a single vippeleje turns into 5 different csv files.
    types = ['.BeatToBeat.csv', '.BPV.csv', '.BRS_BRS0.csv', '.HRV.csv','.OscBP.csv']
    
    allfiles = []
    
    for t in types:
        file = reader(name + t)
    
# skips the first 12 lines, but gets date and patient id from it.
        date = next(file)
        date = date.split(':')[-1].strip()
        for i in range(8):
            next(file)
        
        pid = next(file)
        pid= pid.split(':')[-1].split('_')[0].strip()
        for i in range(2):
            next(file)
    
        units = read_units(next(file))

        l = []
        markings = []
        while True:
# in case of dictionary means we have data of a measurement. In case of a list
# of len 2 means we have a marking, and otherwise they will be empty lines.
            try:
                data = read_data(next(file), units)
                if type(data) == dict:
                    l.append(data)
                elif len(data) == 2:
                    markings.append(list(data))
            except StopIteration:
                print('done')
                break
        allfiles.append(l)

    B2B = pd.DataFrame(allfiles[0])
    BPV = pd.DataFrame(allfiles[1])
    BRS = pd.DataFrame(allfiles[2])
    HRV = pd.DataFrame(allfiles[3])
    OSC = pd.DataFrame(allfiles[4])
        
    B2B_BPV = pd.merge(B2B, BPV, on = ['Time', 'Beat'])
    B2B_BPV_HRV = pd.merge(B2B_BPV, HRV, on= ['Time', 'Beat', 'LF/HF'])
    B2B_BPV_HRV_BRS = pd.merge(B2B_BPV_HRV, BRS, on = 'Time', how = 'left')

# Not including the HR from OSC as it seems suspiciously off compared to actual
# data from the other files. This will result in a HR_x and HR_y 
# (OSC is a different measurement and shall be kept different as well as DBP and SBP)
    B2B_BPV_HRV_BRS_OSC = pd.merge(B2B_BPV_HRV_BRS, OSC, on = 'Time', how = 'outer').sort_values('Time')
    B2B_BPV_HRV_BRS_OSC = B2B_BPV_HRV_BRS_OSC.rename(columns = {'HR_x' : 'HR',
        'SBP' : 'SBP_y', 'DBP' : 'DBP_y'})
    B2B_BPV_HRV_BRS_OSC['date'] = date
    B2B_BPV_HRV_BRS_OSC['ID'] = pid
    return B2B_BPV_HRV_BRS_OSC, markings

# BeatToBeat, BPV and HRV have exact same Time and beat, so those 3 can be merged
# together.
# BeatToBeat and BPV can be merged with 'Time' and 'Beat'. Including
# HRV we will also be merging with 'LF/HF'.
# BRS_BRS0 can be merged with those 3 through a left merge on 'Time'



# An enumerator class to get the expected counts of a Vippeleje experiment.
class expected_counts(Enum):
    StartRecording = 1
    resperation_start = 3
    respitationslut = 3
    Activestanding = 1
    Activestandingdone = 1
    Tiltup = 1
    Tiltdown = 1
    StopRecording = 1
    Nitroglycerine = 1
# I believe this is the marker used, when the patient needed a break so i set it to 1000.
    ContinueDAQafterPause = 1000
    Carotisdxt = 2
    Carotissin = 2
    Carotismassagesin = 2


# There is a single file group that will not be found, to find it without getting
# false positives, I must check the date of the examination.
def check_2_part_file(filepath):
    """

    Inputs:
        filepath -> String
    
    check_2_part_file checks if given filepath is a part 1 of 2 parts examination.
    If true, it will run the assemble data function and return it as a tuple
    # of (True, list of assembled data), otherwise it returns False and a empty list.
    """
    possible = ['del1', '_1', '.1']
    replaces = ['del2', '_2', '.2']
    part_maybe = list(map(lambda x : x in filepath, possible))
    try:
        ind = part_maybe.index(True)
    except Exception:
        return (False, [])
    filepath2 = filepath.replace(possible[ind], replaces[ind])
    file2 = Path(filepath2).is_file()
        
    if file2:
        res2, marks2 = assemble_data(filepath2.replace('.BeatToBeat.csv',''))
        return (True, [res2, marks2])
    else:
        return (False, [])

def get_mark_index(name):
    """
    Parameters
    ----------
    name : string
        mark name to look up.

    Returns
    -------
    i : int or None
        the index for that mark if it exists otherwise None.
    """
    indexes = [i for i,s in enumerate(mark_names) if name in s]
    l = len(indexes)
    if l >= 1:
        if name == 'resperation_slut':
            i = indexes[-1]
        elif 'Carotis' in name:
            return indexes, l
        else:
            i = indexes[0]
    else:
        i = None
    if 'Carotis' in name:
        return [], 0
    return i
        

def make_labels(row):
    """
    Parameters
    ----------
    row : Pandas.Series
        A pandas series gotten through using df.apply(make_labels, axis = 1)

    Returns
    -------
    String
        returns the label given to the dataframe for that row.

    """
    t = row['Time']
    lambd = lambda car: (marks[carotis_i[car]][1], t < (float(marks[carotis_i[car]][0]) + 40) and (
        t > (float(marks[carotis_i[car]][0]) - 10)))
    if t < rest_slut and t > rest_start:
        return 'Rest'
    elif t < active_stand_slut and t > active_stand_start:
        return 'Active standing'
    car_res = list(map(lambd, range(l)))
    car_true = [n for n,b in car_res if b]
    if len(car_true) != 0:
        return car_true[0]
    else:
        return None        

def make_markers(row):
    time = row['Time']
    for i,(t, _) in enumerate(marks):
        if time < float(t):
            return marks[i-1][1]

    
if __name__ == '__main__':
#    name = '/run/user/1000/gvfs/smb-share:server=rghdfsp02,share=dfs/Logget/LovbeskyttetMapper01/Fabry - Collab w KFNM/Data/Vippeleje/011297-0266_schmidt'

#    os.chdir('/run/user/1000/gvfs/smb-share:server=rghdfsp02,share=dfs/Logget/LovbeskyttetMapper01/Fabry - Collab w KFNM/Data/Vippeleje')
    os.chdir('/home/jlar0426/Documents/Temp')
    nameb2b = askopenfilename(filetypes=[('BeattoBeat CSV', '.BeatToBeat.csv')], 
                       initialdir = '/run/user/1000/gvfs/smb-share:server=rghdfsp02,share=dfs/Logget/LovbeskyttetMapper01/Fabry - Collab w KFNM/Data')

    nameb2b = name = os.path.basename(nameb2b)

    name = nameb2b.replace('.BeatToBeat.csv','')
    
    

    res, marks = assemble_data(name)
    
    (boolean, dat) = check_2_part_file(nameb2b)
    if boolean:
        res2 = dat[0]
        marks2 = dat[1]
        res2.Time += res['Time'].iloc[-1]
        marks2 = [[float(t) + res['Time'].iloc[-1], m] for t,m in marks2]
        
        marks = marks + marks2
        res = pd.concat([res, res2], ignore_index=True)

    times, markeds = zip(*marks)
    markeds = list(markeds)

    names, counts = np.unique(markeds, return_counts = True)

    name_add = np.ones(len(counts), dtype = int)

# Adds numbers to those marks that exist duplicates of, and if above expected
# amount gives it a too_many addiction instead.
    for i,n in enumerate(markeds):
        if markeds[i] == 'respitation slut':
            markeds[i] = 'resperation_slut'
        loc = np.where(names == n)[0][0]

        if counts[loc] > 1:
            if expected_counts[n.replace(' ','').replace('-','')].value < name_add[loc]:
                new_name = markeds[i] + '_too_many'
            else:
                new_name = markeds[i] + '_' + str(name_add[loc])
            markeds[i] = new_name
            name_add[loc] += 1


    marks = list(map(list,zip(times, markeds)))

# Drops the nan lines made from the OscBP file, since all OscBP time values
# are different.
    no_nan = res.drop(res[np.isnan(res.HR)].index)
        
# Limiter for maximum y value in graph
    max_y_lim =  max(no_nan[['sBP', 'dBP','HR']].max()) + 5

    fig, ax = plt.subplots(height_ratios = [0.1])

# Plots the 3 different lines: heart rate, systolic- and diatolic blood pressure
    ax.plot(no_nan['Time'], no_nan['HR'], label = 'Heart rate')
    ax.plot(no_nan['Time'], no_nan['sBP'], label = 'systolic')
    ax.plot(no_nan['Time'], no_nan['dBP'], label = 'diatolic')


# Creates the red transparent marker lines, and sets all marker text on the axis.
# Made as a dictionary for easy access when they need to be changed.
    vlines = {}
    for t,n in marks:
        time = float(t)
        vline = ax.vlines(time, 0, max_y_lim, colors = 'r', linestyles = '--')
        
        vtext = ax.text(time, max_y_lim, n, rotation = 20)
        
        vlines[n] = [vline,vtext]


    """
    Creating all of the buttons and it's functions'
    """

# Slider colour, position and size.
    axcolor = 'lightgoldenrodyellow'
    axpos = plt.axes([0.2, 0.001, 0.65, 0.03], facecolor=axcolor)


# Slider that have min val 1 second before first data recording until 500 seconds
# before last recording.
    spos = Slider(axpos, 'Pos', res['Time'][0]-1, res['Time'][res.index[-1]])

# Slider update function, it shows a range of 500 seconds, and just changes
# the area show from the x range.
    def update(val):
        pos = spos.val
        show_ahead = pos + 500
        ax.set_xlim(pos, show_ahead)
        fig.canvas.draw_idle()

    ax.grid(True)
    spos.on_changed(update)
    ax.set_ylim(0, max_y_lim)
    ax.set_xlim(res['Time'][0]-1, res['Time'][0]+500)
    ax.legend(loc = 'upper right')
    ax.set_xlabel('Time / s')
    ax.set_position([0.1,0.45, 0.85, 0.44])
 

# Inserting checkboxes in figure test
    mark_names = [n for _,n in marks]

    cb_pos = plt.axes([0.10, 0.10, 0.15, 0.20])
    fig.text(0.112, 0.32, 'Which marking needs to be changed?')
    check_but = CheckButtons(cb_pos, mark_names)



# When a checkbox is clicked, it will check if there are multiple checked on
# If yes, it will remove everyone except the newest one allowing only 1 checked box.
    def click(label):
        ind = [x for x,y in enumerate(marks) if y[1] == label][0]
        if len(check_but.get_checked_labels()) > 1:
            remove = [i for i, x in enumerate(check_but.get_status()) if x and (ind != i)]
            check_but.set_active(remove[0])
        textbox.set_val(marks[ind][0])

    check_but.on_clicked(click)

    choice_button_pos = plt.axes([0.75, 0.10, 0.05, 0.05])
    fig.text(0.68, 0.125, 'Add number?: ')
    choice_button = CheckButtons(choice_button_pos, ['1','2','3'])

    def click_num(key):
        ind = int(key)-1
        if len(choice_button.get_checked_labels()) > 1:
            remove = [i for i, x in enumerate(choice_button.get_status()) if x and (ind != i)]
            choice_button.set_active(remove[0])

    choice_button.on_clicked(click_num)

# All found marks are:
#array(['', 'Active standing', 'Active standing - done', 'Carotis dxt.',
#       'Carotis sin.', 'Carotismassage sin.', 'Continue DAQ after Pause',
#       'Nitroglycerine', 'Start Recording', 'Stop Recording', 'Tilt down',
#       'Tilt up', 'resperation_start', 'respitation slut'], dtype='<U24')

# seems that dxt means right and sin means left

    markers = ['Start Recording', 'resperation_start', 'resperation_slut', 
              'Active standing', 'Active standing - done', 'Tilt up', 
              'Tilt down', 'Stop Recording', 'Nitroglycerine', 'Carotis sin',
              'Carotis dxt']

    def click_name(name):
        ind = markers.index(name)
        if len(names_button.get_checked_labels()) > 1:
            remove = [i for i, x in enumerate(names_button.get_status()) if x and (ind != i)]
            names_button.set_active(remove[0])

    names_button_pos = plt.axes([0.50, 0.10, 0.15, 0.20])
    fig.text(0.52, 0.32, 'Is it in need of a name change?')
    names_button = CheckButtons(names_button_pos, markers)
    names_button.on_clicked(click_name)


    def add_and_finish(val):
        labels = res.apply(make_labels, axis=1)
        labels = labels.rename('Label')
        res_with_labels = pd.concat([res, labels], axis = 1)
        
        time_marks = res_with_labels.apply(make_markers, axis = 1)
        time_marks = time_marks.rename('Markers')
        res_with_marks = pd.concat([res_with_labels, time_marks], axis = 1)
        print(res_with_marks)
        
        if os.path.isfile('/home/jlar0426/Documents/csv/test.csv'):
            res_with_marks.to_csv('/home/jlar0426/Documents/csv/test.csv', mode='a', index=False, header=False)
        else:
            res_with_marks.to_csv('/home/jlar0426/Documents/csv/test.csv', mode='a', index=False)
        plt.close()
        return 

# textbox position and creation.
    textbox_pos = plt.axes([0.70, 0.20, 0.10, 0.05])
    textbox = TextBox(textbox_pos, 'Time: ')

# apply button position and creation.
    button_pos = plt.axes([0.85, 0.10, 0.10, 0.05])
    finish_button = Button(button_pos, 'Finish and Exit')
    finish_button.on_clicked(add_and_finish)

# Function for a button, which gets the checked name from CheckButtons and
# the text from textbox. If there is a checked name and a possible float in
# the textbox, it will then convert the position of that marker to the new time.
    def change_name(val):

        try:
            time = float(textbox.text)
        except:
            print('What did ya diddily doo?')
#            return

        label = check_but.get_checked_labels()
        if label == []:
            return
        else:
            label = label[0]

# To change text find the dictionary text and use set_x(val)
# To change vline position find the dictionary, remove the vline and insert the new one.
        label_time = [float(t) for t,m in marks if m == label][0]
        if time != label_time:
            vline = vlines[label]
            vline[0].remove()
            vline[0] = ax.vlines(time, 0, max_y_lim, colors = 'r', linestyles = '--')
            vline[1].set_x(time)
            mark_index = [x for x,y in enumerate(marks) if y[1] == label][0]
            marks[mark_index][0] = time
            fig.canvas.draw_idle()
                
        new_name = names_button.get_checked_labels()
        if new_name == []:
            return
        else:
            number = choice_button.get_checked_labels()
            if number == []:
                new_name = new_name[0]
            else:
                new_name = new_name[0] + '_' + number[0]
            
        mark_index = [x for x,y in enumerate(marks) if y[1] == label][0]
        marks[mark_index][1] = new_name
        vlines[label][1].set_text(new_name)
        vlines[new_name] = vlines.pop(label)
        
        ind = [i for i, l in enumerate(check_but.labels) if l.get_text() == label][0]
        check_but.labels[ind].set_text(new_name)

        fig.canvas.draw_idle()        
            

# apply name change.
    button_namechange_pos = plt.axes([0.85, 0.20, 0.10, 0.05])
    name_change_button = Button(button_namechange_pos, 'Apply changes')
    name_change_button.on_clicked(change_name)

    """
    Finished all button creations.
    """
    

    plt.show()

# Starting code for labelling data.

    resp_slut_i = get_mark_index('resperation_slut')
    rest_start = float(marks[resp_slut_i][0]) + 120
    active_stand_i = get_mark_index('Active standing')
    rest_slut = float(marks[active_stand_i][0]) - 120
    
    active_stand_start = float(marks[active_stand_i][0]) - 30
    active_stand_slut = float(marks[active_stand_i][0]) + 240
        
    carotis_i, l = get_mark_index('Carotis')
    



# resperation and respitation keep it as is. I need to create
# a new Label 'Rest' from 60 sec after last resp to 120 sec before Active Standing.
# Active standing needs to be 30 seconds before until 240 seconds after start.
#slet dette: Active standing done, cut the first 120 sec and last 120 sec.
# If Carotis is done, then first 10 seconds, and the first 40 seconds after marker for all Carotis.
#    

# Solution will be to try and detect a 1, del1 or _1 after last name, and if they
# exist look for an identical name but with 2 instead of 1. If it exist run through both
# and append their data together, but add end time from 1 to all of time values from 2.
# otherwise behave as normal.
# In one case there is a 2 part, where name is normal, but there is a del2 which is directly
# connected to the normal name


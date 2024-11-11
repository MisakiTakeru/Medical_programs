#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on Wed Oct  9 2024

@author: Joachim Normann Larsen
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

class stretch_rect:
    
    def __init__(self, rect, allowed_diff):
        self.rect = rect
        self.fig = self.rect.figure
        self.ax = self.rect.axes
        self.allowed_diff = allowed_diff
        self.stretch = False
        
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)

# grabbed need to be one of 3 None, Left or Right to indicate which side that needs to be changed or if neither.
        self.grabbed = 'None'

# will be true when dragging and otherwise False.
        self.dragging = False

# the offset of where the mouse is to the x point dragged.        
        self.offset = None
    
    def on_motion(self, event):
        if event.inaxes != self.ax:
            return
        
        if self.stretch:
            width = self.rect.get_width()
            if self.grabbed == 'Right':
                self.rect.set_width(width + (event.xdata - self.x0))
                self.x0 = event.xdata
# For right side we need to change the width of the rectangle.
            elif self.grabbed == 'Left':
                diff = self.x0 - event.xdata 
                self.rect.set_width(width + diff)
                self.rect.set_x(self.rect.get_x() - diff)
                self.x0 = event.xdata

# For left side we need to change the width of the rectangle and it's x coordinate.
            else:
                return
            self.fig.canvas.draw_idle()
        elif self.dragging:
            new_loc = np.array((event.xdata, event.ydata)) + self.offset
            self.rect.set_xy(new_loc)
            self.fig.canvas.draw_idle()
        return

    def on_click(self, event):
        if event.inaxes != self.ax:
            return
        self.x0, self.y0 = (event.xdata, event.ydata)
        if event.button == 1:
            self.grab_side(event.xdata, event.ydata)
            if self.grabbed == 'Left':
                self.offset = self.rect.get_xy() - np.array([event.xdata, event.ydata])
                self.stretch = True
            elif self.grabbed == 'Right':
                self.offset = self.rect.get_corners()[1] - np.array([event.xdata, event.ydata])
                self.stretch = True
        elif event.button == 3:
            center = self.rect.get_center()
            p = np.array([event.xdata, event.ydata])
            dist = np.linalg.norm(center - p)
            if dist < self.allowed_diff:
                self.dragging = True
                self.offset = self.rect.get_xy() - p
            
        return
    
    def on_release(self, event):
        self.stretch = False
        self.grabbed = 'None'
        self.dragging = False
        return
    
    def grab_side(self, x, y):
        """
        Parameters
        ----------
        x : float
            mouseclick x position.
        y : float
            mouseclick y position.

        Checks if a side of a rectangle is grabbed, and if true saves it onto self.grabbed        
        """
        sides = self.rect.get_corners()
        rect_left_x = sides[0,0]
        rect_right_x = sides[1,0]
        rect_lower_y = sides[0,1]
        rect_upper_y = sides[2,1]

        within_y = ((abs(rect_lower_y - y) < self.allowed_diff) or 
                    (abs(rect_upper_y - y) < self.allowed_diff) 
                    or ((y > rect_lower_y) and (y < rect_upper_y)))

# if there are less than allowed_diff points difference in x values.
        if (abs(rect_left_x  - x) < self.allowed_diff) and within_y:
            self.grabbed = 'Left'
        
        elif (abs(rect_right_x  - x) < self.allowed_diff) and within_y:
            self.grabbed = 'Right'
    def get_area(self):
        x = self.rect.get_x()
        return x, (x + self.rect.get_width())


if __name__ == '__main__':
    fig, ax = plt.subplots()
    rect = Rectangle((0.4,0.4), 0.2, 0.2, fill = False)
    ax.add_patch(rect)
    inter = stretch_rect(rect, 0.01)
    plt.show()
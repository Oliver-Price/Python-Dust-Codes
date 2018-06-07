# -*- coding: utf-8 -*-

import plotly
import plotly.plotly as py
import plotly.figure_factory as ff
import plotly.graph_objs as go

df = [
    dict(Task="STEREO-A HI-1", Start='2007-01-11 00:00:00', Finish='2007-01-12 16:00:00', Resource='Tail & Nucleus in Frame'),
    dict(Task="STEREO-A HI-1", Start='2007-01-12 22:00:00', Finish='2007-01-14 22:00:00', Resource='Tail & Nucleus in Frame'),
    dict(Task="STEREO-A HI-1", Start='2007-01-14 22:00:00', Finish='2007-01-17 10:00:00', Resource='Tail in Frame'),    
    dict(Task="STEREO-A HI-1", Start='2007-01-17 14:00:00', Finish='2007-01-18 18:00:00', Resource='Tail in Frame'),
    dict(Task="STEREO-A HI-1", Start='2007-01-18 22:00:00', Finish='2007-01-19 14:00:00', Resource='Tail in Frame'),
    dict(Task="STEREO-A HI-1", Start='2007-01-19 14:00:00', Finish='2007-01-19 18:20:00', Resource='Test Images'),
    dict(Task="STEREO-A HI-1", Start='2007-01-20 00:00:00', Finish='2007-01-23 22:00:00', Resource='Tail in Frame'),
    #STEREO-A HI-2
    dict(Task="STEREO-A HI-2", Start='2007-01-14 14:00:00', Finish='2007-01-18 18:00:00', Resource='Tail in Frame'),
    dict(Task="STEREO-A HI-2", Start='2007-01-18 22:00:00', Finish='2007-01-19 16:30:00', Resource='Tail in Frame'),
    dict(Task="STEREO-A HI-2", Start='2007-01-20 00:00:00', Finish='2007-01-25 02:00:00', Resource='Tail in Frame'),
    dict(Task="STEREO-A HI-2", Start='2007-01-26 00:00:00', Finish='2007-01-30 00:00:00', Resource='Tail in Frame'),
    dict(Task="STEREO-A HI-2", Start='2007-01-12 19:00:00', Finish='2007-01-12 19:10:00', Resource='Perihelion'),
    #STEREO-B HI-1
    dict(Task="STEREO-B HI-1", Start='2007-01-11 17:59:00', Finish='2007-01-11 19:54:00', Resource='Test Images'), 
    dict(Task="STEREO-B HI-1", Start='2007-01-12 00:01:00', Finish='2007-01-12 01:45:00', Resource='Test Images'),
    dict(Task="STEREO-B HI-1", Start='2007-01-12 10:00:00', Finish='2007-01-12 16:00:00', Resource='Tail & Nucleus in Frame'), 
    dict(Task="STEREO-B HI-1", Start='2007-01-12 20:00:00', Finish='2007-01-12 22:00:00', Resource='Tail & Nucleus in Frame'),
    dict(Task="STEREO-B HI-1", Start='2007-01-13 08:01:00', Finish='2007-01-13 09:37:00', Resource='Test Images'), 
    dict(Task="STEREO-B HI-1", Start='2007-01-14 00:00:00', Finish='2007-01-15 16:00:00', Resource='Tail in Frame'), 
    dict(Task="STEREO-B HI-1", Start='2007-01-16 00:00:00', Finish='2007-01-16 15:00:00', Resource='Tail in Frame'),
    #STEREO-B HI-2  
    dict(Task="STEREO-B HI-2", Start='2007-01-13 18:05:00', Finish='2007-01-13 09:44:00', Resource='Out of Focus'), 
    dict(Task="STEREO-B HI-2", Start='2007-01-14 00:00:00', Finish='2007-01-15 16:00:00', Resource='Out of Focus'), 
    dict(Task="STEREO-B HI-2", Start='2007-01-16 00:00:00', Finish='2007-01-16 21:00:00', Resource='Out of Focus'), 
    dict(Task="STEREO-B HI-2", Start='2007-01-12 19:00:00', Finish='2007-01-12 19:10:01', Resource='Perihelion'),
    #LASCO Clear
    dict(Task="LASCO C3 CLEAR", Start='2007-01-12 20:54:00', Finish='2007-01-15 17:22:00', Resource='Tail & Nucleus in Frame'),
    #Lasco Blue
    dict(Task="LASCO C3 Blue", Start='2007-01-12 20:45:00', Finish='2007-01-15 14:44:00', Resource='Tail & Nucleus in Frame'),
    #Deiries Images
    dict(Task="ESO images from S.Deiries", Start='2007-01-18 00:42:00', Finish='2007-01-23 00:58:00', Resource='Tail & Nucleus in Frame'),
    #Druckmuller Images
    dict(Task="Images from M. Druckmuller", Start='2007-01-19 02:00:00', Finish='2007-01-28 09:30:00', Resource='Tail & Nucleus in Frame')]

colors = {'Test Images': (1, 0.5,0),
          'Tail in Frame': (0,0.8,0),
          'Out of Focus': (0.55,0,0),
          'Tail & Nucleus in Frame': (0,0.5,0.1),
          'Perihelion': (0,0,1)}

fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True, group_tasks=True, showgrid_x=True, showgrid_y=True)
plotly.offline.plot(fig)
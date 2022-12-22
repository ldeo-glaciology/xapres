# -*- coding: utf-8 -*-
"""
Function to plot ApRES data

ApRES_plot('fn',<filename>,'bl',<burstlist>,'cl',<chirplist>,'md',<maxdepth>,'po',<plot options>)
    <plot options> = string containing all or none of 'c','a','r','t'
    'c': plot individual chirps (if absent, plot mean of chirps)
    'a': plot depth profiles, averaged unless "c" included as a plot option
         If no other options, "a" is assumed by default
    'r': plot raw bursts as time series
    't': plot chirp time series, averaged unless "c" included as a plot option
    
    <filename>  filename to plot, including complete path from current folder.
                If no filename given, a file-browser window is displayed.
                
    <burstlist> is either -1 (default) or a Python list of bursts to plot,
                starting at burst 0. -1 means plot all bursts. Each burst is
                plotted in a different figure window
                
    <chirplist> is either -1 (default) or a Python list of chirps within bursts
                to plot, starting at chirp 0. -1 means plot all chirps. Chirps
                from within each burst are either averaged (default) or, if
                plot option "c" is selected, overplotted.
    <maxdepth>  is maximum depth to plot in metres

Created on Sun Oct 18 16:36:29 2020

@author: Keith Nicholls
"""

from .ApRESDefs import *
import matplotlib.pyplot as plt
from tkinter.filedialog import askopenfilename
from os.path import isfile

def ApRES_plot(*arg):
    
     plt.rcParams['figure.max_open_warning'] = 0
     # Unpack argument list
     fn, po, bl, cl, md = GetArguments(arg)

     # Get file containing data to plot
     if len(fn) == 0:
         fn = askopenfilename(initialdir = "/",title = "Select file",
                              filetypes = (("data files","*.dat"),("All files","*.*")))
         
     if isfile(fn):
         fd = ApRESDefs.DataFileObject(fn)
     else:
         print("File", fn, "does not exist\n")
         return()
     
     # Now have number of bursts in file. If bl is -1, that means use them all
     if bl[0] == -1:
         bl = [i for i in range(fd.NoBurstsInFile)]
     else:
         for i in range(len(bl)):
             if bl[i] > fd.NoBurstsInFile - 1:
                 print("Burst ",bl[i]," not in file ",fn )
                 bl.pop(i)
         
     """
     Processing depends on what has been requested. This program supports
     the following plot options:
        "r" Raw (full burst as a time series)
        "t" Average of the chirps selected from the burst, as a time series
        "a" Amplitude plot (return power vs range) from the averaged chirp
        "c" Do not average chirps in each burst
     """
     if po.count("r") > 0:
         for BurstNo in bl:
            b = fd.ExtractBurst(BurstNo)
            plt.figure()
            b.PlotBurst()

     if po.count("t") > 0:
         plt.figure()
         
         for BurstNo in bl:
             b = fd.ExtractBurst(BurstNo)
              
             # Now have number of chirps in burst. If cl[0] is -1, that means use them all
             if cl[0] == -1:
                 cl = [i for i in range(b.Header["NSubBursts"])]
             else:
                 for i in range(len(cl)):
                     if cl[i] > b.Header["NChirps"] - 1:
                         print("Chirp ",cl[i]," not in burst ",BurstNo )
                         cl.pop(i)
                      
             # Check for "c" plot option, meaning do not average chirps
             if po.count("c") == 0:
                 c = b.ExtractChirp(cl)
                 c.PlotChirp()
             else:
                 for i in cl:
                     c = b.ExtractChirp([cl[i]])
                     c.PlotChirp()
              
     if po.count("a") > 0:
         plt.figure()
         
         for BurstNo in bl:
             b = fd.ExtractBurst(BurstNo)
             # Now have number of chirps in burst. If cl[0] is -1, that means use them all
             if cl[0] == -1:
                 cl = [i for i in range(b.Header["NSubBursts"])]
             else:
                 for i in range(len(cl)):
                     if cl[i] > b.Header["NChirps"] - 1:
                         print("Chirp ",cl[i]," not in burst ",BurstNo)
                         cl.pop(i)
              
             if po.count("c") == 0:
                 c = b.ExtractChirp(cl)
                 p = c.FormProfile(b.Header["StartFreq"],b.Header["StopFreq"],2,1)
                 p.PlotProfile(md)
             else:
                 for i in cl:
                     c = b.ExtractChirp([cl[i]])
                     p = c.FormProfile(2e8,4e8,2,1)
                     p.PlotProfile(md)
             
def GetArguments(arg):
    filename = '' # Data filename. Empty string means ask user
    md = 2000    # Default: 2000 m max range to plot 
    po = 'a'     # Default: plot amplitude of mean chirp as a function of range
    bl = [-1]    # Default: plot all bursts in file
    cl = [-1]    # Default: plot all chirps from bursts to be plotted
    
    for ind in range(0,len(arg),2):

        if arg[ind] == "fn" or arg[ind] == "filename":
            filename = arg[ind+1]            

        elif arg[ind] == "md" or arg[ind] == "maxdepth":
            md = arg[ind+1]
            if not isinstance(md,(int, float)):
                print('maxdepth must be positive, real number')
                md = 2000

        elif arg[ind] == "bl" or arg[ind] == "burstlist":
            bl = arg[ind+1]
            if type(bl) != list:
                print('Burst list must be a Python list. Try []')
                bl = [-1]

        elif arg[ind] == "cl" or arg[ind] == "chirplist":
            cl = arg[ind+1]
            if type(cl) != list:
                print('Chirp list must be a Python list. Try []')
                cl = [-1]

        elif arg[ind] == "po" or arg[ind] == "plotoptions":
            po = arg[ind+1]
            if type(po) != str:
                print('plotoptions must be string containing combination of ''c'', ''a'', ''r'', ''t''')
                po = 'a'

        else:
            print("Unrecognised option: ",arg[ind],"\n")

    return filename, po, bl, cl, md

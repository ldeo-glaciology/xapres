Python code
===========

The doc strings in the .py files gives the descriptions of the ApRES classes, and the functionality
of the plotting program.

The classes define the following objects:

FileDescriptorObject
BurstObject
ChirpObject
ProfileObject

The code was written under the Spyder IDE. The version of Python used was (3.8.3), with the Spyder IDE (4.1.4).
The code also uses numPy vs 1.18.5, and matplotlib vs 3.2.2.

To get plots in separate windows, which will allow zoom and pan, "%matplotlib qt" needs to be executed from
the console, otherwise they all appear in the Plots tab of the top-right panel in Spyder.
This is a so-called "magic command", specific to Ipython (I think). And I'm sure someone knows how to set up Spyder
to allow the figure windows to come to the foreground by default, but I am not that person...
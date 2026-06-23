# (M)agnetic (T)weezers Experiment Real-time Bead Tracking

This repository contains the code and resources for real-time bead 
tracking in a horizontal magnetic tweezers experiment. Used with LabView. Call 
Library Function Nodes are used to call the DLL within LabView. Edit 
of previously developed post-processing software to be real-time.

Used by Dr. Sarkar's lab at the Catholic University of America. Part
of REU project.

## Functions
createTracker, trackFrame, getLastFrameTimeMs, destroyTracker

## Dependencies
- LabView 2015 64x or later
- OpenCV 4.12

## Common Problems
If while building the DLL you get an error about not being able to find the OpenCV libraries, 
you can fix this by pasting the opencv_world(4120).dll and opencv_world(4120)d.dll files into 
the same directory as the DLL you are building.
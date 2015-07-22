#!/usr/bin/env python
# 
"""
These functions, when given a magnitude mag between cmin and cmax, return
a colour tuple (red, green, blue). Light blue is cold (low magnitude)
and yellow is hot (high magnitude).

"""
import math

def floatRgb(mag, cmin, cmax):
       """
       Return a tuple of floats between 0 and 1 for the red, green and
       blue amplitudes.

       try:
              # normalize to [0,1]
              x = float(mag-cmin)/float(cmax-cmin)
       except:
              # cmax = cmin
              x = 0.5
       blue = min((max((4*(0.75-x), 0.)), 1.))
       red  = min((max((4*(x-0.25), 0.)), 1.))
       green= min((max((4*math.fabs(x-0.5)-1., 0.)), 1.))
       """
       diff = float(cmax)-float(cmin)
       scale = 255/diff
       num = float(mag)*scale

       red = num
       blue = 255 - num
       green = 5


       return (red, green, blue)

def strRgb(mag, cmin, cmax):
       """
       Return a tuple of strings to be used in Tk plots.
       """

       red, green, blue = floatRgb(mag, cmin, cmax)       
       return "#%02x%02x%02x" % (red*255, green*255, blue*255)

def rgb(mag, cmin, cmax):
       """
       Return a tuple of integers to be used in AWT/Java plots.
       """

       red, green, blue = floatRgb(mag, cmin, cmax)
       return (int(red), int(green), int(blue))

def htmlRgb(mag, cmin, cmax):
       """
       Return a tuple of strings to be used in HTML documents.
       """
       return "#%02x%02x%02x"%rgb(mag, cmin, cmax)

#Create empty list
tempFactors = []
#Run through file and grab temperature factors
file = open('/Users/aaronrenberg/BeakerChord/Main Beaker Notebook/1OTH.pdb','r')
for line in file:
	if line[0:4] == 'ATOM':
		if line[13:15] == 'CA':
			tempFactor = line[61:66]
			tempFactors.append(tempFactor)
print min(tempFactors)
print max(tempFactors)
print rgb(tempFactors[4],min(tempFactors),max(tempFactors))
print 'DONE'
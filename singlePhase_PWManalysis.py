#!/usr/bin/env python3
# file: 1ph_PWManalysis.py
# Author    : Najath Abdul Azeez
# Copyright : http://opelex.net
# License   : See LICENSE file
"""
Analysis of various single-phase PWM technqiues
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from matplotlib.widgets import Slider, RadioButtons

# Function for initializing global variables
def initGlobalVars():
    # Initial fundamental freq, amplitude and frequency modulation
    global f, ma, mf
    f  = 1.0
    ma = 0.9
    mf = 20.0
    # Define initial PWM technique
    global PWMtype
    PWMtype = 'Unipolar'
    # Number of samples, and equivalent time array
    global N, t, zeros
    samples = 1e6
    N = int(samples * f)
    t = np.linspace(0, 1, N)
    zeros = np.zeros_like(t)


##############################################################################
# Functions for creating and initializing plots                              #
##############################################################################
# Function for creating figures and required subplots
def createFig():
    global fig
    fig = plt.figure()
    # Create 3 subplots
    global axFFT, axST, axPWM
    axFFT  = plt.subplot2grid((5, 5), (0, 0), colspan=3, rowspan=5)
    axST   = plt.subplot2grid((5, 5), (0, 3), colspan=2, rowspan=2)
    axPWM  = plt.subplot2grid((5, 5), (2, 3), colspan=2, rowspan=3)


# Plot all waveforms and get the handles
def initPlots():
    # All waveforms are initialized to Zero line
    # sine and triangle waveforms
    global sineW, sine2W, trngW
    axST.plot(t, zeros)
    sineW,  = axST.plot(t, zeros)
    trngW,  = axST.plot(t, zeros)
    sine2W, = axST.plot(t, zeros)
    # PWM waveforms
    global pwmW, pwm1W, pwm2W
    pwm0W,  = axPWM.plot(t, zeros)
    pwm1W,  = axPWM.plot(t, zeros)
    pwmW,   = axPWM.plot(t, zeros)
    pwm2W,  = axPWM.plot(t, zeros)
    # Plot cursors
    global curST, curPWM
    curST,   = axST.plot([0, 0], [-1.3, 1.3], 'k')
    curPWM,  = axPWM.plot([0, 0], [-3.3, 3.3], 'k')


# Function for plot beautification
def plotChromes():
    axPWM.set_yticks(range(-3, 4, 1))
    axPWM.set_yticklabels([])
    axPWM.grid(True, axis='y')


# Function for plotting FFT of PWM output
def plotFFT():
    # Range of frequencies to plot
    h_range = range(100)
    # Calculate FFT
    Y = np.fft.fft(pwm) / N
    # Calculate single-sided amplitude spectrum
    Y1 = 2 * abs(Y[h_range])
    # Correct doubling in DC amplitude
    Y1[0] = Y1[0] / 2
    # Clear, then plot the FFT
    axFFT.cla()
    axFFT.stem(h_range, Y1)
    # Fix Y-axis limits from -0.1 to 1.3
    axFFT.set_ylim(-0.1, 1.3)


##############################################################################
# Functions for updating global variables and plots                          #
##############################################################################
# Function for updating modulating waves
def updateMwave():
    global sine
    sine = ma * np.sin(2 * np.pi * f * t)
    sineW.set_ydata(sine)


# Function for updating carrier waves
def updateCwave():
    global triangle
    triangle = signal.sawtooth(2 * np.pi * mf * f * t + (np.pi / 2), 0.5)
    trngW.set_ydata(triangle)


# Function for updating pwm
def updatePWM():
    global pwm
    pwmH1 = ( sine >= triangle) * 1
    pwmL1 = ( sine <  triangle) * 1
    pwmH2 = (-sine >= triangle) * 1
    # Leg A pwm
    pwmA = pwmH1
    # Leg B pwm
    pwmB = {
        'SingleLeg'  : 0,
        'HalfBridge' : 0.5,
        'Bipolar'    : pwmL1,
        'Unipolar'   : pwmH2,
    }[PWMtype]
    # Output pwm
    pwm  = pwmA - pwmB
    # Modulating wave for second leg
    if PWMtype == 'Unipolar':
        m2 = -1
    else:
        m2 = 0
    # Update pwm plots
    pwmW.set_ydata(pwm - 2)
    pwm1W.set_ydata(pwmA + 2)
    pwm2W.set_ydata(pwmB)
    sine2W.set_ydata(sine * m2)
    # plot FFT of updates pwm waves
    plotFFT()

##############################################################################
# HMI initialization and call back functions                                 #
##############################################################################
def initHMI():
    global sCr, sMa, sMf, rPWM
    # Make space in fig bottom, for placing the sliders and radio button
    fig.subplots_adjust(bottom=0.20)
    # Define colors to be used
    axbgcolor = 'lightgoldenrodyellow'

    # Create slider for cursor
    axCr = plt.axes([0.05, 0.070, 0.63, 0.02], facecolor=axbgcolor)
    sCr  = Slider(axCr, 'time', 0, 1, valinit=0)
    # Connect function for updating change in "Modulation index"
    sCr.on_changed(updateCursor)

    # Create slider for modulation index (ma)
    axMa = plt.axes([0.05, 0.025, 0.25, 0.02], facecolor=axbgcolor)
    sMa  = Slider(axMa, 'ma', 0.1, 1, valinit=ma, valstep = 0.1)
    # Connect function for updating change in "Modulation index"
    sMa.on_changed(updateMa)

    # Create slider for freq modulation (mf)
    axMf = plt.axes([0.43, 0.025, 0.25, 0.02], facecolor=axbgcolor)
    sMf  = Slider(axMf, 'mf', 1, 25, valinit=mf, valstep = 1)
    # Connect function for updating change in mf
    sMf.on_changed(updateMf)

    # Create radio button for selecting PWM technique
    axPWM = plt.axes([0.75, 0.000, 0.20, 0.15], facecolor=axbgcolor)
    rPWM  = RadioButtons(axPWM, ('SingleLeg', 'HalfBridge', 'Bipolar', 'Unipolar'), active=3)
    # Connect function for updating change in PWM technique
    rPWM.on_clicked(updatePWMtype)


# Function for updating cursor
def updateCursor(tNew):
    curST.set_xdata([tNew, tNew])
    curPWM.set_xdata([tNew, tNew])

# Function for updating ma
def updateMa(maNew):
    global ma
    ma = maNew
    updateMwave()
    updatePWM()

# Function for updating mf
def updateMf(mfNew):
    global mf
    mf = mfNew
    updateCwave()
    updatePWM()

# Function for updating PWM
def updatePWMtype(PWMnew):
    global PWMtype
    PWMtype = PWMnew
    updatePWM()

##############################################################################
# Master initialization function                                             #
##############################################################################
def init():
    initGlobalVars()
    createFig()
    initPlots()
    plotChromes()
    updateMwave()
    updateCwave()
    updatePWM()
    # Adjust plot sizes to avoid overlaps
    plt.tight_layout()


##############################################################################
# Call initialization function and show the plot                             #
##############################################################################
init()
initHMI()
plt.show()


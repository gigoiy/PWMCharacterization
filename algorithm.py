#Import libraries for data simulation and signal characterization
import matplotlib.pyplot as plt
import numpy as np 
from scipy import signal

T0 = float(5.59e-6) #Calculated fundamental period
omega = float(1124004.527) #Calculated fundamental angular frequency
n = np.arange(1.00, 7.00, 1.00) #List with values of n (each harmonic)
dutycycle = np.arange(0.00, 1.00, 0.01) #List with each percent of duty cycle
spectrum = [] #Frequency spectrum x-axis values
magnitude = [] #Magnitude values post-transform used in plotting y-axis for magnitude plot
phase = [] #Phase values post-transform used in plotting y-axis for phase plot
a0 = [] #Fourier constants throughout the frequency spectrum
aharmonics = [] #List containing one list per harmonic that contains Fourier coefficients throughout the frequency spectrum
jterms = [] #Calculated values from the imaginay term of the transformed equation
realterms = [] #Calculated values from the real term of the transformed equation
totals = [] #Final values calculated from the transformed equation
t = [] #Time base on the x-axis

#Improper definite integral setup
intervalA1 = 0 #Lower limit of the first integral
intervalB1 = float(T0) #Upper limit of the first integral
intervalA2 = intervalB1 #Lower limit of the second integral 
intervalB2 = float(T0) #Upper limit of the second integral

for j in range(len(n)): #Increment through each harmonic
    a = [] #Initialize a new Fourier coefficient list at the start of each loop
    for i in range(len(dutycycle)): #Increment through the frequency spectrum of the signal
        intervalB1 = T0*dutycycle[i] #Increment the upper limit of the first integral
        intervalB2 = T0*(1.00-dutycycle[i]) #Calculate and increment the upper limit of the second integral
        #Calculate the Fourier coefficients throughout the frequency spectrum and append it to the list for the current harmonic
        a.append(float((2/T0)*(omega*np.sin(n[j]*omega*intervalB1))-(omega*np.sin(n[j]*omega*intervalA1))-(omega*np.sin(n[j]*omega*intervalB2))+(omega*np.sin(n[j]*omega*intervalA2))))
        a[i] = a[i]*(10.00**-9.00) #Convert the exponentials to integers
    aharmonics.append(a) #Append the list of Fourier coefficients of the current harmonic to the finalized list 

for i in range(len(dutycycle)): #Increment through the frequency spectrum
    intervalB1 = T0*dutycycle[i] #Increment the upper limit of the first integral
    intervalB2 = T0*(1.00-dutycycle[i]) #Calculate and increment the upper limit of the second integral
    a0.append(float((1/T0)*(intervalB1-intervalA1-intervalB2+intervalA2))) #Calculate the Fourier constant throughout the frequency spectrum
    t.append(intervalB2) #Creates the time base
    spectrum.append(float(1.00/intervalB2)) #Calculates and creates the frequency spectrum values
    realterms.append(float((3.00*np.sin(omega*t[i]))/omega)) #Calculates the values of the real terms
    jterms.append(float((1.00 + np.cos(omega*t[i]))/omega)) #Calculates the values of the imaginary terms
    totals.append((realterms[i]**2.00)+(jterms[i]**2.00)) #Calculates the sum of both terms, getting the final value
    phase.append(np.arctan(jterms[i]/realterms[i])) #Calculates the phase angle of the resulting vectors
    spectrum[i] = spectrum[i]*(10.00**-6.00) #Convert the exponentials to integers

magnitude = np.sqrt(totals) #Calculate magnitudes of the vectors
for i in range(len(magnitude)): #Increment through the frequency spectrum
    magnitude[i] = magnitude[i]*(10.00**6.00) #Convert the exponentials to integers

#Create a plot showing the harmonics of the signal
fig, harmonics = plt.subplots()

harmonics.bar(x=spectrum, height=magnitude, color='maroon', width=0.1)
plt.xlabel('Frequency (MHz)')
plt.ylabel('Magnitude (udB)')
plt.title('Harmonics of a 179 kHz PWM Square Wave')
plt.grid()
plt.savefig('harmonics.png')

#Create a plot showing the phase of the signal
fig, phasedata = plt.subplots()

phasedata.plot(spectrum, phase, color='green')
phasedata.set(xlabel='Frequency (MHz)', ylabel='Phase (Radians)', title='Phase of a 179 kHz PWM Square Wave on the Frequency Spectrum')
plt.grid()
plt.savefig('phase.png')

#Create a plot showing the Fourier coefficients throughout the frequency spectrum
fig, coefficients = plt.subplots()

coefficients.plot(spectrum, a0, color='purple', label='A0')
coefficients.plot(spectrum, aharmonics[0], color='blue', label='Fundamental')
coefficients.plot(spectrum, aharmonics[1], color='yellow', label='First Harmonic')
coefficients.plot(spectrum, aharmonics[2], color='red', label='Second Harmonic')
coefficients.plot(spectrum, aharmonics[3], color='orange', label='Third Harmonic')
coefficients.plot(spectrum, aharmonics[4], color='brown', label='Fourth Harmonic')
coefficients.plot(spectrum, aharmonics[5], color='gray', label='Fifth Harmonic')
coefficients.set(xlabel='Frequency (MHz)', ylabel='Amplitude (GV)', title='Fourier Coefficients Throughout the Frequency Spectrum')
coefficients.legend()
plt.grid()
plt.savefig('coefficients.png')

#Create a plot of the initial PWM square wave that we are characterizing
fig, squarewave = plt.subplots()
t = np.linspace(0, 3*T0, 10000, endpoint=False)
frequency = 179000
squarewave.plot(t*1e6, signal.square(2*np.pi*t*frequency, duty=.25))
squarewave.set(xlabel='Time (us)', ylabel='Amplitude (V)', title='PWM Square Wave at 179 kHz and 25% Duty Cycle')
plt.grid()
plt.savefig('squarewave.png')

plt.show()
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 18:02:34 2019

@author: rylab
"""

import nidaqmx
import time
import numpy as np
import matplotlib.pyplot as plt
from Chameleon import Chameleon
from wavelength_to_rgb import wavelength_to_rgb

class Calibration:
	def __init__(self, THORpowerMeterRange, verbose=True, wavelengths=np.arange(800,1080+1,10)):
		self.verbose = verbose
		self.THORpowerMeterRange = THORpowerMeterRange
		
		#initialize the connection to the pockels cell and the power meter
		self.setupPockelsCell(channels='PXI1Slot4/ao0:1')
		self.setupPowerMeter(channels='PXI1Slot4/ai13')
		
		#initialize the connection to the chameleon
		self.laser = Chameleon(verbose=self.verbose, com_port='COM1')
		
		#set pockels cell bias, cmd and shutter to safe positions
		self.safeSystem()
		
		self.wavelengths = wavelengths;
		self.wavelengths = np.rint(self.wavelengths) #round wavelengths to integers
		self.wavelengths = self.wavelengths.astype(int) #cast the wavelengths array as integers
		self.wait_time_s = 4 # (s) time to wait after adjusting the bias values (time constant of S175C is approx 1.5-2 s)
		self.bias_params = np.array([-0.005571, 5.522])
		self.bias_estimation = lambda wavelength: self.bias_params[0]*wavelength+self.bias_params[1]
		
		self.reprate = 80e6 # Hz
		
		self.newWavelengths = self.wavelengths.copy()
		np.random.shuffle(self.newWavelengths)
		self.newBiasValues = np.array([])
		
		
		self.record_bias = np.array([])
		self.record_command = np.array([])
		self.record_wavelength = np.array([])
		self.record_power = np.array([])
		self.record_photons_per_pulse = np.array([])
		
		
		
	def setupPockelsCell(self, channels):
		self.pockels_cell = nidaqmx.Task()
		self.pockels_cell.ao_channels.add_ao_voltage_chan(channels)
		self.pockels_cell.timing.cfg_samp_clk_timing(rate=1000,sample_mode=nidaqmx.constants.AcquisitionType.FINITE,samps_per_chan=1)
		
		
		
	def setupPowerMeter(self, channels):
		self.power_meter = nidaqmx.Task()
		self.power_meter.ai_channels.add_ai_voltage_chan(channels)
		self.power_meter.timing.cfg_samp_clk_timing(rate=1000,sample_mode=nidaqmx.constants.AcquisitionType.FINITE,samps_per_chan=100)



	def findBiasMin(self, wavelength):
		#set the wavelength
		self.laser.setWavelengthBlocking(wavelength)

		#find the previous estimation of the bias minimum
		previous_bias = self.bias_estimation(wavelength)
		bias_to_test = np.arange(-2,8,.5)
		bias_to_test_shuffled = bias_to_test.copy()
		np.random.shuffle(bias_to_test_shuffled)
		if(self.verbose):
			print('Previous Bias Estimation: ', previous_bias)
			
		#take a baseline reading
		baseline_value = self.setPockelsCellValuesAndRecordPower(0, 0)
		if(self.verbose):
			print('Baseline Value (Shutter Closed): ', baseline_value)
		
		intensities = np.array([])
		self.laser.openShutterBlocking()
		for bias in bias_to_test_shuffled:
			AI_value = self.setPockelsCellValuesAndRecordPower(bias, 0)
			
			power = AI_value-baseline_value
			
			intensities = np.append(intensities, power)
			
			#record all values
			self.record_bias = np.append(self.record_bias, bias)
			self.record_command = np.append(self.record_command, 0)
			self.record_wavelength = np.append(self.record_wavelength, wavelength)
			self.record_power = np.append(self.record_power, power)
			self.record_photons_per_pulse = np.append(self.record_photons_per_pulse, self.convertPowerToPhotonsPerPulse(power, wavelength))
			
			if(self.verbose):
				print(bias, ',', power)
			
		self.safeSystem()
		
		#fit a function (3rd order) to the data and find the local minimum
		p = np.polyfit(bias_to_test_shuffled, intensities, 3)
		fitline = lambda x,p: p[0]*np.power(x,3) + p[1]*np.power(x,2) + p[2]*x + p[3]
		poly3diff = lambda p: np.array([p[0]*3, p[1]*2, p[2]])
		quadzeros = lambda p: np.array([(-p[1] + np.sqrt(p[1]**2-4*p[0]*p[2]))/(2*p[0]), (-p[1] - np.sqrt(p[1]**2-4*p[0]*p[2]))/(2*p[0])])
		poly2diff = lambda p: np.array([p[0]*2, p[1]])
		poly2 = lambda x,p: p[0]*x + p[1]
		
		#find the local minimum of the 3rd order polynomial
		zeros = quadzeros(poly3diff(p))
		derivative = poly2(zeros, poly2diff(p))
		local_minimum = zeros[np.argmax(derivative)]
		print(local_minimum)
		
		plt.plot(bias_to_test_shuffled, intensities, marker='o', markeredgecolor=wavelength_to_rgb(wavelength/2), markerfacecolor=wavelength_to_rgb(wavelength/2), linestyle='None')
		plt.plot(bias_to_test, fitline(bias_to_test, p), color=wavelength_to_rgb(wavelength/2), linestyle='dashed')
		plt.show()
		plt.pause(.1)
		
		return local_minimum
	
	
	
	def convertPMVoltageToPower(self, voltage):
		return voltage*self.THORpowerMeterRange/2
		
	
	
	def convertPowerToPhotonsPerPulse(self, power_in_mW, wavelength_in_nm):
		power_in_W = power_in_mW/1000
		energy_per_pulse = power_in_W/self.reprate
		wavelength_in_um = wavelength_in_nm/1000
		joules_per_ev = 1.602176634e-19 # joules per electron volt
		energy_per_photon = (1.2398/wavelength_in_um)*joules_per_ev
		photons_per_pulse = energy_per_pulse/energy_per_photon
		return photons_per_pulse
	
	
	
	def calibrateBiasAllWavelengths(self):
		if(self.newBiasValues.size > 0):
			return
		
		plt.figure() #create a new figure
		for wavelength in self.newWavelengths:
			self.newBiasValues = np.append(self.newBiasValues, self.findBiasMin(wavelength))

		# fit a line to the values 
		self.bias_params = np.polyfit(self.newWavelengths, self.newBiasValues, 1)
		
		
		
	def calibrateCommandAllWavelengths(self):
		plt.figure() #create a new figure
		for wavelength in self.newWavelengths:
			self.commandCalibrationAtWavelength(wavelength)
		
		
		
	def setPockelsCellValuesAndRecordPower(self, bias, cmd):
		#set the bias command (bias) and set the control voltage to 0
		self.pockels_cell.write(np.array([[bias], [cmd]]),auto_start=True)
			
		#wait for the value to take and the meter to settle
		time.sleep(.01)
		self.pockels_cell.stop()
		time.sleep(self.wait_time_s)
			
		#read in the value from the power meter
		AI_value = self.convertPMVoltageToPower(np.mean(self.power_meter.read(200)))
		time.sleep(.2)
		
		return AI_value
			
	
	
	def commandCalibrationAtWavelength(self, wavelength):
		#set the wavelength
		self.laser.setWavelengthBlocking(wavelength)
		
		#find the previous estimation of the bias minimum
		bias = self.bias_estimation(wavelength)
		cmd_to_test = np.arange(0, 1.01, .1)
		cmd_to_test_shuffled = cmd_to_test.copy()
		#inp.random.shuffle(cmd_to_test_shuffled)
		
		intensities = np.array([])
		self.laser.openShutterBlocking()
		for cmd in cmd_to_test_shuffled:
			AI_value = self.setPockelsCellValuesAndRecordPower(bias, cmd)
			
			intensities = np.append(intensities, AI_value)
			
			if(self.verbose):
				print(cmd, ',', AI_value)
				
		self.safeSystem()
		
		plt.plot(cmd_to_test_shuffled, intensities, marker='o', markeredgecolor=wavelength_to_rgb(wavelength/2), markerfacecolor=wavelength_to_rgb(wavelength/2), linestyle='None')
		plt.show()
		plt.pause(.1)
		
		return
	
	
	
	def safeSystem(self):
		self.laser.closeShutterBlocking()
		self.pockels_cell.write(np.array([[0], [0]]),auto_start=True)
		time.sleep(.01)
		self.pockels_cell.stop()
		
		
		
	def close(self):
		self.laser.closeShutterBlocking()
		self.laser.close()
		self.pockels_cell.stop()
		self.pockels_cell.close()
		self.power_meter.stop()
		self.power_meter.close()
		
		if(self.verbose):
			print('Calibration Finished')

calib = Calibration(THORpowerMeterRange=5500)
plt.figure()
calib.calibrateBiasAllWavelengths()
calib.close()
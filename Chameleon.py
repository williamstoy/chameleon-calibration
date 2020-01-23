#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 23:44:38 2019

@author: williamstoy
"""

import serial
import time
import re

class Chameleon:
	##########################################################################
	## INITIALIZATION FUNCTIONS
	##########################################################################
	def __init__(self, com_port='COM1', verbose=False):
		self.verbose=verbose
		self.ser = serial.Serial(com_port, 19200, timeout=1)  # open serial port
		# initialize the connection by writing a carriage return
		self.sendCommand(b'')
		time.sleep(0.05)
		
		self.current_wavelength = 1000
		
		
	##########################################################################
	## COMMAND / RESPONSE PRIMITIVE FUNCTIONS
	##########################################################################	
	def sendCommand(self, cmd):
		self.ser.write(cmd + b'\r')
		
	def getResponse(self, cmd):
		response = self.ser.read(50)
		#format the command so that it can be taken out with a regular expression
		cmd = re.sub(b'[\?]', b'\?', cmd)
		cmd = re.sub(b'[\=]', b'\=', cmd)
		#strip out \r and \n, CHAMELEON> and the command
		response = re.sub(b'[\\r|\\n|(CHAMELEON\>)|(' + cmd + b')|\s]', b'', response)
		return response
	
	def sendCmdGetResponse(self, cmd):
		self.sendCommand(cmd)
		time.sleep(0.05)
		return self.getResponse(cmd)
	
	
	##########################################################################
	## WAVELENGTH FUNCTIONS
	##########################################################################
	def setWavelength(self, n):
		self.current_wavelength = n
		if(self.verbose):
			print('Tuning the laser to ', n, 'nm')
		self.sendCmdGetResponse(b'vw=' + bytes(str(n), 'utf-8'))
			
	def queryTunedStatus(self):
		return self.sendCmdGetResponse(b'?ts')
	
	def setWavelengthBlocking(self, n):
		self.setWavelength(n)
		while (self.queryTunedStatus() != b'0'):
			if(self.verbose):
				print('Waiting for tuning to finish...')
			time.sleep(0.25)
		if(self.verbose):
			print('Laser is Tuned')
	
	##########################################################################
	## SHUTTER FUNCTIONS
	##########################################################################
	def openShutter(self):
		if(self.verbose):
			print('Opening the Shutter')
		self.sendCmdGetResponse(b's=1')
			
	def closeShutter(self):
		if(self.verbose):
			print('Closing the Shutter')
		self.sendCmdGetResponse(b's=0')
			
	def queryShutterStatus(self):
		return self.sendCmdGetResponse(b'?s')
			
	def openShutterBlocking(self):
		self.openShutter()
		while (self.queryShutterStatus() != b'1'):
			if(self.verbose):
				print('Waiting for the shutter to open..')
			time.sleep(0.25)
		if(self.verbose):
			print('Shutter is OPEN')
			
	def closeShutterBlocking(self):
		self.closeShutter()
		while (self.queryShutterStatus() != b'0'):
			if(self.verbose):
				print('Waiting for the shutter to close..')
			time.sleep(0.25)
		if(self.verbose):
			print('Shutter is CLOSED')
			
	def queryRelativeHumidity(self):
		if(self.verbose):
			print('Getting Relative Humidity')
		self.sendCmdGetResponse(b'?rh')			
	
	##########################################################################
	## SERIAL FUNCTIONS
	##########################################################################
	def close(self):
		self.ser.close()
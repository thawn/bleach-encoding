# VisiView Macro used to bleach a microtubule with several spatial frequencies.

import time

line_execution_time = 0.005
frap_delay = 0.001
microtubule_velocity = 600.0
bleach_length = 12000.0
wavecounts = [8]#8, 10, 12, 14] #spatial frequencies in terms of number 
				#of bleaching events

#8 corresponds to 1500 nm
#10 corresponds to 1200 nm
#12 corresponds to 1000 nm
#14 corresponds to 857 nm

timepoints=[0]

#In this loop, you convert your space-points when you want to bleach
#in time points you want to beach (i.e. from spatial frequency to 
#time period). The wavecount array contains the number of events. For that
#specific bleach length, the number of events corresponds to the number of times
#you need to bleach in order to draw that particular spatial frequency (ofc the 
#period takes into account the distance you need between the points where you
#bleach). 

for wc in wavecounts:
	period = bleach_length / wc / microtubule_velocity
	print(period)
	print(wc)
	for n in range(wc+1):
		timepoints.append(n*period)
		print(n)
	print(timepoints)

#timepoints = list(dict.fromkeys(timepoints)) #Rachele, there were repetitions in the list. 
timepoints.sort() 
print(timepoints)
print(len(timepoints))

#timeline=[0] #Rachele, because now my list does not start with two zeros anymore 
timeline=[]  #how it was before 

#Devo capire questa parte 

for n,t in enumerate(timepoints):
	if n > 0:
	   timeline.append(timepoints[n] - timepoints[n-1] - line_execution_time) #- 51* frap_delay)
	   #calculates the distance between time points, so it knows after what period of time 
	   #it has to bleach starting from zero 
	   
	   #you subtract line_execution_time and frap_delay because in the time it takes to go
	   #from one time point to the next time point (i.e. the time you need to wait before
	   #bleaching the next timepoint) you have to take into account that the system has an 
	   #intrinsic delay (which means your period of time cannot be shorter than this time 
	   #but also that the system is waiting for that time already no matter what, you don't 
	   #have to count it twice) 
	   
print(timeline)
print(len(timeline))
for t in timeline:
	if t > 0:
		time.sleep(t)
	VV.Frap.Start() # VV is the VisiView API object
	VV.Macro.Control.WaitFor('VV.Frap.IsRunning','==',False)
	time.sleep(frap_delay) # important line, if this is not done, the bleaching intensity of following steps will be wrong
#	VV.Device.SendStringWithoutAnswer('Obis-488_Laser488','SOUR:POW:LEV:IMM:AMPL 0.0004')
#	time.sleep(50*frap_delay)
	

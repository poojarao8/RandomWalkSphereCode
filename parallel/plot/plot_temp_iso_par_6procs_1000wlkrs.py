import numpy as np
import pickle as pk
import matplotlib.pyplot as plt
import glob
import sys
import math
from datetime import datetime

# data from abaqus
time_abq = np.linspace(0,0.2,4001)
babq = np.loadtxt('abaqus_iso_backTemp.txt', usecols=(0,))

# this time is for pulling certain data
ts_abq = np.linspace(0,4001,200)
abq_back_temp = np.zeros(200)

for i in range (0, len(ts_abq)):
  abq_back_temp[i] = babq[min(4000,int(ts_abq[i]))]

#plt.plot(time_abq, abq_back_temp, '-', linewidth=2.5, color='red', label='Back Panel Abaqus')
fileRW = 'temp_iso_par_3procs_1000wlkrs.txt'
timeRW = np.loadtxt(fileRW, usecols=(0,), skiprows=12)

tempRW = np.loadtxt(fileRW, usecols=(1,), skiprows=12)

#plt.plot(timeRW, tempRW, '-', linewidth=2.5, color='blue', label='Back Panel Walker')


#wlkrRW = np.genfromtxt('', delimiter=" ", skip_header=12, comments="?")

#str2date = lambda x: datetime.strptime(x.decode("utf-8"),'%Y-%m-%d.%H:%M:%S')
#dataRW = np.genfromtxt(fileRW, delimiter=" ",converters ={6: str2date},skip_header=12,comments="?")

#print dataRW[15][6].tim

#plt.plot(timeRW, wlkrRW[:, 4])
#plt.plot(timeRW, wlkrRW[:, 4])

plt.plot(time_abq, babq)
plt.plot(timeRW, tempRW)

plt.legend(loc=4)
plt.xlabel('Time (s)', fontsize=14, fontweight='bold')
#plt.ylabel('# Walkers', fontsize=14, fontweight='bold')
plt.ylabel('Temperature (K)', fontsize=14, fontweight='bold')
#plt.ylabel('# Total CPU time', fontsize=14, fontweight='bold')
plt.show()
#plt.show(block=False)
#plt.savefig('temp_iso.png')




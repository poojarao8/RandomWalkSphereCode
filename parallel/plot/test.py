import numpy as np
import pickle as pk
import matplotlib.pyplot as plt
import glob
import sys
import math

dataRW = np.genfromtxt('temp_iso_par.txt', delimiter=" ",
skip_header=12, comments="?")

print dataRW[:, 6]

a = "2017-07-02.23:49:56"

b = "2017-07-02.23:50:22"

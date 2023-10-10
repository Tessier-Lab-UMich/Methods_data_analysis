# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 16:32:36 2019

@author: makow
"""

import pandas as pd
import numpy as np


# defines the filepath of the data, csv or txt format
# filepath eaxmple = "C:\\Users\\makow\\Documents\\Research\\Data Analysis\\10.1.20_AC-SINS.csv"
filepath = '/Users/illogical/Desktop/2022 Research /CS-SINS/6.15/6.15 CS-SINS.csv'

# number of technical replicates performed for each sample, results will be averaged for plasmon wavelength
tech_rep = 2 

# imports data, first column of data is list of wavelengths 450-650, column header = 'Wavelength'
data = pd.read_csv(filepath, sep= ',', header = 0, index_col = 'Wavelength')
# number of antibodies run in the data being processed
mabs = np.arange(0,(len(data.columns)/tech_rep))

# defined variable to add cure fits to
final_nm_max = []
max_au = []
final_nm_ave = []

# takes each column and fits a curve across the 40 points surrounding the max, finds derivative, finds max wavelength by setting derivative = 0
for column in data:
    data_flat = list(data[column])
    abs_max = max(data_flat)
    max_au.append(abs_max)
    idx_max = np.argmax(data_flat)
    points_to_fit = data_flat[idx_max-20:idx_max+20]
    nm_to_fit = list(data.index[idx_max-20:idx_max+20])
    coefs = np.polyfit(nm_to_fit, points_to_fit, 2)
    deriv_coefs = np.polyder(coefs)
    final_nm_max.append((-1*deriv_coefs[1])/deriv_coefs[0])

# averages technical replicates
t = 0
for i in mabs:
    three_cols = final_nm_max[t:t+tech_rep]
    final_nm_ave.append(sum(three_cols)/len(three_cols))
    t = t + tech_rep

results = pd.DataFrame(final_nm_ave, index = mabs, columns = ['Plasmon_Wavelength'])



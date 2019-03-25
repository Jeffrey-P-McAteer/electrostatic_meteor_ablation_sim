#!/usr/bin/python

import readdomain_ftarray as rd
import matimage_plot as plt

arr = rd.readdomain('testArray.bin', (8,32), 4, 'normal', 'ft')
arr1 = rd.readdomain('testArrayFT.bin', (8,32), 4, 'ft', 'ft')

arrs = (arr, arr1)
plt.image_plot(arrs, (2,1))

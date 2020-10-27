import numpy as np
import csv
import pandas

def leff(filter):
	lam_eff = np.sum(filter[:,0] * filter[:,1]) / np.sum(filter[:,1])
	return lam_eff

def ReadFilters(File,convert_to_angstrom):
	#Imports all filter data
	print("Importing Filter Data")
	rowc = 0
	N_filters = 0
	FilterFile = pandas.read_table(File,names=['filters','type'],header=None,dtype={'filters': np.str, 'type': np.str},delim_whitespace=True)
	N_filters = FilterFile.shape[0]
	FilterName = np.array(FilterFile['filters'])
	FilterType = np.array(FilterFile['type'])
	FilterLength = np.int_(np.zeros(N_filters))
	max_FilterLength = 0
	#print('There are ', N_filters, ' filters in the filter file.')
	#print('Finding longest filter...')
	for f in range(0,N_filters):
		filename = FilterName[f]
		data = pandas.read_table(filename,names=['wavelength','transmission'],header=None,skiprows=1,dtype={'wavelength': np.float, 'transmission': np.float},delim_whitespace=True)
		FilterLength[f] = data.shape[0]
	max_len_filters = max(FilterLength)	
	#print('Longest filter has ', max_len_filters, ' wavelength bins')
	Filters = np.empty((N_filters,max_len_filters,2))#Hard coded filter bin number max, max must be 1 < x < 30
	Filters.fill(0)
	FilterWaveEff = np.empty(N_filters)
	for f in range(0,N_filters):
		filename = FilterName[f]
		#print('opening filter '+filename+' of type '+FilterType[f])
		data = pandas.read_table(filename,names=['wavelength','transmission'],header=None,skiprows=1,dtype={'wavelength': np.float, 'transmission': np.float},delim_whitespace=True)
		for i in range(0,FilterLength[f]): 
			Filters[f,i,0] = np.float_(data['wavelength'][i]) * convert_to_angstrom[f]
			Filters[f,i,1] = np.float_(data['transmission'][i])
		FilterWaveEff[f] = leff(Filters[f,:,:])
	return FilterName,Filters,FilterType,FilterLength,N_filters,FilterWaveEff
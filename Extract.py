####################################################################################################################
### title: "Extract"
### author: "Kzzhu"
### time: "2022.06.04"
####################################################################################################################
python -m pip install pyradiomics
conda install -c radiomics pyradiomics

from __future__ import print_function
import sys
import os
import logging
import six
from radiomics import featureextractor, getFeatureClasses
import radiomics
import pandas as pd

dataDir = "E:/CT"
folderList = pd.read_csv('E:/list.csv', sep = '\t')
extractor = featureextractor.RadiomicsFeatureExtractor()
featureClasses = getFeatureClasses()
settings = {}
# settings['binWidth'] = 25
# settings['resampledPixelSpacing'] = [1.25, 1.25, 1.25]  
# settings['interpolator'] = 'sitkBSpline'
# settings['verbose'] = True
extractor.enableImageTypeByName('Wavelet')
extractor.enableImageTypeByName('LoG', customArgs={'sigma':[0.5, 1, 1.5, 2, 2.5]})

df = pd.DataFrame()
for folder in folderList:
	imageName = dataDir + folder + '/data.nii'
	maskName = dataDir + folder + '/mask.nii.gz'
	featureVector = extractor.execute(imageName, maskName)
	df_add = pd.DateFrame.from_dict(featureVector.values()).T
	df_add.columns = featureVector.keys()
	df = pd.concat([df, df_add])

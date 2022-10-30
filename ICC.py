####################################################################################################################
### title: "ICC"
### author: "Kzzhu"
### time: "2022.06.04"
####################################################################################################################
#! pip install pingouin -i https://pypi.tuna.tsinghua.edu.cn/simple
import pingouin as pg
import pandas as pd
import numpy as np
import os
import csv
from pathlib import Path
from pandas import Series, DataFrame

folderPath = r'E:/ICC'
data_1 = pd.read_csv('E:/intra.csv', sep = '\t')
data_2 = pd.read_csv('E:/inter.csv', sep = '\t')
#print(data_2)

data_1.insert(0, 'reader', np.ones(data_1.shape[0]))
data_2.insert(0, 'reader', np.ones(data_2.shape[0])*2)
#print(data_2)

data_1.insert(0, 'target', range(data_1.shape[0]))
data_2.insert(0, 'target', range(data_2.shape[0]))

data = pd.concat([data_1, data_2])

excel_file = r"E:/ICC.xlsx"
data = pd.read_excel(excel_file)
data.shape
# print(data)

rows = data.shape[0]
cols = data.shape[1]

target = data.iloc[:, [1]]
reader = data.iloc[:, [2]]

ICC = np.zeros((6, cols-3))
for i in range(3, cols):
    feature = data.iloc[:, [i]]
    data0 = np.hstack((target, reader, feature))
    data0 = DataFrame(data0, columns = ["target", "reader", "feature"])
    icc = pg.intraclass_corr(data = data0, targets = 'target', raters = 'reader', ratings = 'feature')
    ICC[:, [i-3]] = icc.iloc[:, [2]]
ICC1 = DataFrame(ICC)

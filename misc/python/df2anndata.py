import os
import pyreadr
import time
import anndata
import sys
import pandas as pd

os.chdir("/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer")

print('loading data ...')
print(time.strftime('>> %Y-%m-%d %H:%M:%S', time.localtime(time.time())))
data = pyreadr.read_r('data/12T_harmony_celltype.rds')
data = data[None].transpose()
print('data loaded')
print(time.strftime('>> %Y-%m-%d %H:%M:%S', time.localtime(time.time())))

# 如果是一个df
if isinstance(data, pd.DataFrame):
	print('data is a DataFrame')
    print(data.head())
    print(data.shape)
	# try to convert it to anndata
	print('try to convert it to anndata')
	print('>> %Y-%m-%d %H:%M:%S', time.localtime(time.time()))
	try:
		data = anndata.AnnData(data)
	except:
		print('data is a DataFrame, but can not be converted to anndata')
		sys.exit(1)
	print('data is converted to anndata')
	print('>> %Y-%m-%d %H:%M:%S', time.localtime(time.time()))
	# save it to h5ad
	print('save it to h5ad')
	print('>> %Y-%m-%d %H:%M:%S', time.localtime(time.time()))
	data.write('data/reference_gene.h5ad')
	print('data is saved')
	print('>> %Y-%m-%d %H:%M:%S', time.localtime(time.time()))
else:
	print('data is not a df')

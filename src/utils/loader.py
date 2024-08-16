
import src.utils.config as config
import sys

# Define a function to load data
def load_data(path_name, file_type):
	paths = {
		"tissue": config.tissue_path,
		"genus": config.genus_path,
		"species": config.species_path,
		"ref_h5": config.reference_h5_path,
		"ref_dict": config.reference_dict_path,
		"tissue_temp": config.tissue_temp_path,
		"reference": config.reference_path
	}

	if path_name not in paths:
		print("path_name error")
		sys.exit(1)

	path = paths[path_name]

	if file_type == "h5ad":
		import stereo as st
		return st.io.read_h5ad(
			file_path = path,
			bin_type="bins",
			bin_size = config.bin_size
		)
	elif file_type == "r":
		import pyreadr
		return pyreadr.read_r(path)[None]
	elif file_type == "gem":
		import stereo as st
		return st.io.read_gem(
			file_path = path,
			bin_type="bins",
			sep="\t",
			bin_size = config.bin_size
		)
	elif file_type == "txt.gz":
		import pandas as pd
		return pd.read_csv(
				path,
				compression = 'gzip',
				sep = '\t',
				header = 0
			)
	else:
		print("file_type error")
		sys.exit(1)

def save_temp(data, file):
	import stereo as st
	st.io.stereo_to_anndata(
	    data,
	    output = file
	)


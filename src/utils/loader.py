import src.utils.config as config
import sys
import stereo as st

# define a function to load data
def load_data(path_name):
    if path_name == "tissue":
        path = config.tissue_path
    elif path_name == "genus":
        path = config.genus_path
    elif path_name == "species":
        path = config.species_path
    else:
        print("path_name error")
        sys.exit(1)
    return st.io.read_gem(
        file_path = path,
        bin_type = "bins",
        sep = "\t",
        bin_size = config.bin_size
    )
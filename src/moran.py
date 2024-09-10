import numpy as np
import libpysal
import pysal.lib as ps
from esda.moran import Moran_Local
import geopandas as gpd


class Moran:
    def __init__(self, data, row, col):
        self.data = data
        self.w = libpysal.weights.lat2W(row, col)
        self.moran = Moran_Local(data, self.w)

    def I(self):
        return self.moran.Is

    def p_sim(self):
        return self.moran.p_sim

    def save(self, path):
        self.moran.to_file(path)

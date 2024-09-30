#!/bin/env python
"""
A wrapper script that shows how the solardemo library can be used
"""
from glob import glob
import sys
import os
from solardemo import solardemo

if __name__=="__main__":
    aia_131_dir = "/home/linn/july/solar/data/1690534254"
    aia_131_paths = glob(os.path.join(aia_131_dir, "*.fits"))


    sf = solardemo()

    sf.read_aia(key=131, file_paths = aia_131_paths)
    sf.aia[131] = solardemo.downscale_map(sf.aia[131], dim=[512, 512])

    flux_datapath = '/home/linn/july/solar/data/XRS/sci_gxrs-l2-irrad_g15_d20140107_v0-0-0.nc'
    flare_start = "2014-01-07T18:04:00"
    flare_end = "2014-01-07T18:58:00"

    aia_l15_maps = solardemo.l1_to_l15(sf.aia[131])
    solardemo.aia_to_png(aia_l15_maps, "aia_flare_512")


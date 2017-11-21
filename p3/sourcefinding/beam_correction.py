import numpy as np


def read_fitted_beam():
    beamfile = '/home/erichards/work/data/FITBEAM_FINAL_DOMEGA.txt'
    with open(beamfile, 'r') as f:
        lines = f.readlines()

    pbdir = {'angle' : [], 'power' : [], 'domega' : []}
    for line in lines:
        data = line.split()
        pbdir['angle'].append(float(data[0]))
        pbdir['power'].append(float(data[1]))
        pbdir['domega'].append(float(data[2]))

    return pbdir


def find_nearest_pbcorr(angle):
    pbdata = read_fitted_beam()
    idx = (np.abs(np.array(pbdata['angle']) - angle)).argmin()
    return pbdata['power'][idx]

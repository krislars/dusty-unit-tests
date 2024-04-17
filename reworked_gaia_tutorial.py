pip inimport numpy as np
from matplotlib import pyplot as plt
import speclite.filters
from dust_extinction.parameter_averages import F19 # more available!
import astropy.units as u
from astropy.io import fits

def get_type_modelname(type):
    # some text goes here
    stardict = { # from the STScI page referenced above
    'O3V'   : 'ckp00_45000[g45]',
    'O4V'   : 'ckp00_43000[g45]',
    'O5V'   : 'ckp00_41000[g45]',
    #'O5_5V': 'ckp00_40000[g40]', # yields NANs?
    'O6V'   : 'ckp00_39000[g40]',
    'O6_5V' : 'ckp00_38000[g40]',
    'O7V'   : 'ckp00_37000[g40]',
    'O7_5V' : 'ckp00_36000[g40]',
    'O8V'   : 'ckp00_35000[g40]',
    'O8_5V' : 'ckp00_34000[g40]',
    'O9V'   : 'ckp00_33000[g40]',
    'O9_5V' : 'ckp00_32000[g40]',
    'B0V'   : 'ckp00_30000[g40]',
    'B1V'   : 'ckp00_25000[g40]',
    'B3V'   : 'ckp00_19000[g40]',
    'B5V'   : 'ckp00_15000[g40]',
    'B8V'   : 'ckp00_12000[g40]',
    'A0V'   : 'ckp00_9500[g40]',
    'A1V'   : 'ckp00_9250[g40]',
    'A5V'   : 'ckp00_8250[g40]',
    'F0V'   : 'ckp00_7250[g40]',
    'F2V'   : 'ckp00_7000[g40]',
    'F5V'   : 'ckp00_6500[g40]',
    'F8V'   : 'ckp00_6250[g45]',
    'G0V'   : 'ckp00_6000[g45]',
    'G5V'   : 'ckp00_5750[g45]',
    'G8V'   : 'ckp00_5500[g45]',
    'K0V'   : 'ckp00_5250[g45]',
    'K2V'   : 'ckp00_4750[g45]',
    'K4V'   : 'ckp00_4500[g45]',
    'K5V'   : 'ckp00_4250[g45]',
    'K7V'   : 'ckp00_4000[g45]',
    'M0V'   : 'ckp00_3750[g45]',
    'M2V'   : 'ckp00_3500[g45]',
    'M5V'   : 'ckp00_3500[g50]',
    'B0III' : 'ckp00_29000[g35]',
    'B5III' : 'ckp00_15000[g35]',
    'G0III' : 'ckp00_5750[g30]',
    'G5III' : 'ckp00_5250[g25]',
    'K0III' : 'ckp00_4750[g20]',
    'K5III' : 'ckp00_4000[g15]',
    'M0III' : 'ckp00_3750[g15]',
    'O5I'   : 'ckp00_40000[g45]',
    'O6I'   : 'ckp00_39000[g40]',
    'O8I'   : 'ckp00_34000[g40]',
    'B0I'   : 'ckp00_26000[g30]',
    'B5I'   : 'ckp00_14000[g25]',
    'A0I'   : 'ckp00_9750[g20]',
    'A5I'   : 'ckp00_8500[g20]',
    'F0I'   : 'ckp00_7750[g20]',
    'F5I'   : 'ckp00_7000[g15]',
    'G0I'   : 'ckp00_5500[g15]',
    'G5I'   : 'ckp00_4750[g10]',
    'K0I'   : 'ckp00_4500[g10]',
    'K5I'   : 'ckp00_3750[g05]',
    'M0I'   : 'ckp00_3750[g00]',
    'M2I'   : 'ckp00_3500[g00]'}
    return stardict[type]

def test_get_type_modelname():
    # unit test
    assert get_type_modelname('A1V') == 'ckp00_9250[g40]'


def get_type_dat(type):
    # some text goes here
    modelname=get_type_modelname(type)
    filename=modelname[:-5]
    # Read file for the spectral type from STScI
    url = 'https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/ck04models/ckp00/'+filename+'.fits'
    star = fits.open(url)
    # Read in star model data
    wavedat = star[1].data['WAVELENGTH'] * u.AA
    fluxdat = star[1].data[modelname[-4:-1]] * u.erg / (u.cm ** 2 * u.AA * u.s)
    wav = wavedat[(wavedat > 2900*u.AA) & (wavedat < 12000*u.AA)] # clipped for dust_ext
    flux = fluxdat[(wavedat > 2900*u.AA) & (wavedat < 12000*u.AA)] # clipped for dust_ext
    return wav,flux

def test_get_type_dat():
    # unit test
    wav,flux = get_type_dat('A1V')
    assert len(flux) > 0


def get_type_temp(type):
    # some text goes here
    modelname = get_type_modelname(type)
    temp = get_type_modelname(type).split('_')[1].split('[')[0]
    return float(temp)

def test_get_type_temp():
    # unit test
    assert get_type_temp('A1V') == 9250.0


# Choose extinction curve!  Already in function and object form
ext = F19(Rv=3.1)


def get_norm_ext(type,band,ext,av):
    # some text goes here
    filters = speclite.filters.load_filters('gaiadr2-*')
    shortnames = [filters.names[0][-2:],filters.names[1][-1:],filters.names[2][-2:]]
    bandname = 'gaiadr2-'+band
    wav,flux =  get_type_dat(type)
    mags = filters.get_ab_magnitudes(flux, wav)
    flux_ext = flux*ext.extinguish(wav, Av=av)
    mags_ext = filters.get_ab_magnitudes(flux_ext, wav)
    return mags_ext[bandname][0]-mags[bandname][0]

def test_get_norm_ext():
    # unit test, although this only works if you already have the band
    photext = get_norm_ext('A1V','G',ext,1.0)
    assert 0.1 <= photext <= 2


Av = 2.0
ext = F19(Rv=3.1)
band = 'G'

types = ['B0V','A0V','F0V','G0V','K0V','M0V']

for sptype in types:
    plt.plot(get_type_temp(sptype),get_norm_ext(sptype,band,ext,Av)/Av,'ro')
plt.xlabel('temperature')
plt.ylabel('A_'+band+' / A_V')
plt.xscale('log')
plt.xlim(3000,50000)
plt.title('A_V = '+ str(Av))
plt.savefig("myfig.pdf", format="pdf")
plt.show()


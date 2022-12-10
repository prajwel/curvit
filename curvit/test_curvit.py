import curvit
import numpy as np

from astropy.io import fits

# TODO
def pseudo_data(source_coo = 2400,
               point_source_sigma = 3 / 2.355, 
               frames = 10000,
               framerate = 28.7185):
    
    Fx = np.random.normal(source_coo, point_source_sigma, frames)
    Fy = Fx[::-1]
    MJD_L2 = 2E8 + (np.arange(1, frames + 1) / framerate)
    ENP = np.ones(frames) * framerate
    
    col1 = fits.Column(name='Fx', format='D', array=Fx)
    col2 = fits.Column(name='Fy', format='D', array=Fy)
    col3 = fits.Column(name='MJD_L2', format='D', array=MJD_L2)
    col4 = fits.Column(name='EFFECTIVE_NUM_PHOTONS', format='D', array=ENP)
    
    hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4])
    
    return hdu

# TODO
def test_create_file(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "psuedo_data.fits"
    hdu = pseudo_data()
    hdu.writeto(p)
    assert curvit.curvit.makecurves(p) == 1
    assert len(curvit.curvit.read_columns(p)) == 4


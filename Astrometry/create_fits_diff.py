from astropy.io import fits
import os

fits_folder = r"C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B\HI-1"
fits_out_folder = r"C:\PhD\Comet_data\Comet_PanSTARRS_C2011L4\Gallery\Stereo_B\HI-1-diff"

fits_list = os.listdir(fits_folder)
fits_list = [s for s in fits_list if 'fits' in s]

for fits_no in xrange(1,len(fits_list)):
    
    fits_string = fits_list[fits_no].split(".")[0]
    
    fits_current = os.path.join(fits_folder, fits_list[fits_no])
    fits_previous = os.path.join(fits_folder, fits_list[fits_no-1])
    
    data_cur, header_cur = fits.getdata(fits_current, header=True)
    data_pre, header_pre = fits.getdata(fits_previous, header=True)

    data_out = data_cur - data_pre
       
    fits_outfile = os.path.join(fits_out_folder, fits_string + "_diff.fits")    
    fits.writeto(fits_outfile, data_out, header_cur, clobber=True)
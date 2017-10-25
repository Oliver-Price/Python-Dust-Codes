from astropy.io import fits
import os

fits_small_folder = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_B\HI-1-MGN\mixedup_header'
fits_big_folder = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Gallery\Stereo_B\HI-1\originals\headerless\extras'

fits_small_list = os.listdir(fits_small_folder)
fits_small_list = [s for s in fits_small_list if 'fits' in s]

fits_big_list = os.listdir(fits_big_folder)
fits_big_list = [s for s in fits_big_list if 'fts' in s]

fits_left_out = []
savemode = False
if savemode == True:
    fits_diff_folder = os.path.join(fits_big_folder,"extras")
    os.makedirs(fits_diff_folder)

for fits_no in range(0,len(fits_big_list)):
    
    fits_string = fits_big_list[fits_no].split(".")[0]
    if (len([s for s in fits_small_list if fits_big_list[fits_no].split('.')[0] in s])==0):
        fits_left_out.append(fits_big_list[fits_no])
        if savemode == True:
            data, header = fits.getdata(os.path.join(fits_big_folder,fits_big_list[fits_no]), header=True)
            fits_outfile = os.path.join(fits_diff_folder, fits_big_list[fits_no])    
            fits.writeto(fits_outfile, data, header, clobber=True)
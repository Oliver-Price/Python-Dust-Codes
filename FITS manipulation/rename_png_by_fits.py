import os

fitsdir = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Druckmuller\fits'
pngdir = r'C:\PhD\Comet_data\Comet_McNaught_C2006P1\Druckmuller\cropped'

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)

for fits_id in range(1,fits_total):
    
    try:  
        basename = fits_list[fits_id].split('.')[0]
        notime = basename[24:].lower()
        pngin = os.path.join(pngdir,notime+".png")
        pngout = os.path.join(pngdir,basename + ".png")
        print(pngin,pngout)
        os.rename(pngin, pngout)
    except:
        pass
import os

fitsdir = r'C:\PhD\Comet_data\Comet for Deep Learning\Images\Colour\fits'

dir_list = sorted(os.listdir(fitsdir))
fits_list = [s for s in dir_list if ".f" in s]
fits_total = len(fits_list)

for fits_id in range(0,fits_total):
    
    fileold = os.path.join(fitsdir,fits_list[fits_id])

    stringbits = fits_list[fits_id].split('_')
    if len(stringbits[1]) == 4:
        pass
    else:
        stringbits[1] = stringbits[1][0:4] + '_' + stringbits[1][4:6] + '_' + stringbits[1][6:8]
        
    filename = os.path.join(fitsdir,'_'.join(stringbits))
    
    os.rename(fileold,filename)
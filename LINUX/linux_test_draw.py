comimg = Image.open(imgsav)
d = ImageDraw.Draw(comimg)

rafiles = r"C:\PhD\Python\test_directory\ra.txt"
decfiles = r"C:\PhD\Python\test_directory\dec.txt"
boolfiles = r"C:\PhD\Python\test_directory\inimg.txt"

ra_lin = np.loadtxt(rafiles)
dec_lin = np.loadtxt(decfiles)
point_bools = np.loadtxt(boolfiles)

ra_lin_pix = ra2xpix(ra_lin,border,pixwidth,rafmin,scale)
dec_lin_pix = dec2ypix(dec_lin,border,pixheight,decmin,scale)

for x in xrange(0, 100):
    for y in xrange(0,100):
        b = d.line([(ra_lin_pix[x,y]+1,dec_lin_pix[x,y]+1), \
                (ra_lin_pix[x,y]-1,dec_lin_pix[x,y]-1)],\
                fill = (255,0,int(255*point_bools[x,y]),255))
        b = d.line([(ra_lin_pix[x,y]-1,dec_lin_pix[x,y]+1), \
                (ra_lin_pix[x,y]+1,dec_lin_pix[x,y]-1)],\
                fill = (255,0,int(255*point_bools[x,y]),255))

comimg.show()
non_zeros_0 = simres[:,:,14].nonzero()[0]
non_zeros_1 = simres[:,:,14].nonzero()[1]
no_points = non_zeros_0.size
nu_radecs = np.zeros((no_points,2))

for x in xrange(0,no_points):
    nu_radecs[x,0] = simres[non_zeros_0[x],non_zeros_1[x],10]
    nu_radecs[x,1] = simres[non_zeros_0[x],non_zeros_1[x],11]

pixlocs = w.wcs_world2pix(nu_radecs,0)

simres_pix= np.zeros((tno,bno,2))
for x in xrange(0,no_points):
    simres_pix[non_zeros_0[x],non_zeros_1[x],0] = pixlocs[x][0]
    simres_pix[non_zeros_0[x],non_zeros_1[x],1] = pixlocs[x][1]
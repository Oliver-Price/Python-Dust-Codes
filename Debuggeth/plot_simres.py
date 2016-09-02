
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 35})
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
locs = np.where(simres[:,:,14] == 1)
ts = ax.scatter(simres[locs[0],locs[1],0],simres[locs[0],locs[1],1], s=10)
plt.xlabel('Time of Ejection (days)')
plt.ylabel('Beta')


#theta = lat, phi = long
def rcolatlon2xyz(r,theta,phi):
    x = r * sind(theta) * cosd(phi)
    y = r * sind(theta) * sind(phi)
    z = r * cosd(theta)
    return (x,y,z)

def xyz2rcolatlon(x,y,z):
    r = np.sqrt(x*x + y*y + z*z)
    phi = np.arctan2(y,x)
    theta = np.arccos(z/r)
    return (r,np.degrees(theta)%360,np.degrees(phi)%360)

#%%

def rlatlon2xyz(r,theta,phi):
    theta = 90 - theta
    x = r * sind(theta) * cosd(phi)
    y = r * sind(theta) * sind(phi)
    z = r * cosd(theta)
    return (x,y,z)

def xyz2rlatlon(x,y,z):
    r = np.sqrt(x*x + y*y + z*z)
    phi = np.arctan2(y,x)
    theta = np.arccos(z/r)
    return (r,(90-np.degrees(theta)%360),np.degrees(phi)%360)
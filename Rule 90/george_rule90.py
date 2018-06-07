from PIL import Image, ImageDraw
import numpy as np
import time

t0 = time.clock()

#Calculate
steps = 100
initial = [50,35,32,79,7,67,40]#indices of cells in first row that =1
arr = np.zeros([steps, steps])
for index in initial:
    arr[0][index] = 1
for n in range(steps - 1):
    for i in range(steps):
        arr[n+1][i] = int(arr[n][i-1])^int(arr[n][(i+1)%steps])
        
t1 = time.clock()

#Draw
size = 10
im = Image.new('RGBA', (steps * size, steps * size))
draw = ImageDraw.Draw(im)
for y, row in enumerate(arr):
    for x, cell in enumerate(row):
        f = int((1 - cell) * 255)
        draw.rectangle((x * size, y * size, (x + 1) * size, (y + 1) * size), fill = (f, f, f))
#im.show()
#im.save('rule90.png')

t2 = time.time()
print(t1-t0)
#print(t2-t0)
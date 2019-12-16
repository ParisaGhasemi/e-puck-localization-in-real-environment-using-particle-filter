from Gui import *

sync = "sync.txt"
data = "datafile.txt"
numparticles = 3000
df = open(data, 'r')
sf = open(sync, 'r')

count = 0
particles = list()
GS = 3

env = open("env.txt", 'r')
dims = env.readline().strip().split()
HEIGHT = int(dims[0])
WIDTH = int(dims[1])
SCALE = float(dims[2])
print( "height: ", HEIGHT, "\n", "width: ", WIDTH, "\nSCALE: ", SCALE, "\n")
env.close()

mapp = []
for i in range(HEIGHT):
    mapp.append([])

mp = open("map.txt", 'r')
i = HEIGHT - 1
while i >= 0:
    l = mp.readline().strip().split(' ')
    mapp[i] = [int(X) for X in l]
    i -= 1
mp.close()

g = Gui()
g.title('particles')

canvas = g.ca(width=WIDTH * 3, height=HEIGHT * 3)
canvas.config(bg='white')

points = []
for i in range(numparticles):
    points.append(None)
xoffset = (WIDTH / 2) * GS
yoffset = (HEIGHT / 2) * GS

sum = [0.0,0.0,0.0]
def start():
    global count, points
    cline = sf.readline()
    while not cline:
        cline = sf.readline()
    newcount = int(cline.strip())
    if newcount > count:
        particles.append(list())
        line = df.readline()
        while (line[0:3] != 'END'):
            pos = [float(x) for x in line.strip().split()]
            sum[0]+=pos[0]
            sum[1]+=pos[1]
            sum[2]+=pos[2]
            t = tuple(pos)  #
            particles[count].append(t)
            line = df.readline()

        count += 1
    elif newcount == -1:
        return

    for i in range(numparticles):
        if points[i]:
            points[i].delete()
        points[i] = canvas.circle(
            (particles[count - 1][i][0] * GS - xoffset, particles[count - 1][i][1] * GS - yoffset), 2, fill='yellow')
    sum[0]/=numparticles
    sum[1]/=numparticles
    sum[2]/=numparticles

    g.after(20, start)

button1 = g.bu('press', command=start)
for i in range(HEIGHT):
    for j in range(WIDTH):
        if mapp[i][j] == 1:
            x = j * GS
            y = i * GS
            sq = canvas.rectangle([[x - xoffset, y - yoffset], [x + GS - xoffset, y + GS - yoffset]], outline='green')

g.after(100, start)
g.mainloop()

sf.close()
df.close()

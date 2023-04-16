import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

nx, ny = (201, 201)
n = 2
pi = 3.1415926
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
print(nx)
print(x[2] - x[1])
xv, yv = np.meshgrid(x, y)
print(xv)
u = np.multiply(np.cos(n*pi*xv), np.sin(n*pi*yv))
v = -np.multiply(np.sin(n*pi*xv), np.cos(n*pi*yv))
p = -0.25*(np.cos(2*n*pi*xv) + np.cos(2*pi*n*yv))
px = 0.5*n*pi*np.sin(2*pi*n*xv)
py = 0.5*n*pi*np.sin(2*n*pi*yv)

# case 1
k1 = 0.1
m = 0.5
epx1 = k1*np.cos(2*m*pi*xv)
epy1 = k1*np.cos(2*m*pi*yv)

p1 = np.zeros((nx, ny))
erp1 = np.zeros((nx, ny))
with open('case_1_pressure.txt', 'r') as f:
    f.readline()
    f.readline()
    f.readline()
    for i in range(nx):
        for j in range(ny):
            line = f.readline()
            vars = line.split(" ")
            p1[j, i] = vars[4]
            erp1[j, i] = p1[j, i] - p[j, i]
f.close()

# case 2
k2 = k1*math.sqrt(2);
k = 2*m;
epx2 = k2*np.multiply(np.cos(k*pi*xv), np.sin(k*pi*yv))
epy2 = -k2*np.multiply(np.sin(k*pi*xv), np.cos(k*pi*yv))
p2 = np.zeros((nx, ny))
erp2 = np.zeros((nx, ny))

with open('case_2_pressure.txt', 'r') as f:
    f.readline()
    f.readline()
    f.readline()
    for i in range(nx):
        for j in range(ny):
            line = f.readline()
            vars = line.split(" ")
            p2[j, i] = vars[4]
            erp2[j, i] = p2[j, i] - p[j, i]
f.close()

# case 3
k3 = k1*math.sqrt(2)/2;
epx3 = k3*np.ones((nx, ny))
epy3 = k3*np.ones((nx, ny))

p3 = np.zeros((nx, ny))
erp3 = np.zeros((nx, ny))

with open('case_3_pressure.txt', 'r') as f:
    f.readline()
    f.readline()
    f.readline()
    for i in range(nx):
        for j in range(ny):
            line = f.readline()
            vars = line.split(" ")
            p3[j, i] = vars[4]
            erp3[j, i] = p3[j, i] - p[j, i]
f.close()


fig, axs = plt.subplots(2, 2)
cmap = cm.get_cmap('viridis')
levels = np.linspace(-1, 1, 21)
cs = axs[0, 0].contourf(xv, yv, p, levels, cmap=cmap)
axs[0, 0].set_title('Real pressure')
axs[0, 0].set_xlim(0, 1)
axs[0, 0].set_ylim(0, 1)
fig.colorbar(cs, ax=axs[0, 0])
cs = axs[0, 1].contourf(xv, yv, p1, levels, cmap=cmap)
axs[0, 1].set_title('Case1 pressure')
axs[0, 1].set_xlim(0, 1)
axs[0, 1].set_ylim(0, 1)
fig.colorbar(cs, ax=axs[0, 1])
cs = axs[1, 0].contourf(xv, yv, p2, levels, cmap=cmap)
axs[1, 0].set_title('Case2 pressure')
axs[1, 0].set_xlim(0, 1)
axs[1, 0].set_ylim(0, 1)
fig.colorbar(cs, ax=axs[1, 0])
cs = axs[1, 1].contourf(xv, yv, p3, levels, cmap=cmap)
axs[1, 1].set_title('Case3 pressure')
axs[1, 1].set_xlim(0, 1)
axs[1, 1].set_ylim(0, 1)
fig.colorbar(cs, ax=axs[1, 1])
fig.set_size_inches(16, 12, forward=True)
fig.savefig("pressure.png")


fig, axs = plt.subplots(2, 2)
cmap = cm.get_cmap('viridis')
levels = np.linspace(-1, 1, 21)
cs = axs[0, 0].contourf(xv, yv, p, levels, cmap=cmap)
axs[0, 0].set_title('Real pressure')
axs[0, 0].set_xlim(0, 1)
axs[0, 0].set_ylim(0, 1)
fig.colorbar(cs, ax=axs[0, 0])
levels = np.linspace(-0.2, 0.2, 21)
cs = axs[0, 1].contourf(xv, yv, epx1, levels, cmap=cmap)
axs[0, 1].set_title('Case1 epx')
axs[0, 1].set_xlim(0, 1)
axs[0, 1].set_ylim(0, 1)
fig.colorbar(cs, ax=axs[0, 1])
cs = axs[1, 0].contourf(xv, yv, epx2, levels, cmap=cmap)
axs[1, 0].set_title('Case2 epx')
axs[1, 0].set_xlim(0, 1)
axs[1, 0].set_ylim(0, 1)
fig.colorbar(cs, ax=axs[1, 0])
cs = axs[1, 1].contourf(xv, yv, epx3, levels, cmap=cmap)
axs[1, 1].set_title('Case3 epx')
axs[1, 1].set_xlim(0, 1)
axs[1, 1].set_ylim(0, 1)
fig.colorbar(cs, ax=axs[1, 1])
fig.set_size_inches(16, 12, forward=True)
fig.savefig("Dpdx.png")


fig, axs = plt.subplots(2, 2)
cmap = cm.get_cmap('viridis')
levels = np.linspace(-1, 1, 21)
cs = axs[0, 0].contourf(xv, yv, p, levels, cmap=cmap)
axs[0, 0].set_title('Real pressure')
axs[0, 0].set_xlim(0, 1)
axs[0, 0].set_ylim(0, 1)
fig.colorbar(cs, ax=axs[0, 0])
levels = np.linspace(-0.08, 0.08, 17)
rms = math.sqrt(np.square(erp1).mean())
avg = np.abs(erp1).mean()
cs = axs[0, 1].contourf(xv, yv, erp1, levels, cmap=cmap)
axs[0, 1].set_title(f'Case 1, pDiff. RMS={round(rms, 4)}, Avg={round(avg, 4)}')
norm = mpl.colors.Normalize(vmin=-1, vmax=1)
pcm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
fig.colorbar(cs, ax=axs[0, 1])

cs = axs[1, 0].contourf(xv, yv, erp2, levels, cmap=cmap)
rms = math.sqrt(np.square(erp2).mean())
avg = np.abs(erp2).mean()
axs[1, 0].set_title(f'Case 2, pDiff. RMS={round(rms, 4)}, Avg={round(avg, 4)}')
norm = mpl.colors.Normalize(vmin=-1, vmax=1)
pcm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
fig.colorbar(cs, ax=axs[1, 0])

cs = axs[1, 1].contourf(xv, yv, erp3, levels, cmap=cmap)
rms = math.sqrt(np.square(erp3).mean())
avg = np.abs(erp3).mean()
axs[1, 1].set_title(f'Case 3, pDiff. RMS={round(rms, 4)}, Avg={round(avg, 4)}')
norm = mpl.colors.Normalize(vmin=-1, vmax=1)
pcm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
fig.colorbar(cs, ax=axs[1, 1])
fig.set_size_inches(16, 12, forward=True)
fig.savefig("Perror.png")


#plt.show()



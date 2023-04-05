import numpy as np
import math
import matplotlib.pyplot as plt

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
epx1 = 0.5*k1*n*pi*np.cos(2*m*pi*xv)
epy1 = 0.5*k1*n*pi*np.cos(2*m*pi*yv)

with open('case_1_m0.5.txt', 'w') as f:
    f.write("TITLE = \"Case 1: irrotational error\"\n")
    f.write("VARIABLES=\"X\" \"Y\" \"DuDt\" \"DvDt\" \"PReal\"\n")
    f.write(f"ZONE I={nx}, J={ny}\n")
    for i in range(nx):
        for j in range(ny):
            f.write(f"{xv[j, i]} {yv[j, i]} {-px[j, i] - epx1[j, i]} {-py[j, i] - epy1[j, i]} {p[j, i]}\n")
f.close()

# case 2
k2 = 0.1;
m = 4;
epx2 = k2*0.5*n*pi*np.multiply(np.cos(m*pi*xv), np.sin(m*pi*yv))
epy2 = k2*0.5*n*pi*np.multiply(np.sin(m*pi*xv), np.cos(m*pi*yv))

with open('case_2.txt', 'w') as f:
    f.write("TITLE = \"Case 2: an solenoidal error in grad P\"\n")
    f.write("VARIABLES=\"X\" \"Y\" \"DuDt\" \"DvDt\" \"PReal\"\n")
    f.write(f"ZONE I={nx}, J={ny}\n")
    for i in range(nx):
        for j in range(ny):
            f.write(f"{xv[j, i]} {yv[j, i]} {-px[j, i] - epx2[j, i]} {-py[j, i] - epy2[j, i]} {p[j, i]}\n")
f.close()

# case 3
k3 = 0.1;
epx3 = k3*0.5*n*pi*np.ones((nx, ny))
epy3 = k3*0.5*n*pi*np.ones((nx, ny))

with open('case_3.txt', 'w') as f:
    f.write("TITLE = \"Case 3: an harmonic error in grad P\"\n")
    f.write("VARIABLES=\"X\" \"Y\" \"DuDt\" \"DvDt\" \"PReal\"\n")
    f.write(f"ZONE I={nx}, J={ny}\n")
    for i in range(nx):
        for j in range(ny):
            f.write(f"{xv[j, i]} {yv[j, i]} {-px[j, i] - epx3[j, i]} {-py[j, i] - epy3[j, i]} {p[j, i]}\n")
f.close()

### Parallel line Omni2D integration

arguments.txt file
Imax=32
Jmax=32
rho=1.0
lineSpacing=1.0
eps=1e-10
numOfAngles=1000
numOfIterations=3
dx=0.006
dy=0.006
dataFile=data.txt
numOfHeaderLines=0

contains the arguments
lineSpacing=1.0 means parallel line spacing is 1.0dx i.e. 1 grid spacing
numOfAngles: how many angles are there, follow which parallel line group will be generated.
dataFile=data.txt, is the acceleration, each column is
x y dudt dvdt
J dimension comes first
numOfHeaderLines: header lines in data.txt

outFile=out.txt, output file, each column is
x y dudt dvdt p pCount

pCount is how many times this point is crossed during integration
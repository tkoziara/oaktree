simu = SIMULATION ('out/well', 1, 0.001, 0.003)

width = 10.0

cube = CUBE ((0, 0, 0), width, width, width, (1, 1, 1, 1, 1, 1))

pipe = CYLINDER ((0.5*width, 0.5*width, -width), 3*width, 0.15, (2, 2, 2))

ROTATE (pipe, (0.5*width, 0.5*width, 0.5*width), (1, 1, 1), 15)

well = DIFFERENCE (cube, pipe)

SOLID (simu, well, 'well')

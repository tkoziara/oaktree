simu = SIMULATION ('out/cubes', 1.0, 0.001, 1, 0.01, (-1, -1, -1, 2, 2, 2))

cube = CUBE ((0, 0, 0), 1, 1, 1, 1, (1, 2, 3, 4, 5, 6))

ROTATE (cube, (0, 0, 0), (1, 1, 1), 10)

SOLID (simu, cube, "cube")

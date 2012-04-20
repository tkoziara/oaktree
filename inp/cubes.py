simu = SIMULATION ('out/cubes', 1.0, 0.001, 1.0/4.0, 0.001, (0, 0, 0, 1, 1, 1))

cube = CUBE ((0, 0, 0), 1, 1, 1, 1, (1, 2, 3, 4, 5, 6))

SOLID (simu, cube, "cube")

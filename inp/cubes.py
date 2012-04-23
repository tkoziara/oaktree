simu = SIMULATION ('out/cubes', 1.0, 0.001, 1, 0.01, (-1, -1, -1, 2, 2, 2))

cube1 = CUBE ((0, 0, 0), 1, 1, 1, 1, (1, 2, 3, 4, 5, 6))

cube2 = CUBE ((0.5, 0.5, 0.5), 1, 1, 1, 1, (1, 2, 3, 4, 5, 6))

ROTATE (cube2, (0, 0, 0), (1, 1, 1), 10)

cube = DIFFERENCE (cube1, cube2)

cube3 = CUBE ((0.75, 0.5, 0.25), 0.1, 1, 1, 1, (1, 2, 3, 4, 5, 6))

cube = DIFFERENCE (cube, cube3)

sphere = SPHERE ((1, 0, 1), 0.4, 1, 1)

cube = DIFFERENCE (cube, sphere)

SOLID (simu, cube, "cube")

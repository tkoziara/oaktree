simu = SIMULATION ('out/cube', 1.0, 0.001, 0.001)

a = CUBE ((0, 0, 0), 1, 1, 1, (1, 2, 3, 4, 5, 6))
b = CUBE ((0.5, 0.5, 0.5), 1, 1, 1, (1, 2, 3, 4, 5, 6))
ROTATE (b, (0, 0, 0), (1, 1, 1), 10)
c = DIFFERENCE (a, b)
a = CUBE ((0.75, 0.5, 0.25), 0.1, 1, 1, (1, 2, 3, 4, 5, 6))
c = DIFFERENCE (c, a)
a = SPHERE ((1, 0, 1), 0.4, 1)
c = DIFFERENCE (c, a)

SOLID (simu, c)

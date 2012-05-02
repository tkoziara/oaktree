simu = SIMULATION ('out/cylinders', 1.0, 0.001, 0.001, (-1, -1, -1, 1, 1, 2))

a = CYLINDER ((0, 0, 0), 1, 1, 1, (1, 1, 1))
b = CYLINDER ((0, 0, 0), 1, 0.8, 1, (1, 1, 1))
c = DIFFERENCE (a, b)
a = CYLINDER ((-1, 0, 0.5), 2, 0.25, 1, (1, 1, 1))
ROTATE (a, (-1, 0, 0.5), (0, 1, 0), 90)
c = DIFFERENCE (c, a)
a = CYLINDER ((0, -1, 0.5), 2, 0.25, 1, (1, 1, 1))
ROTATE (a, (0, -1, 0.5), (1, 0, 0), -90)
c = DIFFERENCE (c, a)

SOLID (simu, c, "c")

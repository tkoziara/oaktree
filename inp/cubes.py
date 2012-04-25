simu = SIMULATION ('out/cubes', 1.0, 0.001, 1, 0.001, (-1, -1, -1, 2, 2, 2))

x = 0
y = 0
z = 0

a = CUBE ((x, y, z), 1, 1, 1, 1, (1, 2, 3, 4, 5, 6))
b = CUBE ((x+0.5, y+0.5, z+0.5), 1, 1, 1, 1, (1, 2, 3, 4, 5, 6))
ROTATE (b, (x, y, z), (1, 1, 1), 10)
c = DIFFERENCE (a, b)
a = CUBE ((x+0.75, y+0.5, z+0.25), 0.1, 1, 1, 1, (1, 2, 3, 4, 5, 6))
c = DIFFERENCE (c, a)
a = SPHERE ((x+1, y+0, z+1), 0.4, 1, 1)
c = DIFFERENCE (c, a)

a = CYLINDER ((0, 0, 0), 1, 1, 1, (1, 1, 1))
b = CYLINDER ((0, 0, 0), 1, 0.8, 1, (1, 1, 1))
c = DIFFERENCE (a, b)
a = CYLINDER ((0.4, 0, 0), 1, 0.8, 1, (1, 1, 1))
c = DIFFERENCE (c, a)

SOLID (simu, c, "c")

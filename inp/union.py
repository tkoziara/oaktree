simu = SIMULATION ('out/union', 1.0, 0.001, 0.001)

a = CUBE ((1, -1, 0), 1, 2, 3, (1, 1, 1, 1, 1, 1))
b = CYLINDER ((1, 0, 0), 1, 1, (1, 1, 1))
ROTATE (b, (1, 0, 0), (0, 1, 0), 90)
c = UNION (a, b)
a = COPY (b)
MOVE (a, (0, 0, 3))
c = UNION (a, c)
a = COPY (b)
MOVE (a, (0, 1, 1.5))
c = UNION (a, c)
a = COPY (b)
MOVE (a, (0, -1, 1.5))
c = UNION (a, c)

ROTATE (c, (0, 0, 0), (1, 1, 1), 45)

SOLID (simu, c, "c")

# Rev 23 added more filtering to octree_insert_solid to avoid refinement
# in unions on internal boundaries of hidden primitives; this is a test for that

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

SOLID (simu, c, "c")

simu = SIMULATION ('out/implant', 1.0, 0.001, 0.0001)

a = SPHERE ((0, 0, 0), 0.025, 1)
b = SPHERE ((0, 0, 0), 0.02, 1)
c = DIFFERENCE (a, b)
a = CUBE ((-0.025, -0.025, -0.025), 0.05, 0.025, 0.025, (1, 1, 1, 1, 1, 1))
c = INTERSECTION (a, c)

a = SPHERE ((0, 0, 0), 0.018, 2)
b = CYLINDER ((0, 0, 0), 0.05, 0.01, (2, 2, 2))
d = UNION (a, b)

DOMAIN (simu, c)
DOMAIN (simu, d)

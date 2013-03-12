simu = SIMULATION ('out/mls', 1.0, 0.001, 0.005)

pv = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0),
      (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]

op = []

for p in pv:
  nx = p[0] - 0.5
  ny = p[1] - 0.5
  nz = p[2] - 0.5
  op.append ((p[0], p[1], p[2], nx, ny, nz))

a = MLS (op, 0.5, 1)

DOMAIN (simu, a)

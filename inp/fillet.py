simu = SIMULATION ('out/fillet', 1.0, 0.001, 0.001)

a = CUBE ((0, 0, 0), 1, 1, 1, (1, 2, 3, 4, 5, 6))
b = CUBE ((0.5, 0.5, 0.5), 1, 1, 1, (1, 2, 3, 4, 5, 6))
c = DIFFERENCE (a, b)

FILLET (c, (0, 0, 0), 0.1, 0.1, 1)
FILLET (c, (0.5, 0.5, 0.5), 0.1, 0.1, 1)

#FILLET (c, (0, 0, 1.0), 0.1, 0.1, 1) #FIXME: handling corners (temporary common fillet creation)
#FILLET (c, (0.5, 0.5, 1.0), 0.1, 0.1, 1) #FIXME: mixed case and putting the fillet in the right place in the tree

SOLID (simu, c)

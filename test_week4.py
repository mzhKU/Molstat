from forces import lennardjones
import numpy as np

U = [[0.0,1.4],[0.0,0.0],[0.0,0.0]]
box = [10.0,10.0,10.0]
# the code should print:
#
# (-0.44437003107014544, array([[ 1.6719969, -1.6719969]
#
# as the first line if everything is done correctly
#
F = lennardjones(U,box)[1]

F0 = F

U1 = [[0.0,2.4],[0.0,0.0],[0.0,0.0]]
box1 = [10.0,10.0,10.0]

F = lennardjones(U1, box1)[1]

print F0
print F


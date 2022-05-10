from geom import *
import numpy as np
from math import *

def initialize(coords):
    points = [OPoint(x, y) for x, y in coords]
    for i in range(len(points)):
        points[i].index = i
        points[i].neighbors = [[points[j], distance(points[i], points[j])] for j in range(len(points)) if j != i]
        points[i].update_distances(999, 999) # force recalc
    return points

def uniform_start(n):
    # uniformly sample n points from the unit square
    return np.random.random((n,2))

def grid_start(n):
    # sample n points from a grid
    maxcoord = int(ceil(sqrt(n)))
    return [((i+.5)/maxcoord, (j+.5)/maxcoord) for i in range(maxcoord) for j in range(maxcoord)][:n]

def hexagonal_start(n):

    def nrealized(nrow, ncol):
        return (nrow//2)*(ncol + ncol-1) + (nrow%2)*ncol
    def radius(nrow, ncol):
        return min(1/(2*ncol), 1/(sqrt(3)*(nrow-1)+2))

    # of rows is sqrt(3)/2 * n
    # of columns is just n
    numrow, numcol = None, None
    bestrad = 0
    for nrow in range(1, int(ceil(sqrt(n*2)))):
        for ncol in range(2, int(ceil(sqrt(n*2)))):
            if nrealized(nrow, ncol) >= n:
                r = radius(nrow, ncol)
                if r > bestrad:
                    bestrad = r
                    numrow, numcol = nrow, ncol
                break
    coords = []
    rowindex, colindex = 0, 0
    colrad = 1 / (2 * numcol)
    while len(coords) < n:
        if colindex == 0 and n-len(coords) < numcol:
            # we are starting the last line. Let's put it in the middle.
            colindex += (numcol - (n-len(coords))) // 2
        if rowindex % 2 == 0:
            coords.append(((2*colindex+1)*colrad, rowindex*sqrt(3)*bestrad + bestrad))
        else:
            coords.append(((2*colindex+2)*colrad, rowindex*sqrt(3)*bestrad + bestrad))
        colindex += 1
        if (rowindex%2 == 0 and colindex == numcol) or (rowindex%2 == 1 and colindex == numcol-1):
            # start a new row
            rowindex += 1
            colindex = 0
    return coords

import display
from collections import namedtuple
import numpy as np

def distance(p1, p2):
    return ((p1.x - p2.x)**2 + (p1.y - p2.y)**2)**0.5

def walldist(p):
    return min(p.x, p.y, 1 - p.x, 1 - p.y)*2

def add_jitter(p, magnitude):
    p.x += np.random.random() * magnitude - magnitude/2
    p.y += np.random.random() * magnitude - magnitude/2

def fix_bounds(p):
    # these bounds correspond to >250000 points,
    # and it makes our life a bit easier to not
    # have to worry about them being on the edge
    p.x = max(p.x, .001)
    p.x = min(p.x, .999)
    p.y = max(p.y, .001)
    p.y = min(p.y, .999)

def min_distance(points, xp, yp, split='x', basecase=5):
    # Find the minimum distance between any two points in the base case:
    if basecase is None or len(xp) < basecase:
        d = float('inf')
        for i in range(len(xp)):
            for j in range(i + 1, len(xp)):
                d = min(d, distance(points[xp[i]], points[xp[j]]))
        return d
    # Find the minimum distance between any two points in the recursive case:
    elif split == 'x':
        # Partition points by median x-coordinate:
        midx = len(xp) // 2
        xp1, xp2 = xp[:midx], xp[midx:]
        # To efficiently partition y-coordinates, we need to set-ify. Fortunately, we can do this in linear time.
        xp1s, xp2s = set(xp1), set(xp2)
        yp1, yp2 = [y for y in yp if y in xp1s], [y for y in yp if y in xp2s]
        # Recursive calls to find min distance
        d1, d2 = min_distance(points, xp1, yp1, 'y', basecase), min_distance(points, xp2, yp2, 'y', basecase)
        d = min(d1, d2)
        # We now need to scan the strip.
        strip = [points[xp[i]] for i in range(len(xp)) if abs(points[xp[i]].x - points[xp[midx]].x) < d]
        for i in range(len(strip)):
            for j in range(i + 1, len(strip)):
                d = min(d, distance(strip[i], strip[j]))
        return d
    elif split == 'y':
        # Partition points by median y-coordinate:
        midy = len(yp) // 2
        yp1, yp2 = yp[:midy], yp[midy:]
        # To efficiently partition x-coordinates, we need to set-ify. Fortunately, we can do this in linear time.
        yp1s, yp2s = set(yp1), set(yp2)
        xp1, xp2 = [x for x in xp if x in yp1s], [x for x in xp if x in yp2s]
        # Recursive calls to find min distance
        d1, d2 = min_distance(points, xp1, yp1, 'x', basecase), min_distance(points, xp2, yp2, 'x', basecase)
        d = min(d1, d2)
        # We now need to scan the strip.
        strip = [points[yp[i]] for i in range(len(yp)) if abs(points[yp[i]].y - points[yp[midy]].y) < d]
        for i in range(len(strip)):
            for j in range(i + 1, len(strip)):
                d = min(d, distance(strip[i], strip[j]))
        return d

def max_radius(points):
    # find sorted indices of points
    sorted_points_x = sorted([(i, p) for i, p in enumerate(points)], key=lambda p: p[1].x)
    sorted_points_y = sorted([(i, p) for i, p in enumerate(points)], key=lambda p: p[1].y)
    sorted_indices_x = [x[0] for x in sorted_points_x]
    sorted_indices_y = [x[0] for x in sorted_points_y]
    minx = sorted_points_x[0][1].x
    maxx = sorted_points_x[-1][1].x
    miny = sorted_points_y[0][1].y
    maxy = sorted_points_y[-1][1].y
    min_wall_dist = min(minx, miny, 1-maxx, 1-maxy) # ensures far from edges
    min_p2p_dist = min_distance(points, sorted_indices_x, sorted_indices_y)
    return min(min_wall_dist, min_p2p_dist/2)

LPoint = namedtuple('LPoint', ['x', 'y', 'is_wall'], defaults=(None, None, False))

class OPoint: # for optimization, as opposed to drawing.
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.index = None
        self.is_wall = False
        self.neighbors = []
        self.self_dist_since_update = 0
        self.neighbor_dist_since_update = 0
        self.num_recheck = None
    def update_distances(self, selfupdate, neighborupdate, num_interest=1):
        if len(self.neighbors) < 10:
            self.num_recheck = len(self.neighbors)
            return
        # num_interest represents how many nearest neighbors we may
        # want to query in the future. If it's 1, we only care about the closest.
        # if it's 3, we might care about the next 2 as well.
        self.self_dist_since_update += selfupdate
        self.neighbor_dist_since_update += neighborupdate
        if self.num_recheck is None:
            self.num_recheck = len(self.neighbors)-1
        if num_interest > len(self.neighbors):
            num_interest = len(self.neighbors)
        # Find the greatest neighbor distance which could have become the new smallest.
        while self.neighbors[self.num_recheck][1]-self.neighbors[num_interest-1][1] < 2*(self.self_dist_since_update + self.neighbor_dist_since_update):
            self.num_recheck += 1
            # If this has gotten too large, we need to recalculate the distances.
            if self.num_recheck > 10:
                self.num_recheck = min(num_interest, len(self.neighbors))
                self.self_dist_since_update = 0
                self.neighbor_dist_since_update = 0
                for neighbor in self.neighbors:
                    neighbor[1] = distance(self, neighbor[0])
                self.neighbors.sort(key=lambda x: x[1])
                break
    def closest_neighbors(self, n):
        relevant_distances = [(self.neighbors[i][0], distance(self, self.neighbors[i][0])) for i in range(self.num_recheck)]
        relevant_distances.extend([
            (LPoint(0, self.y, True), self.x*2),
            (LPoint(1, self.y, True), (1-self.x)*2),
            (LPoint(self.x, 0, True), (self.y)*2),
            (LPoint(self.x, 1, True), (1-self.y)*2)
        ])
        relevant_distances.sort(key=lambda x: x[1])
        return relevant_distances[:n]
    def tangent_neighbors(self, global_min_dist):
        TANGENCY_THRESHOLD = 1.01 # when is a neighbor tangent? versus merely nearby.
        relevant_distances = [(self.neighbors[i][0], distance(self, self.neighbors[i][0])) for i in range(len(self.neighbors))]
        relevant_distances.extend([
            (LPoint(0, self.y, True), self.x*2),
            (LPoint(1, self.y, True), (1-self.x)*2),
            (LPoint(self.x, 0, True), (self.y)*2),
            (LPoint(self.x, 1, True), (1-self.y)*2)
        ])
        relevant_distances.sort(key=lambda x: x[1])
        return [x for x in relevant_distances if x[1] < global_min_dist*TANGENCY_THRESHOLD]

    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        return "({}, {}, i={})".format(self.x, self.y, self.index)


def read_opt(fn):
    with open(fn, 'r') as f:
        lines = f.readlines()
    points = []
    for line in lines:
        _, x, y = tuple(x.strip() for x in line.strip().split())
        points.append(OPoint(float(x)+.5, float(y)+.5))
    for i in range(len(points)):
        points[i].neighbors = [[points[j], distance(points[i], points[j])] for j in range(len(points)) if j != i]
        points[i].update_distances(999, 999) # force recalc
    return points

def save_points(points, fn):
    with open(fn, 'w') as f:
        for i, p in enumerate(points):
            f.write("{} {} {}\n".format(i, p.x-0.5, p.y-0.5))

def copy_points(points):
    newpoints = [OPoint(p.x, p.y) for p in points]
    for i in range(len(newpoints)):
        newpoints[i].neighbors = [newpoints[j] for j in range(len(newpoints)) if j != i]
    return newpoints

if __name__ == "__main__":
    points = [OPoint(x/10, y/10) for x in range(1,10) for y in range(1,10)]
    display.draw_circles(points, max_radius(points)/2)
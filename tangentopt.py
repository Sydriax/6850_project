from geom import *
from math import *
from display import *
from initialize import *
from tqdm import tqdm
from forceopt import get_linear_force, opt

def project_parabola(p, wall, r):
    # print('entered parabola land for wall', wall)
    WP = 1
    # first, determine which wall is in question. Assumes the first point is a wall.
    if (wall.x == 0):
        dy2 = (r.y - p.y)*(r.y - p.y) # vertical sq distance between points
        rx2 = r.x*r.x
        # want x s.t. 2*sqrt(x^2) = sqrt(dy^2 + (r.x-x)^2)
        return (WP*(sqrt(3*dy2 + 4*rx2) - r.x) / 3, p.y)
    elif (wall.y == 0):
        dx2 = (r.x - p.x)*(r.x - p.x) # horizontal sq distance between points
        ry2 = r.y*r.y
        # want y s.t. 2*sqrt(y^2) = sqrt(dx^2 + (r.y-y)^2)
        return (p.x, WP*(sqrt(3*dx2 + 4*ry2) - r.y) / 3)
    elif (wall.x == 1):
        dy2 = (r.y - p.y)*(r.y - p.y)
        rx2 = (1-r.x)*(1-r.x)
        return (1 - (WP*(sqrt(3*dy2 + 4*rx2) - (1-r.x)) / 3), p.y)
    elif (wall.y == 1):
        dx2 = (r.x - p.x)*(r.x - p.x)
        ry2 = (1-r.y)*(1-r.y)
        return (p.x, 1 - (WP*(sqrt(3*dx2 + 4*ry2) - (1-r.y)) / 3))
    else:
        raise Exception("Unknown wall")

def project_walls(p, wall1, wall2):
    if wall1.y == 0 or wall1.y == 1:
        # swap walls -- we want the vertical wall to be first
        wall1, wall2 = wall2, wall1
    # now actually solve
    if wall1.x == wall2.y:
        # we want to project onto y=x
        d = (p.x+p.y)/2
        return (d, d)
    else:
        # we want to project onto y=1-x
        if wall1.x == 0:
            d = (p.x+1-p.y)/2
            return (d, 1-d)
        else:
            d = (p.y+1-p.x)/2
            return (1-d, d)

def project_bisector(p, r1, r2):
    # print('Entered bisector land')
    linevec = (r2.y-r1.y, r1.x-r2.x) # vector along the line
    # normalize normal
    l2_norm = sqrt(linevec[0]**2 + linevec[1]**2)
    linevec = (linevec[0]/l2_norm, linevec[1]/l2_norm)
    # find a point on the line -- in this case, the midpoint of r1, r2
    midpoint = ((r1.x+r2.x)/2, (r1.y+r2.y)/2)
    # this is our new basis, so shift the point to the origin
    np = (p.x-midpoint[0], p.y-midpoint[1])
    # project the point onto the line
    proj = np[0]*linevec[0] + np[1]*linevec[1]
    # shift back to the original position
    return (proj*linevec[0]+midpoint[0], proj*linevec[1]+midpoint[1])

def project_tangency(p, r1, r2):
    if r1.is_wall and r2.is_wall:
        ans = project_walls(p, r1, r2)
    elif r1.is_wall:
        ans = project_parabola(p, r1, r2)
    elif r2.is_wall:
        ans = project_parabola(p, r2, r1)
    else:
        ans = project_bisector(p, r1, r2)
    # print(r1.is_wall, r2.is_wall, ans)
    return ans


# idea: repeatedly move points partways towards perpendicular bisector
# of their two nearest neighbors.
def get_force(p, dt, min_dist):
    tangents = [t[0] for t in p.tangent_neighbors(min_dist)]
    lt = len(tangents)
    # print(lt)
    if lt > 6:
        raise Exception("Too many tangent neighbors")
    if lt <= 1:
        return get_linear_force(p, dt*min_dist/100)
        # return (0,0)
    elif lt == 2:
        np = project_tangency(p, tangents[0], tangents[1])
    elif lt > 2:
        # we compute them in a ring (to save combinatorial time)
        # and take the average
        np_sum = [0,0]
        for i in range(lt-1):
            for j in range(i+1, lt):
                np = project_tangency(p, tangents[i], tangents[j])
                np_sum[0] += np[0]
                np_sum[1] += np[1]
        norm = lt*(lt-1)/2
        np = (np_sum[0]/norm, np_sum[1]/norm)
        # print('final:', np)
    # now we have a target, so we need to move it
    return ((p.x-np[0])*dt, (p.y-np[1])*dt)

def update_points_tangent(points, forces, addjitter):
    distmoved = []
    # update the points based on the forces
    for i in range(len(points)):
        old = LPoint(points[i].x, points[i].y)
        points[i].x -= forces[i][0]
        points[i].y -= forces[i][1]
        # add jitter proportional to dt to prevent points from being stuck on top of each other
        if addjitter:
            add_jitter(points[i], 0.0001)
        fix_bounds(points[i])
        distmoved.append(distance(old, points[i]))
    # update the distances
    maxdistmoved = max(distmoved)
    # print('regopt:', maxdistmoved)
    for i in range(len(points)):
        points[i].update_distances(distmoved[i], maxdistmoved, 6)
    return points

def tangent_opt(points, dt=0.1, niter=800, dt_ratio=1):
    # iterate niter times
    for i in tqdm(range(niter)):
        min_dist = max_radius(points)*2
        forces = []
        for p in points:
            # print('starting point at (', p.x, p.y, ')')
            forces.append(get_force(p, dt, min_dist))
            # print('force is', forces[-1])
        addjitter = i < niter//2
        points = update_points_tangent(points, forces, addjitter=addjitter)
        # decrease dt
        dt *= dt_ratio
    return points

if __name__ == "__main__":
    # NPOINTS = 407
    np.random.seed(21) # reproduceable results
    NPOINTS = 55
    quality = 0
    coords = uniform_start(NPOINTS)
    #coords = [(0.12,0.12),(0.3,0.3), (0.4,0.75),(0.75,0.4)]
    points = initialize(coords)
    points = opt(points, 0.1, 800, 0.99)
    prepoints = [LPoint(p.x, p.y) for p in points]
    points = tangent_opt(points, .01, 1, 1)
    points = opt(points, 0.1*.99**800, 800, 0.99)
    circwin = []
    circwin.append(draw_circles(prepoints, max_radius(prepoints), 'init'))
    circwin.append(draw_circles(points, max_radius(points), 'final'))
    # optpoints = read_opt('31opt.txt')
    # circwin.append(draw_circles(optpoints, max_radius(optpoints), 'optimal'))
    print('init  radius', max_radius(prepoints))
    print('final radius', max_radius(points))
    # print('optim radius', max_radius(optpoints))
    maintain_windows(circwin)

from random import uniform
from geom import *
from tqdm import tqdm
from display import *
from initialize import *
from pprint import pprint
from math import *
import numpy as np

def get_linear_force(p, dt):
    closest2 = p.closest_neighbors(2)
    f = (closest2[0][0].x - p.x, closest2[0][0].y - p.y)
    if f[0] == 0 and f[1] == 0:
        f = (1, 1)
    fnorm = 1 / ((f[0]**2 + f[1]**2)**0.5)
    return (f[0]*fnorm*dt, f[1]*fnorm*dt)

def calculate_forces(points):
    # find the point closest to each point
    closestpoints = []
    distdiffnorms = []
    for i in range(len(points)):
        closest2 = points[i].closest_neighbors(2)
        closestpoints.append(closest2[0][0])
        distdiffnorms.append(closest2[1][1] - closest2[0][1])
    forces = [(closestpoints[i].x - points[i].x, closestpoints[i].y - points[i].y) for i in range(len(points))]
    for i in range(len(forces)):
        # normalize to strength 1
        if forces[i][0] == 0 and forces[i][1] == 0:
            forces[i] = (1, 1)
        fnorm = 1 / ((forces[i][0]**2 + forces[i][1]**2)**0.5)
        forces[i] = (forces[i][0]*fnorm, forces[i][1]*fnorm)
    return forces, distdiffnorms

def update_points(points, forces, dt, addjitter):
    distmoved = []
    # update the points based on the forces
    for i in range(len(points)):
        old = LPoint(points[i].x, points[i].y)
        points[i].x -= forces[i][0] * dt
        points[i].y -= forces[i][1] * dt
        # add jitter proportional to dt to prevent points from being stuck on top of each other
        if addjitter:
            add_jitter(points[i], dt/5)
        fix_bounds(points[i])
        distmoved.append(distance(old, points[i]))
    # update the distances
    maxdistmoved = max(distmoved)
    # print('regopt:', maxdistmoved)
    for i in range(len(points)):
        points[i].update_distances(distmoved[i], maxdistmoved, 2)
    return points

def opt(points, dt, niter, dt_ratio=.99):
    # iterate niter times
    for i in tqdm(range(niter)):
        forces, _ = calculate_forces(points)
        addjitter = (dt > 0.001/(len(points)**.5))
        # addjitter = i < niter//2
        points = update_points(points, forces, dt, addjitter=addjitter)
        # decrease dt
        dt *= dt_ratio
    return points

def final_update_points(points, forces, distdiffnorms):
    distmoved = []
    # update the points based on the forces
    for i in range(len(points)):
        old = LPoint(points[i].x, points[i].y)
        points[i].x -= forces[i][0] * distdiffnorms[i] / 10
        points[i].y -= forces[i][1] * distdiffnorms[i] / 10
        fix_bounds(points[i])
        distmoved.append(distance(old, points[i]))
    # update the distances
    maxdistmoved = max(distmoved)
    # print('finopt:', maxdistmoved)
    for i in range(len(points)):
        points[i].update_distances(distmoved[i], maxdistmoved, 2)
    return points

def final_opt(points, niter):
    # iterate niter times
    for _ in tqdm(range(niter)):
        forces, distdiffnorms = calculate_forces(points)
        points = final_update_points(points, forces, distdiffnorms)
    return points

if __name__ == "__main__":
    # NPOINTS = 407
    NPOINTS = 3
    quality = 0
    finalprepoints = None
    finalpoints = None
    # for i in tqdm(range(20)):
    #grid = hexagonal_start(NPOINTS)
    grid = uniform_start(NPOINTS)
    points = initialize(grid)
    init = [LPoint(p.x, p.y) for p in points]
    points = opt(points, .25/(NPOINTS**.5), 250, 0.99)
    prepoints = [LPoint(p.x, p.y) for p in points]
    points = opt(points, .25/(NPOINTS**.5) * 0.99**250, 2000, 0.9975)
    newquality = max_radius(points)
    if newquality >= quality:
        quality = newquality
        finalprepoints = prepoints
        finalpoints = points
    # print('circles:')
    # pprint(finalpoints)
    print('radius', max_radius(finalpoints))
    circwin = []
    circwin.append(draw_circles(init, max_radius(init), 'init'))
    circwin.append(draw_circles(finalprepoints, max_radius(finalprepoints), 'hot'))
    circwin.append(draw_circles(finalpoints, max_radius(finalpoints), 'cold'))
    #optpoints = read_opt(str(NPOINTS)+'opt.txt')
    #circwin.append(draw_circles(optpoints, max_radius(optpoints), 'optimal'))
    maintain_windows(circwin)

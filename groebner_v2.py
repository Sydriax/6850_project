from numbers import Rational
from geom import *
from math import *
from display import *
from initialize import *
from tqdm import tqdm
from pprint import pprint
from forceopt import opt

import sympy

def wall_polynomial(pi, wall, symbols, symrad):
    if wall.x == 0: # left vertical wall
        return symbols[pi][0] - symrad
    elif wall.x == 1: # right vertical wall
        return 1 - symbols[pi][0] - symrad
    elif wall.y == 0: # bottom horizontal wall
        return symbols[pi][1] - symrad
    elif wall.y == 1: # top horizontal wall
        return 1 - symbols[pi][1] - symrad
    raise Exception("Unknown wall")

def circle_polynomial(pi, ri, symbols, symrad):
    dx2 = (symbols[ri][0] - symbols[pi][0])**2
    dy2 = (symbols[ri][1] - symbols[pi][1])**2
    return dx2 + dy2 - (2*symrad)**2

def get_tangent_polynomial(p, r, symbols, symrad):
    if r.is_wall:
        return wall_polynomial(p.index, r, symbols, symrad)
    return circle_polynomial(p.index, r.index, symbols, symrad)

def stay_polynomials(p, symbols):
    return [
        symbols[p.index][0] - p.x,
        symbols[p.index][1] - p.y
    ]

def get_tangent_polynomials(p, min_dist, symbols, symrad):
    tangents = [t[0] for t in p.tangent_neighbors(min_dist)]
    print(p.index, 'tangent to', [t if t.is_wall else t.index for t in tangents])
    lt = len(tangents)
    # print(lt)
    if lt > 6:
        raise Exception("Too many tangent neighbors")
    if lt < 1:
        return  [] #stay_polynomials(p, symbols)
    constraints = []
    for i in range(lt):
        constraints.append(get_tangent_polynomial(p, tangents[i], symbols, symrad))
    return constraints

def detect_symmetries(points, symbols, symrad):
    EPSILON = 0.001
    constraints = []
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            if abs(points[i].x+points[j].x-1) < EPSILON:
                print('X symmetry detected:', i, j)
                constraints.append(symbols[i][0] + symbols[j][0] - 1)
            if abs(points[i].y+points[j].y-1) < EPSILON:
                print('Y symmetry detected:', i, j)
                constraints.append(symbols[i][1] + symbols[j][1] - 1)
            if abs(points[i].x-points[j].x) < EPSILON:
                print('X alignment detected:', i, j)
                constraints.append(symbols[i][0] - symbols[j][0])
            if abs(points[i].y-points[j].y) < EPSILON:
                print('Y alignment detected:', i, j)
                constraints.append(symbols[i][1] - symbols[j][1])
    return constraints

def get_all_constraints(points, symbols, symrad):
    constraints = []
    min_dist = max_radius(points)*2
    for p in points:
        constraints.extend(get_tangent_polynomials(p, min_dist, symbols, symrad))
    return constraints

def initialize_symbols(points):
    s = ' '.join(['x{} y{}'.format(i, i) for i in range(len(points))])
    symbols = sympy.symbols(s, postitive=True)
    symrad = sympy.symbols('rad', positive=True)
    return ([[symbols[i], symbols[i+1]] for i in range(0, len(points)*2, 2)], symrad)

def eval_solution(coords):
    for i in range(len(coords)):
        if coords[i].x < 0 or coords[i].x > 1:
            return None
        if coords[i].y < 0 or coords[i].y > 1:
            return None
    return max_radius(coords)

def solve_groebner(points):
    symbols, symrad = initialize_symbols(points)
    print(symbols, symrad)
    # update all point distances. Here, we can afford an expensive operation.
    for p in points:
        p.update_distances(999, 999, 6)
    constraints = get_all_constraints(points, symbols, symrad)
    constraints.extend(detect_symmetries(points, symbols, symrad))
    pprint(constraints)
    xsym = [x[0] for x in symbols]
    ysym = [x[1] for x in symbols]
    # poly = sympy.Poly(constraints, *xsym, *ysym)
    g = sympy.groebner(constraints, *xsym, *ysym, symrad, order='lex')
    print('Groebner', g)
    print('Starting solve')
    solutions = sympy.solve(g, *(xsym + ysym), symrad, Rational=False)
    print('Finished solve!')
    print(len(solutions), 'solutions:', solutions)
    current_best_rad = max_radius(points)
    for soln in solutions:
        newcoords = soln[:-1]
        # rad = soln[-1]
        xsol = newcoords[:len(xsym)]
        ysol = newcoords[len(xsym):]
        proposal = [LPoint(sympy.N(x), sympy.N(y)) if len(x.free_symbols) == 0 and len(y.free_symbols) == 0 else LPoint(p.x, p.y) for x, y, p in zip(xsol, ysol, points)]
        proprad = eval_solution(proposal)
        if proprad is None:
            print('Invalid solution :(')
            continue
        elif proprad >= current_best_rad:
            print('New best solution!')
            current_best_rad = proprad
            for i in range(len(points)):
                points[i].x = proposal[i].x
                points[i].y = proposal[i].y
        else:
            print('Solution is valid not better than current best')
    print('Finished!')
    pprint(points)
    return points

if __name__ == "__main__":
    # NPOINTS = 407
    np.random.seed(20) # reproduceable results
    NPOINTS = 7
    circwin = []
    while True:
        coords = uniform_start(NPOINTS)
        points = initialize(coords)
        points = opt(points, 0.1, 2000, 0.9925)
        print(max_radius(points))
        s = input()
        if len(s) > 0:
            break
    prepoints = [LPoint(p.x, p.y) for p in points]
    circwin.append(draw_circles(prepoints, max_radius(prepoints), 'init'))
    circwin[0].update()
    print('points')
    pprint(points)
    points = solve_groebner(points)
    circwin.append(draw_circles(points, max_radius(points), 'final'))
    # optpoints = read_opt('31opt.txt')
    # circwin.append(draw_circles(optpoints, max_radius(optpoints), 'optimal'))
    print('init  radius', max_radius(prepoints))
    print('final radius', max_radius(points))
    # print('optim radius', max_radius(optpoints))
    maintain_windows(circwin)


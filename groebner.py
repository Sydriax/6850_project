from geom import *
from math import *
from display import *
from initialize import *
from tqdm import tqdm
from pprint import pprint
from forceopt import opt

import sympy

def wall_polynomial(pi, wall1, wall2, symbols):
    assert(wall1.x == 0 or wall1.x == 1) # wall1 must be vertical.
    if wall1.x == 0: # left vertical wall
        if wall2.y == 0: # bottom horizontal wall
            return symbols[pi][0] - symbols[pi][1] # enforce x-y=0
        else: # top horizontal wall
            return symbols[pi][0] + symbols[pi][1] - 1 # enforce x+y-1=0
    else:
        if wall2.y == 0:
            return symbols[pi][0] + symbols[pi][1] - 1 # enforce x+y-1=0
        else:
            return symbols[pi][0] - symbols[pi][1] # enforce x-y=0

def parabola_polynomial(pi, wall, ri, symbols):
    dx2 = symbols[ri][0]**2 + symbols[pi][0]**2 - 2*symbols[ri][0]*symbols[pi][0]
    dy2 = symbols[ri][1]**2 + symbols[pi][1]**2 - 2*symbols[ri][1]*symbols[pi][1]
    if wall.x == 0:
        return 4*symbols[pi][0]**2 - dy2 - dx2
    elif wall.y == 0:
        return 4*symbols[pi][1]**2 - dy2 - dx2
    elif wall.x == 1:
        return 4-8*symbols[pi][0]+4*symbols[pi][0]**2 - dy2 - dx2
    elif wall.y == 1:
        return 4-8*symbols[pi][1]+4*symbols[pi][1]**2 - dy2 - dx2

def bisector_polynomial(pi, ri1, ri2, symbols):
    dx2_1 = symbols[ri1][0]**2 + symbols[pi][0]**2 - 2*symbols[ri1][0]*symbols[pi][0]
    dy2_1 = symbols[ri1][1]**2 + symbols[pi][1]**2 - 2*symbols[ri1][1]*symbols[pi][1]
    dx2_2 = symbols[ri2][0]**2 + symbols[pi][0]**2 - 2*symbols[ri2][0]*symbols[pi][0]
    dy2_2 = symbols[ri2][1]**2 + symbols[pi][1]**2 - 2*symbols[ri2][1]*symbols[pi][1]
    return dx2_1 - dx2_2 + dy2_1 - dy2_2

def get_tangent_polynomial(p, r1, r2, symbols):
    if r1.is_wall and r2.is_wall:
        if r1.x == 0 or r1.x == 1:
            ans = wall_polynomial(p.index, r1, r2, symbols)
        else:
            ans = wall_polynomial(p.index, r2, r1, symbols)
    elif r1.is_wall:
        ans = parabola_polynomial(p.index, r1, r2.index, symbols)
    elif r2.is_wall:
        ans = parabola_polynomial(p.index, r2, r1.index, symbols)
    else:
        ans = bisector_polynomial(p.index, r1.index, r2.index, symbols)
    # print(r1.is_wall, r2.is_wall, ans)
    return ans

def stay_polynomials(p, symbols):
    return [
        symbols[p.index][0] - p.x,
        symbols[p.index][1] - p.y
    ]

def get_tangent_polynomials(p, min_dist, symbols):
    tangents = [t[0] for t in p.tangent_neighbors(min_dist)]
    print(p, 'tangent to', tangents)
    lt = len(tangents)
    # print(lt)
    if lt > 6:
        raise Exception("Too many tangent neighbors")
    if lt <= 1:
        return stay_polynomials(p, symbols)
    elif lt == 2:
        return [get_tangent_polynomial(p, tangents[0], tangents[1], symbols)]
    elif lt > 2:
        # we compute them all
        poly = []
        for i in range(1, lt):
            poly.append(get_tangent_polynomial(p, tangents[0], tangents[i], symbols))
        return poly

def get_all_constraints(points, symbols):
    constraints = []
    min_dist = max_radius(points)*2
    for p in points:
        constraints.extend(get_tangent_polynomials(p, min_dist, symbols))
    return constraints

def initialize_symbols(points):
    s = ' '.join(['x{} y{}'.format(i, i) for i in range(len(points))])
    symbols = sympy.symbols(s, postitive=True)
    return [[symbols[i], symbols[i+1]] for i in range(0, len(points)*2, 2)]

def solve_groebner(points):
    symbols = initialize_symbols(points)
    print(symbols)
    # update all point distances. Here, we can afford an expensive operation.
    for p in points:
        p.update_distances(999, 999, 6)
    constraints = get_all_constraints(points, symbols)
    pprint(constraints)
    xsym = [x[0] for x in symbols]
    ysym = [x[1] for x in symbols]
    # poly = sympy.Poly(constraints, *xsym, *ysym)
    print('Starting solve')
    solutions = sympy.solve(constraints, *(xsym + ysym))[0]
    print('Finished solve!')
    print('Solutions:', solutions)
    xsol = solutions[:len(xsym)]
    ysol = solutions[len(xsym):]
    for i in range(len(points)):
        points[i].x = sympy.N(xsol[i])
        points[i].y = sympy.N(ysol[i])
    return points

if __name__ == "__main__":
    # NPOINTS = 407
    np.random.seed(42) # reproduceable results
    NPOINTS = 5
    quality = 0
    coords = uniform_start(NPOINTS)
    points = initialize(coords)
    points = opt(points, 0.1, 1000, 0.99)
    prepoints = [LPoint(p.x, p.y) for p in points]
    print('points', points)
    points = solve_groebner(points)
    circwin = []
    circwin.append(draw_circles(prepoints, max_radius(prepoints), 'init'))
    circwin.append(draw_circles(points, max_radius(points), 'final'))
    # optpoints = read_opt('31opt.txt')
    # circwin.append(draw_circles(optpoints, max_radius(optpoints), 'optimal'))
    print('init  radius', max_radius(prepoints))
    print('final radius', max_radius(points))
    # print('optim radius', max_radius(optpoints))
    maintain_windows(circwin)


from forceopt import *
from display import *
from initialize import *

def small_uniform_start(n, s):
    return np.random.random((n,2))*s + (1-s)/2

def close_resample_start(n, md):
    coords = [(np.random.random(),np.random.random())]
    while len(coords) < n:
        resample = True
        while resample:
            x = np.random.random()
            y = np.random.random()
            min_dist = min([distance(LPoint(x,y), LPoint(c[0],c[1])) for c in coords])
            min_dist = min(min_dist, walldist(LPoint(x,y)))
            if min_dist > md:
                resample = False
        coords.append((x,y))
    return coords

def truncated_normal_start(n):
    coords = []
    while len(coords) < n:
        x = np.random.normal()/5+0.5
        y = np.random.normal()/5+0.5
        while abs(x-0.5) >= 0.5 or abs(y-0.5) >= 0.5:
            x = np.random.normal()+0.5
            y = np.random.normal()+0.5
        coords.append((x,y))
    return coords

if __name__ == "__main__":
    NPOINTS = 55
    # coords = small_uniform_start(NPOINTS, 0.1)
    # coords = hexagonal_start(NPOINTS)
    # coords = close_resample_start(NPOINTS, 0.11)
    coords = truncated_normal_start(NPOINTS)
    points = initialize(coords)
    prepoints = [LPoint(p.x, p.y) for p in points]
    points = opt(points, .0001, 1000, dt_ratio=1)
    points = opt(points, .25/(NPOINTS**.5), 4000, dt_ratio=0.995)
    windows = []
    windows.append(draw_circles(prepoints, max_radius(prepoints), 'pre'))
    windows.append(draw_circles(points, max_radius(points), 'final'))
    print('original radius', max_radius(prepoints))
    print('our      radius', max_radius(points))
    print('ratio: {:.6f}%'.format(100*max_radius(points)/max_radius(prepoints)))
    maintain_windows(windows)
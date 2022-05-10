from forceopt import *
from display import *
from tangentopt import tangent_opt

if __name__ == "__main__":
    np.random.seed(42) # reproduceable results
    NPOINTS = 55
    points = read_opt(str(NPOINTS)+'opt.txt')
    prepoints = [LPoint(p.x, p.y) for p in points]
    if os.path.exists(str(NPOINTS)+'jitteropt.txt'):
        points = read_opt(str(NPOINTS)+'jitteropt.txt')
    else:
        for p in points:
            add_jitter(p, 0.002)
        # jitteredpoints = [LPoint(p.x, p.y) for p in points]
        points = opt(points, .0001, 2000, dt_ratio=0.9975)
        save_points(points, str(NPOINTS)+'jitteropt.txt')
    points = tangent_opt(points, 0.1, 0)
    windows = []
    windows.append(draw_circles(prepoints, max_radius(prepoints), 'opt'))
    # windows.append(draw_circles(jitteredpoints, max_radius(jitteredpoints), 'jittered'))
    windows.append(draw_circles(points, max_radius(points), 'final'))
    # print('jittered radius', max_radius(jitteredpoints))
    print('original radius', max_radius(prepoints))
    print('our      radius', max_radius(points))
    print('ratio: {:.6f}%'.format(100*max_radius(points)/max_radius(prepoints)))
    maintain_windows(windows)
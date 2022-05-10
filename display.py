from graphics import *

WINSIZE = 1000
MARGINSIZE = 100

def transform_coords(x, y):
    # x and y are between 0 and 1
    x = x * (WINSIZE - 2*MARGINSIZE) + MARGINSIZE
    y = y * (WINSIZE - 2*MARGINSIZE) + MARGINSIZE
    return (x, y)


def draw_circles(circlecoords, radius, label=None):
    TEXTSIZE=36
    win = GraphWin("My Packing", WINSIZE, WINSIZE, autoflush=False)
    win.setBackground("white")
    bounds = Rectangle(Point(MARGINSIZE, MARGINSIZE), Point(WINSIZE-MARGINSIZE, WINSIZE-MARGINSIZE))
    bounds.draw(win)
    for circle in circlecoords:
        c = Circle(Point(*transform_coords(circle.x, circle.y)), radius*(WINSIZE - 2*MARGINSIZE))
        c.setOutline("red")
        c.draw(win)
    t = Text(Point(WINSIZE/2,WINSIZE-MARGINSIZE/2), "Radius: {}".format(radius))
    t.setSize(TEXTSIZE)
    t.draw(win)
    if label is not None:
        t = Text(Point(WINSIZE/2,MARGINSIZE/2), label)
        t.setSize(TEXTSIZE)
        t.draw(win)
    win.flush()
    return win

def maintain_windows(windows):
    wincopy = windows[:]
    while len(wincopy) > 0:
        for win in wincopy:
            win.update()
        wincopy = [w for w in wincopy if not w.isClosed()]
        time.sleep(.1)

if __name__ == "__main__":
    from collections import namedtuple
    OCircle = namedtuple('Circle', ['x', 'y'])
    draw_circles([OCircle(.25, .25), OCircle(.5, .5), OCircle(.75, .75)], .1)
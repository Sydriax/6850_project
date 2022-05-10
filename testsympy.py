import sympy

x, y = sympy.symbols('x y')
s = sympy.solve(x**2 + y**2 - 1, x, y)
print(s)
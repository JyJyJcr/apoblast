import sympy
from apoblast import Model, Library, Transformation, Constrainer, collect_follower
from IPython.display import display
from utils import listorder

import logging

logging.basicConfig(level=logging.DEBUG)

# define coordinate
x = sympy.Symbol("x")
y = sympy.Symbol("y")
z = sympy.Symbol("z")

# define field
ux = sympy.Function("u_x")
uy = sympy.Function("u_y")
uz = sympy.Function("u_z")

# define model
m = Model(
    [
        x,
        y,
        z,
    ],
    [
        ux,
        uy,
        uz,
    ],
)

# define library
monoterms = [
    (ux(x, y, z), 0),
    (uy(x, y, z), 0),
    (uz(x, y, z), 0),
    (sympy.diff(ux(x, y, z), x), 1),
    (sympy.diff(ux(x, y, z), y), 1),
    (sympy.diff(ux(x, y, z), z), 1),
    (sympy.diff(uy(x, y, z), x), 1),
    (sympy.diff(uy(x, y, z), y), 1),
    (sympy.diff(uy(x, y, z), z), 1),
    (sympy.diff(uz(x, y, z), x), 1),
    (sympy.diff(uz(x, y, z), y), 1),
    (sympy.diff(uz(x, y, z), z), 1),
    (sympy.diff(ux(x, y, z), x, 2), 2),
    (sympy.diff(ux(x, y, z), x, y), 2),
    (sympy.diff(ux(x, y, z), x, z), 2),
    (sympy.diff(ux(x, y, z), y, 2), 2),
    (sympy.diff(ux(x, y, z), y, z), 2),
    (sympy.diff(ux(x, y, z), z, 2), 2),
    (sympy.diff(uy(x, y, z), x, 2), 2),
    (sympy.diff(uy(x, y, z), x, y), 2),
    (sympy.diff(uy(x, y, z), x, z), 2),
    (sympy.diff(uy(x, y, z), y, 2), 2),
    (sympy.diff(uy(x, y, z), y, z), 2),
    (sympy.diff(uy(x, y, z), z, 2), 2),
    (sympy.diff(uz(x, y, z), x, 2), 2),
    (sympy.diff(uz(x, y, z), x, y), 2),
    (sympy.diff(uz(x, y, z), x, z), 2),
    (sympy.diff(uz(x, y, z), y, 2), 2),
    (sympy.diff(uz(x, y, z), y, z), 2),
    (sympy.diff(uz(x, y, z), z, 2), 2),
]
# r => s r
# t => s^zet t
# u => s^chi u
zet = 2
chi = -0.5
l = Library(
    m, listorder(monoterms, 0, 0, 4, 3, lambda o, d: o * chi - d >= 3 * chi - 2)
)

# define transformations

# x reflection
trrefl = Transformation(
    m,
    (
        -x,
        y,
        z,
    ),
    (-ux(x, y, z), uy(x, y, z), uz(x, y, z)),
)

# rotations
th = sympy.Symbol("th")
# x rotation
trxrot = Transformation(
    m,
    (x, sympy.cos(th) * y + sympy.sin(th) * z, sympy.cos(th) * z - sympy.sin(th) * y),
    (
        ux(x, y, z),
        sympy.cos(th) * uy(x, y, z) + sympy.sin(th) * uz(x, y, z),
        sympy.cos(th) * uz(x, y, z) - sympy.sin(th) * uy(x, y, z),
    ),
    parameter=[th],
)
# y rotation
tryrot = Transformation(
    m,
    (
        sympy.cos(th) * x - sympy.sin(th) * z,
        y,
        sympy.cos(th) * z + sympy.sin(th) * x,
    ),
    (
        sympy.cos(th) * ux(x, y, z) - sympy.sin(th) * uz(x, y, z),
        uy(x, y, z),
        sympy.cos(th) * uz(x, y, z) + sympy.sin(th) * ux(x, y, z),
    ),
    parameter=[th],
)
# z rotation
trzrot = Transformation(
    m,
    (sympy.cos(th) * x + sympy.sin(th) * y, sympy.cos(th) * y - sympy.sin(th) * x, z),
    (
        sympy.cos(th) * ux(x, y, z) + sympy.sin(th) * uy(x, y, z),
        sympy.cos(th) * uy(x, y, z) - sympy.sin(th) * ux(x, y, z),
        uz(x, y, z),
    ),
    parameter=[th],
)

# define LHS
lhs = [ux(x, y, z), uy(x, y, z), uz(x, y, z)]
# define LHS library
lhs_library = Library(m, (ux(x, y, z), uy(x, y, z), uz(x, y, z)))

# for invariant terms
# lhs=[sympy.Number(1)]
# lhs_library=Library(m,(sympy.Number(1),))

# execute alogirithm
follower = collect_follower(
    lhs,
    l,
    # reflection
    Constrainer(trrefl),
    # rotation
    Constrainer(
        trxrot,
        parameter_fix=[
            0,
        ],
        differential=th,
    ),
    Constrainer(
        tryrot,
        parameter_fix=[
            0,
        ],
        differential=th,
    ),
    Constrainer(
        trzrot,
        parameter_fix=[
            0,
        ],
        differential=th,
    ),
    lhs_library=lhs_library,
)

# display results
# display(follower_mat)
for f in follower:
    display(sympy.Matrix(f))

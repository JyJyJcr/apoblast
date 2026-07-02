import sympy
from apoblast import Model, Library, Transformation, Constrainer, collect_follower
from IPython.display import display

import logging

logging.basicConfig(level=logging.INFO)

# define coordinate
x = sympy.Symbol("x")
y = sympy.Symbol("y")

# define field (no field in this case)
# u=sympy.Function('u')
# v=sympy.Function('v')

# define model
m = Model(
    [x, y],
    [
        # u,
        # v
    ],
)

# define library
monoterms = [x, y]


def listorder(monoterms, order):
    if len(monoterms) == 0:
        return [1]
    else:
        arr = []
        for i in range(0, order + 1):
            arr += map(
                lambda t: (monoterms[0] ** i) * t, listorder(monoterms[1:], order - i)
            )
        return arr


l = Library(m, listorder(monoterms, 10))

# define transformations

# reflections
# x reflection
trrefl = Transformation(
    m,
    (-x, y),
    (
        # -u(x,y),
        # v(x,y)
    ),
)
# y reflection
trrefl2 = Transformation(
    m,
    (x, -y),
    (
        # u(x,y),
        # -v(x,y)
    ),
)
# 45 degree rotation
trrefl3 = Transformation(
    m,
    (y, x),
    (
        # v(x,y),
        # u(x,y)
    ),
)


# rotation
def trrot(th):
    return Transformation(
        m,
        (sympy.cos(th) * x + sympy.sin(th) * y, sympy.cos(th) * y - sympy.sin(th) * x),
        (
            #     sympy.cos(th)*u(x,y)+sympy.sin(th)*v(x,y),
            #     sympy.cos(th)*v(x,y)-sympy.sin(th)*u(x,y)
        ),
    )


# rotation with exponential form, which results in a longer computation time. do not use.
def expsin(th):
    return (sympy.exp(sympy.I * th) - sympy.exp(-sympy.I * th)) / 2 * sympy.I


def expcos(th):
    return (sympy.exp(sympy.I * th) + sympy.exp(-sympy.I * th)) / 2


def exptrrot(th):
    return Transformation(
        m,
        (expcos(th) * x + expsin(th) * y, expcos(th) * y - expsin(th) * x),
        (
            #      expcos(th)*u(x,y)+expsin(th)*v(x,y),
            #      expcos(th)*v(x,y)-expsin(th)*u(x,y)
        ),
    )


# modern way to define rotation
th = sympy.Symbol("th")
trrot_free = Transformation(
    m,
    [sympy.cos(th) * x + sympy.sin(th) * y, sympy.cos(th) * y - sympy.sin(th) * x],
    [
        #      sympy.cos(th)*u(x,y)+sympy.sin(th)*v(x,y),
        #      sympy.cos(th)*v(x,y)-sympy.sin(th)*u(x,y)
    ],
    parameter=[th],
)

# define LHS
lhs = [sympy.Number(1)]
# define LHS library
lhs_library = Library(m, [sympy.Number(1)])

# execute alogirithm
follower = collect_follower(
    lhs,
    l,
    # Constrainer(trrefl),
    Constrainer(
        trrefl2
    ),  # notice: with 5-fold rotational symmetry, trrefl and trrefl2 are not same.
    # Constrainer(trrefl3),
    Constrainer(trrot(sympy.pi / 5 * 2)),
    # Constrainer(trrot_free,parameter_fix=[sympy.pi/5*2,]), # modern way
    # Constrainer(trrot(1)), # old way to express cotinuous rotation
    lhs_library=lhs_library,
)

# display results
# display(follower_mat)
for f in follower:
    display(sympy.Matrix(f))

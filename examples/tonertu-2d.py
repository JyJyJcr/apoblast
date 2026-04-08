import sympy
from apoblast import Model, Library, Transformation,Constrainer, collect_follower, collect_follower_raw
from IPython.display import display
from utils import listdorder, listorder

# define coordinate
x=sympy.Symbol('x')
y=sympy.Symbol('y')

#define field
ux=sympy.Function('u_x')
uy=sympy.Function('u_y')

# define model
m=Model(
    [
        x,
        y,
    ],[
        ux,
        uy,
    ]
)

# define library
monoterms=[
    (ux(x,y),0),

    (uy(x,y),0),

    (sympy.diff(ux(x,y),x),1),
    (sympy.diff(ux(x,y),y),1),

    (sympy.diff(uy(x,y),x),1),
    (sympy.diff(uy(x,y),y),1),

    (sympy.diff(ux(x,y),x,2),2),
    (sympy.diff(ux(x,y),x,y),2),
    (sympy.diff(ux(x,y),y,2),2),

    (sympy.diff(uy(x,y),x,2),2),
    (sympy.diff(uy(x,y),x,y),2),
    (sympy.diff(uy(x,y),y,2),2),
]
# r => s r
# t => s^zet t
# u => s^chi u
d=2
zet=2
chi=0
l=Library(m,listorder(monoterms,0,0,3,2,lambda o,d: o*chi-d >= chi-2))

# define transformations

# x reflection
trrefl=Transformation(m,
    (
        -x,
        y,
    ),(
        -ux(x,y),
        uy(x,y),
    )
)

# rotation
th=sympy.Symbol('th')
trrot=Transformation(m,
    (
        sympy.cos(th)*x+sympy.sin(th)*y,
        sympy.cos(th)*y-sympy.sin(th)*x,
    ),(
        sympy.cos(th)*ux(x,y)+sympy.sin(th)*uy(x,y),
        sympy.cos(th)*uy(x,y)-sympy.sin(th)*ux(x,y),
    ),
    parameter=[th]
)

# field only rotation
trfrot=Transformation(m,
    (
        x,
        y,
    ),(
        sympy.cos(th)*ux(x,y)+sympy.sin(th)*uy(x,y),
        sympy.cos(th)*uy(x,y)-sympy.sin(th)*ux(x,y),
    ),
    parameter=[th]
)

# define LHS
lhs=[ux(x,y),uy(x,y)]
# define LHS library
lhs_library=Library(m,(ux(x,y),uy(x,y)))

# for invariant terms
#lhs=[sympy.Number(1)]
#lhs_library=Library(m,(sympy.Number(1),))

# execute alogirithm
follower=collect_follower(lhs,
    l,
    # reflection
    Constrainer(trrefl),
    # rotation
    Constrainer(trrot,parameter_fix=[0,],differential=th),
    #Constrainer(trfrot,parameter_fix=[0,],differential=th),
    lhs_library=lhs_library,
    print_progress=True)

# display results
# display(follower_mat)
for f in follower:
    display(sympy.Matrix(f))

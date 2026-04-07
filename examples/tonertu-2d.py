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


# 110
#print(len(listorder(monoterms,0,0,3,1)))
# 0
#print(len(listorder(monoterms,0,0,3,1,filter=lambda o,d: False)))

# r => s r
# t => s^zet t
# u => s^chi u
d=2
zet=2
chi=0
l=Library(m,listorder(monoterms,0,0,3,2,lambda o,d: o*chi-d >= chi-2))


# define transformation

# y axis reflection
trrefl=Transformation(m,
    (
        -x,
        y,
    ),(
        -ux(x,y),
        uy(x,y),
    )
)

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

lhs=[ux(x,y),uy(x,y)]
lhs_library=Library(m,(ux(x,y),uy(x,y)))
#lhs_library=Library(m,(u(x,y),v(x,y)))

#lhs=[sympy.Number(1)]
#lhs_library=Library(m,(sympy.Number(1),))

follower=collect_follower(lhs,
    l,
    Constrainer(trrefl),
    Constrainer(trrot,parameter_fix=[0,],differential=th),
    #Constrainer(trfrot,parameter_fix=[0,],differential=th),
    #exptrrot(1),
    lhs_library=lhs_library,
    print_progress=True)
# print("follower")
# display(follower_mat)
for f in follower:
    display(sympy.Matrix(f))

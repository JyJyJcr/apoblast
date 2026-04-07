import sympy
from apoblast import Model, Library, Transformation,Constrainer, collect_follower, collect_follower_raw
from IPython.display import display
from utils import listdorder, listorder

# define coordinate
x=sympy.Symbol('x')
y=sympy.Symbol('y')
t=sympy.Symbol('t')

#define field
h=sympy.Function('h')

# define model
m=Model(
    [
        x,
        y,
        t,
    ],[
        h,
    ]
)

# 110
#print(len(listorder(monoterms,0,0,3,1)))
# 0
#print(len(listorder(monoterms,0,0,3,1,filter=lambda o,d: False)))

# r => s r
# t => s^zet t
# u => s^chi u
# zet=2
# chi=0
monoterms=listdorder(h(x,y,t),[x,y,t],0,4)
l=Library(m,listorder(monoterms,0,0,5,4))


# define transformation

dx=sympy.Symbol('dx')
dy=sympy.Symbol('dy')
dt=sympy.Symbol('dt')

trshift=Transformation(m,
    (
        x+dx,
        y+dy,
        t+dt,
    ),(
        h(x,y,t),
    ),
    parameter=[dx,dy,dt]
)

# x axis reflection
trrefl=Transformation(m,
    (
        -x,
        y,
        t,
    ),(
        h(x,y,t),
    )
)

th=sympy.Symbol('th')


trrot=Transformation(m,
    (
        sympy.cos(th)*x+sympy.sin(th)*y,
        sympy.cos(th)*y-sympy.sin(th)*x,
        t,
    ),(
        h(x,y,t),
    ),
    parameter=[th]
)

#

dh=sympy.Symbol('dh')

trhshift=Transformation(m,
    (
        x,
        y,
        t,
    ),(
        h(x,y,t)+dh,
    ),
    parameter=[dh]
)
trhrefl=Transformation(m,
    (
        x,
        y,
        t,
    ),(
        -h(x,y,t),
    )
)

ex=sympy.Symbol('ex')
ey=sympy.Symbol('ey')
lam=sympy.Symbol('lambda')

trtilt=Transformation(m,
    (
        x-lam*ex*t,
        y-lam*ey*t,
        t,
    ),(
        h(x,y,t) + ex*x + ey*y,
    ),
    parameter=[ex,ey,lam]
)

lhs=[sympy.diff(h(x,y,t),x,2)+sympy.diff(h(x,y,t),y,2)]
lhs_library=Library(m,[sympy.diff(h(x,y,t),x,2),sympy.diff(h(x,y,t),y,2)])
print("lhs")

#lhs=[sympy.Number(1)]
#lhs_library=Library(m,(sympy.Number(1),))

follower=collect_follower(lhs,
    l,
    # reflection
    Constrainer(trrefl),
    # rotation
    Constrainer(trrot,parameter_fix=[0,],differential=th),
    # space time shift
    Constrainer(trshift,parameter_fix=[0,0,0,],differential=dx),
    #Constrainer(trshift,parameter_fix=[0,0,0,0,],differential=dy),
    #Constrainer(trshift,parameter_fix=[0,0,0,0,],differential=dz),
    Constrainer(trshift,parameter_fix=[0,0,0,],differential=dt),
    # h shift
    Constrainer(trhshift,parameter_fix=[0,],differential=dh),
    # h reflection
    # Constrainer(trhrefl),
    # tilt
    Constrainer(trtilt,parameter_fix=[0,0,1,],differential=ex),
    #Constrainer(trtilt,parameter_fix=[0,0,0,1,],differential=ey),
    #Constrainer(trtilt,parameter_fix=[0,0,0,1,],differential=ez),
    
    #exptrrot(1),
    lhs_library=lhs_library,
    print_progress=True)
# print("follower")
# display(follower_mat)
for f in follower:
    display(sympy.Matrix(f))

import sympy
from apoblast import Model, Library, Transformation,Constrainer, collect_follower, collect_follower_raw
from IPython.display import display
from utils import listdorder, listorder

# define coordinate
x=sympy.Symbol('x')
y=sympy.Symbol('y')
z=sympy.Symbol('z')
t=sympy.Symbol('t')

#define field
h=sympy.Function('h')

# define model
m=Model(
    [
        x,
        y,
        z,
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
monoterms=listdorder(h(x,y,z,t),[x,y,z,t],0,4)
l=Library(m,listorder(monoterms,0,0,5,4))


# define transformation

dx=sympy.Symbol('dx')
dy=sympy.Symbol('dy')
dz=sympy.Symbol('dz')
dt=sympy.Symbol('dt')

trshift=Transformation(m,
    (
        x+dx,
        y+dy,
        z+dz,
        t+dt,
    ),(
        h(x,y,z,t),
    ),
    parameter=[dx,dy,dz,dt]
)

# x axis reflection
trrefl=Transformation(m,
    (
        -x,
        y,
        z,
        t,
    ),(
        h(x,y,z,t),
    )
)

th=sympy.Symbol('th')

trxrot=Transformation(m,
    (
        x,
        sympy.cos(th)*y+sympy.sin(th)*z,
        sympy.cos(th)*z-sympy.sin(th)*y,
        t,
    ),(
        h(x,y,z,t),
    ),
    parameter=[th]
)

tryrot=Transformation(m,
    (
        sympy.cos(th)*x-sympy.sin(th)*z,
        y,
        sympy.cos(th)*z+sympy.sin(th)*x,
        t,
    ),(
        h(x,y,z,t),
    ),
    parameter=[th]
)

trzrot=Transformation(m,
    (
        sympy.cos(th)*x+sympy.sin(th)*y,
        sympy.cos(th)*y-sympy.sin(th)*x,
        z,
        t,
    ),(
        h(x,y,z,t),
    ),
    parameter=[th]
)

#

dh=sympy.Symbol('dh')

trhshift=Transformation(m,
    (
        x,
        y,
        z,
        t,
    ),(
        h(x,y,z,t)+dh,
    ),
    parameter=[dh]
)
trhrefl=Transformation(m,
    (
        x,
        y,
        z,
        t,
    ),(
        -h(x,y,z,t),
    )
)

ex=sympy.Symbol('ex')
ey=sympy.Symbol('ey')
ez=sympy.Symbol('ez')
lam=sympy.Symbol('lambda')

trtilt=Transformation(m,
    (
        x-lam*ex*t,
        y-lam*ey*t,
        z-lam*ez*t,
        t,
    ),(
        h(x,y,z,t) + ex*x + ey*y + ez*z,
    ),
    parameter=[ex,ey,ez,lam]
)

lhs=[sympy.diff(h(x,y,z,t),x,2)+sympy.diff(h(x,y,z,t),y,2)+sympy.diff(h(x,y,z,t),z,2)]
lhs_library=Library(m,[sympy.diff(h(x,y,z,t),x,2),sympy.diff(h(x,y,z,t),y,2),sympy.diff(h(x,y,z,t),z,2)])
print("lhs")

#lhs=[sympy.Number(1)]
#lhs_library=Library(m,(sympy.Number(1),))

follower=collect_follower(lhs,
    l,
    # reflection
    Constrainer(trrefl),
    # rotation
    Constrainer(trxrot,parameter_fix=[0,],differential=th),
    Constrainer(tryrot,parameter_fix=[0,],differential=th),
    Constrainer(trzrot,parameter_fix=[0,],differential=th),
    # space time shift
    Constrainer(trshift,parameter_fix=[0,0,0,0,],differential=dx),
    #Constrainer(trshift,parameter_fix=[0,0,0,0,],differential=dy),
    #Constrainer(trshift,parameter_fix=[0,0,0,0,],differential=dz),
    Constrainer(trshift,parameter_fix=[0,0,0,0,],differential=dt),
    # h shift
    Constrainer(trhshift,parameter_fix=[0,],differential=dh),
    # h reflection
    # Constrainer(trhrefl),
    # tilt
    Constrainer(trtilt,parameter_fix=[0,0,0,1,],differential=ex),
    #Constrainer(trtilt,parameter_fix=[0,0,0,1,],differential=ey),
    #Constrainer(trtilt,parameter_fix=[0,0,0,1,],differential=ez),
    
    #exptrrot(1),
    lhs_library=lhs_library,
    print_progress=True)
# print("follower")
# display(follower_mat)
for f in follower:
    display(sympy.Matrix(f))

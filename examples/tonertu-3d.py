import sympy
from apoblast import Model, Library, Transformation,Constrainer, collect_follower, collect_follower_raw
from IPython.display import display
from utils import listdorder, listorder

# define coordinate
x=sympy.Symbol('x')
y=sympy.Symbol('y')
z=sympy.Symbol('z')


#define field
ux=sympy.Function('u_x')
uy=sympy.Function('u_y')
uz=sympy.Function('u_z')

# define model
m=Model(
    [
        x,
        y,
        z,
    ],[
        ux,
        uy,
        uz,
    ]
)

monoterms=[
    (ux(x,y,z),0),
    (uy(x,y,z),0),
    (uz(x,y,z),0),
    (sympy.diff(ux(x,y,z),x),1),
    (sympy.diff(ux(x,y,z),y),1),
    (sympy.diff(ux(x,y,z),z),1),
    (sympy.diff(uy(x,y,z),x),1),
    (sympy.diff(uy(x,y,z),y),1),
    (sympy.diff(uy(x,y,z),z),1),
    (sympy.diff(uz(x,y,z),x),1),
    (sympy.diff(uz(x,y,z),y),1),
    (sympy.diff(uz(x,y,z),z),1),

    (sympy.diff(ux(x,y,z),x,2),2),
    (sympy.diff(ux(x,y,z),x,y),2),
    (sympy.diff(ux(x,y,z),x,z),2),
    (sympy.diff(ux(x,y,z),y,2),2),
    (sympy.diff(ux(x,y,z),y,z),2),
    (sympy.diff(ux(x,y,z),z,2),2),

    (sympy.diff(uy(x,y,z),x,2),2),
    (sympy.diff(uy(x,y,z),x,y),2),
    (sympy.diff(uy(x,y,z),x,z),2),
    (sympy.diff(uy(x,y,z),y,2),2),
    (sympy.diff(uy(x,y,z),y,z),2),
    (sympy.diff(uy(x,y,z),z,2),2),

    (sympy.diff(uz(x,y,z),x,2),2),
    (sympy.diff(uz(x,y,z),x,y),2),
    (sympy.diff(uz(x,y,z),x,z),2),
    (sympy.diff(uz(x,y,z),y,2),2),
    (sympy.diff(uz(x,y,z),y,z),2),
    (sympy.diff(uz(x,y,z),z,2),2),
]


# 110
#print(len(listorder(monoterms,0,0,3,1)))
# 0
#print(len(listorder(monoterms,0,0,3,1,filter=lambda o,d: False)))

# r => s r
# t => s^zet t
# u => s^chi u
zet=2
chi=-0.5
l=Library(m,listorder(monoterms,0,0,4,3,lambda o,d: o*chi-d >= 3*chi-2))


# define transformation

# y axis reflection
trrefl=Transformation(m,
    (
        -x,
        y,
        z,
    ),(
        -ux(x,y,z),
        uy(x,y,z),
        uz(x,y,z)
    )
)

th=sympy.Symbol('th')

trxrot=Transformation(m,
    (
        x,
        sympy.cos(th)*y+sympy.sin(th)*z,
        sympy.cos(th)*z-sympy.sin(th)*y
    ),(
        ux(x,y,z),
        sympy.cos(th)*uy(x,y,z)+sympy.sin(th)*uz(x,y,z),
        sympy.cos(th)*uz(x,y,z)-sympy.sin(th)*uy(x,y,z)
    ),
    parameter=[th]
)

tryrot=Transformation(m,
    (
        sympy.cos(th)*x-sympy.sin(th)*z,
        y,
        sympy.cos(th)*z+sympy.sin(th)*x,
    ),(
        sympy.cos(th)*ux(x,y,z)-sympy.sin(th)*uz(x,y,z),
        uy(x,y,z),
        sympy.cos(th)*uz(x,y,z)+sympy.sin(th)*ux(x,y,z)
    ),
    parameter=[th]
)

trzrot=Transformation(m,
    (
        sympy.cos(th)*x+sympy.sin(th)*y,
        sympy.cos(th)*y-sympy.sin(th)*x,
        z
    ),(
        sympy.cos(th)*ux(x,y,z)+sympy.sin(th)*uy(x,y,z),
        sympy.cos(th)*uy(x,y,z)-sympy.sin(th)*ux(x,y,z),
        uz(x,y,z)
    ),
    parameter=[th]
)

lhs=[ux(x,y,z),uy(x,y,z),uz(x,y,z)]
lhs_library=Library(m,(ux(x,y,z),uy(x,y,z),uz(x,y,z)))
#lhs_library=Library(m,(u(x,y),v(x,y)))

#lhs=[sympy.Number(1)]
#lhs_library=Library(m,(sympy.Number(1),))

follower=collect_follower(lhs,
    l,
    Constrainer(trrefl),
    #Constrainer(trxrot,parameter_fix=[sympy.pi/3,]),
    Constrainer(trxrot,parameter_fix=[0,],differential=th),
    Constrainer(tryrot,parameter_fix=[0,],differential=th),
    Constrainer(trzrot,parameter_fix=[0,],differential=th),
    #exptrrot(1),
    lhs_library=lhs_library,
    print_progress=True)
# print("follower")
# display(follower_mat)
for f in follower:
    display(sympy.Matrix(f))

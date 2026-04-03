import sympy
from apoblast import Model, Library, Transformation,Constrainer, collect_follower, collect_follower_raw
from IPython.display import display

# define coordinate
x=sympy.Symbol('x')
y=sympy.Symbol('y')


#define field
u=sympy.Function('u')
v=sympy.Function('v')

# define model
m=Model(
    [
        x,
        y
    ],[
        #u,
        #v
    ]
)

monoterms=[
    x,y,
    # u(x,y),
    # v(x,y),
    # sympy.diff(u(x,y),x),
    # sympy.diff(u(x,y),y),
    # sympy.diff(v(x,y),x),
    # sympy.diff(v(x,y),y),
    # sympy.sin(u(x,y)),
    # sympy.cos(u(x,y)),
    # sympy.sin(v(x,y)),
    # sympy.cos(v(x,y)),
    # sympy.diff(u(x,y),x,2),
    # sympy.diff(u(x,y),x,y),
    # sympy.diff(u(x,y),y,2),
    # sympy.diff(v(x,y),x,2),
    # sympy.diff(v(x,y),x,y),
    # sympy.diff(v(x,y),y,2),
]


def listorder(monoterms,order):
    if len(monoterms) == 0:
        return [1]
    else:
        arr=[]
        for i in range(0,order+1):
            arr+=map(lambda t: (monoterms[0]**i)*t,listorder(monoterms[1:],order-i))
        return arr

#print(listorder(monoterms,3))

# l=Library(m,
#     map(lambda t: t[0]*t[1],itertools.product(poly, diff))
# )
l=Library(m,listorder(monoterms,10))


# define transformation

# y axis reflection
trrefl=Transformation(m,
    (
        -x,
        y
    ),(
       # -u(x,y),
       # v(x,y)
    )
)

# x axis reflection
trrefl2=Transformation(m,
    (
        x,
        -y
    ),(
       # u(x,y),
       # -v(x,y)
    )
)
trrefl3=Transformation(m,
    (
        y,
        x
    ),(
       # v(x,y),
       # u(x,y)
    )
)


#trrefl.apply(sympy.diff(v(x,y),x))

#t=TransformationMatrix(l,trrefl)
#t.mat
# # rotation
def trrot(th):
    return Transformation(m,
        (
            sympy.cos(th)*x+sympy.sin(th)*y,
            sympy.cos(th)*y-sympy.sin(th)*x
        ),(
       #     sympy.cos(th)*u(x,y)+sympy.sin(th)*v(x,y),
       #     sympy.cos(th)*v(x,y)-sympy.sin(th)*u(x,y)
        )
    )

def expsin(th):
    return (sympy.exp(sympy.I*th)-sympy.exp(-sympy.I*th))/2*sympy.I
def expcos(th):
    return (sympy.exp(sympy.I*th)+sympy.exp(-sympy.I*th))/2

def exptrrot(th):
    return Transformation(m,
        (
            expcos(th)*x+expsin(th)*y,
            expcos(th)*y-expsin(th)*x
        ),(
      #      expcos(th)*u(x,y)+expsin(th)*v(x,y),
      #      expcos(th)*v(x,y)-expsin(th)*u(x,y)
        )
    )

# # due to the computation power limit, need to choose rotation angles carefully: they are better to be rational * pi


# lhs=[u(x,y),v(x,y)]
# lhs_library=Library(m,(u(x,y),v(x,y)))

th=sympy.Symbol("th")



lhs=[sympy.Number(1)]
lhs_library=Library(m,[sympy.Number(1)])

rot=Transformation(m,
        [
            sympy.cos(th)*x+sympy.sin(th)*y,
            sympy.cos(th)*y-sympy.sin(th)*x
        ],[
      #      sympy.cos(th)*u(x,y)+sympy.sin(th)*v(x,y),
      #      sympy.cos(th)*v(x,y)-sympy.sin(th)*u(x,y)
        ],
        parameter=[th]
    )

X=sympy.Symbol("X")
Y=sympy.Symbol("Y")

shift=Transformation(m,
        [
            x+X,
            y+Y,
        ],[
       #     u(x,y),
       #     v(x,y)
        ],
        parameter=[X,Y]
    )

follower=collect_follower(lhs,
    l,
    # Constrainer(trrefl),
    Constrainer(trrefl2),
    # Constrainer(trrefl3),
    Constrainer(trrot(sympy.pi/5*2)),
    # Constrainer(trrot(sympy.pi/4)),
    # Constrainer(trrot(sympy.pi/6)),
    # Constrainer(trrot(sympy.pi/5)),



    #Constrainer(shift,parameter_fix=[0,0],differential=X),
    #Constrainer(shift,parameter_fix=[0,0],differential=Y),
    #Constrainer(rot,parameter_fix=[0],differential=th),
    
    # trrefl3,
    # trrot(sympy.pi/2),
    # trrot(sympy.pi/4),
    # trrot(sympy.pi/6),
    # trrot(sympy.pi/5),
    #exptrrot(1),
    lhs_library=lhs_library,
    print_progress=True)
print("follower")
for f in follower:
   display(sympy.Matrix(f))

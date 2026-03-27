# all the test functions start with "test_"
# and the test files start with "test_"

import logging
import sympy
from apoblast import Model, Library, Transformation, Constrainer, ConstraintMatrix, collect_follower
from IPython.display import display

def test_model():
    # define coordinate
    t=sympy.Symbol('t')
    x=sympy.Symbol('x')
    y=sympy.Symbol('y')

    # define field
    u=sympy.Function('u')
    v=sympy.Function('v')

    # define model
    m=Model(
        [
            t,
            x,
            y
        ],[
            u,
            v
        ]
    )

    assert m.coord == (t,x,y)
    assert m.field == (u,v)
    assert m.applied_field == (u(t,x,y),v(t,x,y))

    terms = [\
        sympy.Number(1),
        x,y,
        u(t,x,y),v(t,x,y),\
        u(t,x,y)*u(t,x,y),u(t,x,y)*v(t,x,y),v(t,x,y)*v(t,x,y),\
        sympy.cos(u(t,x,y)),sympy.sin(u(t,x,y)),\
        5*u(t,x,y)+6*v(t,x,y),
        sympy.diff(u(t,x,y),t,x),
    ]

    for term in terms:
        assert m.verify(term)

def test_library():
    # define coordinate
    t=sympy.Symbol('t')
    x=sympy.Symbol('x')
    y=sympy.Symbol('y')

    # define field
    u=sympy.Function('u')
    v=sympy.Function('v')

    th=sympy.Symbol('th')

    # define model
    m=Model(
        [
            t,
            x,
            y
        ],[
            u,
            v
        ]
    )

    l=Library(m,[\
        x,y,
        u(t,x,y),v(t,x,y),\
        u(t,x,y)*u(t,x,y),u(t,x,y)*v(t,x,y),v(t,x,y)*v(t,x,y),\
        sympy.cos(u(t,x,y)),sympy.sin(u(t,x,y)),\
        5*u(t,x,y)+6*v(t,x,y),
        sympy.diff(u(t,x,y),t,x),
    ])

    arr,dic=l.classify(x+sympy.cos(th+u(t,x,y))+sympy.cos(th+v(t,x,y)))
    assert arr == [1,0,0,0,0,0,0,0,0,0,0]
    #logging.info(dic)
    assert dic == {sympy.cos(th+u(t,x,y)):1,sympy.cos(th+v(t,x,y)):1}

def test_transformation():
    # define coordinate
    t=sympy.Symbol('t')
    x=sympy.Symbol('x')
    y=sympy.Symbol('y')

    # define field
    u=sympy.Function('u')
    v=sympy.Function('v')

    # define model
    m=Model(
        [
            t,
            x,
            y
        ],[
            u,
            v
        ]
    )

    m.verify(x,with_field=False,parameter=())
    m.verify(y,with_field=False,parameter=())
    m.verify(t,with_field=False,parameter=())

    th=sympy.Symbol('th')

    tr=Transformation(m,
        [
            t,
            sympy.cos(th)*x+sympy.sin(th)*y,
            sympy.cos(th)*y-sympy.sin(th)*x
        ],[
            sympy.cos(th)*u(t,x,y)+sympy.sin(th)*v(t,x,y),
            sympy.cos(th)*v(t,x,y)-sympy.sin(th)*u(t,x,y)
        ],
        parameter=[th]
    )

    assert tr.apply(sympy.diff(u(t,x,y),x)+sympy.diff(v(t,x,y),y)).simplify() == sympy.diff(u(t,x,y),x)+sympy.diff(v(t,x,y),y)

def test_transformation_matrix():
    # define coordinate
    x=sympy.Symbol('x')
    y=sympy.Symbol('y')

    # define field
    u=sympy.Function('u')
    v=sympy.Function('v')

    # define model
    m=Model(
        [
            x,
            y
        ],[
            u,
            v
        ]
    )

    monoterms=[
        u(x,y),
        v(x,y),
        sympy.diff(u(x,y),x),
        sympy.diff(u(x,y),y),
        sympy.diff(v(x,y),x),
        sympy.diff(v(x,y),y),
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
    l=Library(m,listorder(monoterms,2))



    th=sympy.Symbol('th')

    tr=Transformation(m,
        [
            sympy.cos(th)*x+sympy.sin(th)*y,
            sympy.cos(th)*y-sympy.sin(th)*x
        ],[
            sympy.cos(th)*u(x,y)+sympy.sin(th)*v(x,y),
            sympy.cos(th)*v(x,y)-sympy.sin(th)*u(x,y)
        ],
        parameter=[th]
    )

    c=Constrainer(tr,parameter_fix=[0,],differential=th)
    t=ConstraintMatrix(l,c)
    print(t.mat)

def test_apoblast():

    # define coordinate
    x=sympy.Symbol('x')
    y=sympy.Symbol('y')
    z=sympy.Symbol('z')


    #define field
    ux=sympy.Function('ux')
    uy=sympy.Function('uy')
    uz=sympy.Function('uz')

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
        ux(x,y,z),
        uy(x,y,z),
        uz(x,y,z),
        sympy.diff(ux(x,y,z),x),
        sympy.diff(ux(x,y,z),y),
        sympy.diff(ux(x,y,z),z),
        sympy.diff(uy(x,y,z),x),
        sympy.diff(uy(x,y,z),y),
        sympy.diff(uy(x,y,z),z),
        sympy.diff(uz(x,y,z),x),
        sympy.diff(uz(x,y,z),y),
        sympy.diff(uz(x,y,z),z),
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
    l=Library(m,listorder(monoterms,2))


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

    # m,d=lhs_library.vector_to_matrix(lhs)
    # print(d)
    # if d:
    #     raise ValueError("lhs is not in the library: {}",d)

    #lhs=[1]
    #lhs_library=Library(m,[1])

    follower=collect_follower(lhs,
        l,
        Constrainer(trrefl),
        Constrainer(trxrot,parameter_fix=[0,],differential=th),
        Constrainer(tryrot,parameter_fix=[0,],differential=th),
        Constrainer(trzrot,parameter_fix=[0,],differential=th),
        #exptrrot(1),
        lhs_library=lhs_library,
        print_progress=True)
    print("follower")
    for f in follower:
        display(sympy.Matrix(f))

if __name__ == "__main__":
    test_apoblast()

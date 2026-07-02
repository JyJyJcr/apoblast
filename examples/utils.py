import sympy


def listdorder(atom, cords, dorder, dorder_max, filter=(lambda d: True)):
    # print("listdorder",atom,cords,dorder,dorder_max)
    if dorder_max < dorder:
        return []
    elif len(cords) == 0:
        if filter(dorder):
            return [(atom, dorder)]
        else:
            print([])
            return []
    else:
        arr = []
        for i in range(0, dorder_max - dorder + 1):
            arr += map(
                lambda t: (sympy.diff(t[0], (cords[0], i)), t[1]),
                listdorder(atom, cords[1:], dorder + i, dorder_max, filter=filter),
            )
        # print(arr)
        return arr


def listorder(
    monoterms, order, dorder, order_max, dorder_max, filter=(lambda o, d: True)
):
    if order_max < order or dorder_max < dorder:
        return []
    elif len(monoterms) == 0:
        if filter(order, dorder):
            return [sympy.Number(1)]
        else:
            return []
    else:
        arr = []
        for i in range(0, order_max - order + 1):
            arr += map(
                lambda t: (monoterms[0][0] ** i) * t,
                listorder(
                    monoterms[1:],
                    order + i,
                    dorder + i * monoterms[0][1],
                    order_max,
                    dorder_max,
                    filter=filter,
                ),
            )
        return arr

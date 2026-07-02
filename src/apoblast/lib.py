# %env SYMPY_DEBUG=True

# memo
# nontrivial functions: sympy.simplify, sympy.expand sympy.diff as_coefficients_dict subs
# connected_components nullspace extract

# import sys
# import mpmath
# sys.modules['sympy.mpmath'] = mpmath
import sympy
from types import MappingProxyType
import logging

# slience pylance
from sympy.core.symbol import Symbol as Symbol_sp
from sympy.core.function import UndefinedFunction as UndefinedFunction_sp
from sympy.core.function import Function as Function_sp
from sympy.core.function import AppliedUndef as AppliedUndef_sp
from sympy.core.function import Derivative as Derivative_sp
from sympy.core.expr import AtomicExpr as AtomicExpr_sp
from sympy.core.add import Add as Add_sp
from sympy.core.mul import Mul as Mul_sp
from sympy.core.power import Pow as Pow_sp

logger = logging.getLogger(__name__)

# describe some field system
# - continuous index -> coordinate
# - discrete index -> finite dimention vector field


class Model:
    # coord
    # field
    # applied_field
    def __init__(self, coord, field):
        self.coord = tuple(coord)
        self.field = tuple(field)
        if not all(isinstance(e, Symbol_sp) for e in self.coord):
            print(self.coord)
            raise TypeError("coord: into<tuple<sympy.core.symbol.Symbol>>")
        if not all(isinstance(e, UndefinedFunction_sp) for e in self.field):
            raise TypeError("field: into<tuple<sympy.core.function.UndefinedFunction>>")
        # self.coord=coord
        # self.field=field
        self.applied_field = tuple(map(lambda e: e(*coord), field))

    # @property
    # def applied_field(self):
    #     return tuple(self.__applied_field)

    def verify(self, expr, with_field=True, parameter=()):
        if isinstance(expr, AtomicExpr_sp):
            if isinstance(expr, Symbol_sp):
                return expr in self.coord or expr in parameter
            else:
                return True
        elif isinstance(expr, Derivative_sp):
            orig, *derivs = expr.args
            return self.verify(orig, with_field, parameter) and all(
                axis in self.coord for axis, _order in derivs
            )
        elif isinstance(expr, AppliedUndef_sp):
            return with_field and self.coord == expr.args and (expr.func in self.field)
        elif isinstance(
            expr,
            (
                Function_sp,
                Add_sp,
                Mul_sp,
                Pow_sp,
            ),
        ):
            return all(self.verify(e, with_field, parameter) for e in expr.args)
        else:
            return False

    def __repr__(self):
        return f"Model(coord={self.coord}, field={self.field})"


# stock valid terms with a specific model


class Library:
    # model
    # term
    def __init__(self, model, term):
        if not isinstance(model, Model):
            raise TypeError("model: Model")
        self.model = model
        # raise TypeError("term: iter<[valid term]>")
        self.term = tuple(term)
        if not all(self.model.verify(e) for e in self.term):
            raise TypeError("term: into<tuple<<[valid term]>>")

    def classify(self, expr):
        arr = []
        dic = sympy.expand(
            expr,
        ).as_coefficients_dict(*self.model.coord, *self.model.field)
        for i in range(0, len(self.term)):
            fun = self.term[i]
            if fun in dic:
                arr.append(dic.pop(fun))
            else:
                arr.append(0)
        return arr, dic

    # extract expression matrix of some vector field

    def vector_to_cmatrix(self, field):
        mat = []
        dic_arr = []
        for p in field:
            arr, dic = self.classify(p)
            # if dic:
            #     raise RuntimeError("field "+ str(field) +" is not expressible in the library")
            mat.append(arr)
            dic_arr.append(dic)
        return sympy.Matrix(mat), dic_arr

    # extract sympy expression of some vector field

    def cmatrix_to_vector(self, field):
        d, l = sympy.shape(field)
        vec = []
        for i in range(0, d):
            tmp = 0
            for j in range(0, l):
                tmp += field[i, j] * self.term[j]
            vec.append(tmp)
        return tuple(vec)


def cmatrix_to_cdvector(matrix):
    return matrix.transpose().reshape(matrix.cols * matrix.rows)


def cdvector_to_cmatrix(vector, dim):
    return vector.reshape(((vector.rows) // dim), dim).transpose()


# describe transformation procedure:
# - coordinate rules
# - field rewrite rules
# these two are enough to define the trans if the two conditon below is satisfied:
# - trans is local
# - govering eq is local and expressed only with derivatives and various f: R^k->R s


class Transformation:
    # model
    # coord_replace
    # field_replace
    # jacovian_inverse
    # parameters
    def __init__(
        self, model, coord_replace, field_replace, parameter=(), jacovian_inverse=None
    ):
        if not isinstance(model, Model):
            raise TypeError("model: Model")
        # if not isinstance(coord_replace,tuple):
        #     raise TypeError("coord_replace: list<[new coord]>")
        # if not isinstance(field_replace,tuple):
        #     raise TypeError("field_replace: list<[new field]>")
        if not all(
            (isinstance(p, Symbol_sp) and p not in model.coord) for p in parameter
        ):
            raise TypeError("parameter: list<[fresh]>")
        coord_replace = tuple(map(lambda t: sympy.simplify(t), coord_replace))
        field_replace = tuple(map(lambda t: sympy.simplify(t), field_replace))
        if not all(
            model.verify(v, with_field=False, parameter=parameter)
            for v in coord_replace
        ) or len(coord_replace) != len(model.coord):
            # print(all(model.verify(v,with_field=False,parameter=parameter) for v in coord_replace),coord_replace,model.coord)
            raise TypeError("coord_replace: list<[new coord]>")
        if not all(model.verify(v, parameter=parameter) for v in field_replace) or len(
            field_replace
        ) != len(model.field):
            raise TypeError("field_replace: list<[new field]>")

        self.model = model
        self.parameter = tuple(parameter)
        self.coord_replace = coord_replace
        self.field_replace = field_replace
        self.coord_replace_dict = MappingProxyType(
            {model.coord[i]: coord_replace[i] for i in range(0, len(model.coord))}
        )
        self.field_replace_dict = MappingProxyType(
            {model.field[i]: field_replace[i] for i in range(0, len(model.field))}
        )
        m_coord = sympy.Matrix(model.coord)
        m_recoord = sympy.Matrix(coord_replace)
        if jacovian_inverse != None:
            self.jacovian_inverse = jacovian_inverse
        else:
            self.jacovian_inverse = sympy.simplify(m_recoord.jacobian(m_coord).inv())
        self.jacovian_inverse_dict = MappingProxyType(
            {
                model.coord[i]: self.jacovian_inverse.col(i)
                for i in range(0, len(model.coord))
            }
        )

    def apply(self, expr, cache=None):
        # print(len(cache))
        if cache == None:
            return self.apply_internal(expr, cache)
        elif expr in cache:
            # print("cache hit")
            return cache[expr]
        else:
            new = self.apply_internal(expr, cache)
            cache[expr] = new
            return new

    def apply_internal(self, expr, cache):
        # print(expr,expr.func)
        if isinstance(expr, AtomicExpr_sp):
            # atomic element: coord or const or literal const
            if isinstance(expr, Symbol_sp):
                if expr in self.coord_replace_dict:
                    return self.coord_replace_dict[expr]
                else:
                    raise RuntimeError("Unrecognized coord found: " + expr)
            else:
                return expr
        elif isinstance(expr, Derivative_sp):
            # derivative
            orig, *internal_deriv, deriv = expr.args
            axis, order = deriv
            # print("deriv by",axis)
            if order > 1:
                internal = sympy.Derivative(orig, *internal_deriv, (axis, order - 1))
            else:
                if len(internal_deriv) > 0:
                    internal = sympy.Derivative(orig, *internal_deriv)
                else:
                    internal = orig
            transed_internal = self.apply(internal, cache)

            if axis in self.jacovian_inverse_dict:
                return sympy.Add(
                    *(
                        self.jacovian_inverse_dict[axis][i]
                        * sympy.diff(transed_internal, self.model.coord[i])
                        for i in range(0, len(self.model.coord))
                    )
                )
            else:
                raise RuntimeError("Unrecognized derivative found: " + expr)

        elif isinstance(expr, AppliedUndef_sp):
            # field
            # print("field")
            if expr.func in self.field_replace_dict:
                return self.field_replace_dict[expr.func]
            else:
                raise RuntimeError("Unrecognized function symbol found: " + expr)
        elif isinstance(
            expr,
            (
                Function_sp,
                Add_sp,
                Mul_sp,
                Pow_sp,
            ),
        ):
            # other all single point function
            # print("normal func")
            return expr.func(*map(lambda e: self.apply(e, cache), expr.args))
        else:
            raise RuntimeError("Unrecognized func found: " + expr)

    def __repr__(self):
        coord_replace_str = ""
        for c_old, c_new in zip(self.model.coord, self.coord_replace):
            coord_replace_str += str(c_old) + " |-> " + str(c_new) + ", "

        field_replace_str = ""
        for f_old, f_new in zip(self.model.field, self.field_replace):
            field_replace_str += str(f_old) + " |-> " + str(f_new) + ", "

        return f"Transformation(coord_replace=[{coord_replace_str[:-2]}], field_replace=[{field_replace_str[:-2]}], parameter={self.parameter})"


class Constrainer:
    # library
    # transformation
    # parameter_fix
    # differential
    def __init__(self, transformation, parameter_fix=(), differential=None):
        if not len(parameter_fix) == len(transformation.parameter):
            raise TypeError("parameter_fix: into<tuple<[number]>>")
        self.transformation = transformation
        self.parameter_fix = parameter_fix
        self.differential = differential

    def __repr__(self):
        return f"Constrainer(transformation={self.transformation}, parameter_fix={self.parameter_fix}, differential={self.differential})"


# expression of transformation in a specific library
# note: this convertion is parallelizable, but currently not impl


class ConstraintMatrix:
    # library
    # mat
    # exmat
    # exmat_key
    # is_diff
    def __init__(self, library, constrainer):
        if library.model != constrainer.transformation.model:
            raise RuntimeError("library should use same model as constrainer")
        self.library = library

        self.is_diff = constrainer.differential != None

        # cache={}

        # out=joblib.Parallel(n_jobs=1)(joblib.delayed(lambda t: library.classify(constrainer.transformation.apply(t,None)))(t) for t in library.term)
        # print(out)

        def trans(t):
            t = constrainer.transformation.apply(t, None)
            if constrainer.differential != None:
                t = sympy.diff(t, constrainer.differential)
            for i in range(0, len(constrainer.parameter_fix)):
                # print(i,constrainer.transformation.parameter[i],constrainer.parameter_fix[i])
                t = t.subs(
                    constrainer.transformation.parameter[i],
                    constrainer.parameter_fix[i],
                )
            return t

        out = [library.classify(trans(t)) for t in library.term]

        mat = []
        exmat_dic = []
        for arr, dic in out:
            mat.append(arr)
            exmat_dic.append(dic)

        mat = sympy.Matrix(mat)
        self.__mat = mat

        exmat_all_key = list(set().union(*exmat_dic))
        exmat = []
        for i in range(0, len(exmat_dic)):
            arr = []
            for k in exmat_all_key:
                if k in exmat_dic[i]:
                    arr.append(exmat_dic[i][k])
                else:
                    arr.append(0)
            exmat.append(arr)

        if len(exmat_all_key) == 0:
            self.__exmat = None
            self.__exmat_key = None
        else:
            exmat = sympy.Matrix(exmat)
            self.__exmat = exmat
            self.__exmat_key = tuple(exmat_all_key)

    @property
    def mat(self):
        return sympy.Matrix(self.__mat)

    @property
    def exmat(self):
        if self.__exmat == None:
            return None
        else:
            return sympy.Matrix(self.__exmat)

    @property
    def exmat_key(self):
        if self.__exmat_key == None:
            return None
        else:
            return tuple(self.__exmat_key)


# default nullspace solver (very slow with irrational,large dim!!!)


def find_nulls(kernel):
    k_h, k_w = sympy.shape(kernel)
    K = sympy.Matrix.vstack(
        sympy.Matrix.hstack(sympy.Matrix.zeros(k_w), kernel.transpose()),
        sympy.Matrix.hstack(kernel, sympy.Matrix.zeros(k_h)),
    )

    # k_h,k_w=sympy.shape(kernel)
    # if k_w < k_h:
    #     kernel=kernel.row_join(sympy.Matrix.zeros(k_h,k_h-k_w))

    nulls = []

    for conn in K.connected_components():
        col = []
        row = []
        for idx in conn:
            if idx < k_w:
                col.append(idx)
            else:
                row.append(idx - k_w)

        b = kernel.extract(row, col)
        # if len(col) > 0 and len(row) > 0:
        #   display(b)
        # print(b,row,col)
        for nv in b.nullspace():
            vec = sympy.Matrix.zeros(k_w, 1)
            for i in range(len(col)):
                idx = col[i]
                vec[idx, 0] = nv[i, 0]
            nulls.append(vec)
    # totaldim,_d=sympy.shape(B)
    # topdim=0
    # for b in blist:
    #     display(b)
    #     thisdim,_d=sympy.shape(b)
    #     for i in range(topdim,topdim+thisdim):
    #         print(P[i])

    #     bottomdim=totaldim-thisdim-topdim
    #     print("sh",b,topdim,thisdim,totaldim)

    #     nullmat_part=sympy.Matrix.hstack(*nulls_part)

    #     _d,w=sympy.shape(nullmat_part)
    #     if w > 0:
    #         arr=[]
    #         if topdim > 0:
    #             arr.append(sympy.Matrix.zeros(topdim,w))
    #         arr.append(nullmat_part)
    #         if bottomdim > 0:
    #             arr.append(sympy.Matrix.zeros(bottomdim,w))

    #         nullmats.append(sympy.Matrix.vstack(*arr))
    #     topdim+=thisdim

    return sympy.Matrix.hstack(*nulls)
    # mat=sympy.matrices.sparse.SparseMatrix(kernel)
    # nulls=mat.nullspace()


def find_eigen(kernel, filter=None):
    # print("eigen",kernel)
    vs = []
    r = kernel.eigenvects()
    for eigenvalue, multiplicity, vect in r:
        if filter == None or filter(eigenvalue):
            # print("eigenvalue",eigenvalue)
            for v in vect:
                # print("eigenvector",v)
                vs.append(sympy.Matrix(v))
    return sympy.Matrix.hstack(*vs)


class Solver:
    # __kernel_matrix
    # __library
    # __dim
    def __init__(self, library, dim):
        self.__library = library
        self.__dim = dim
        self.__kernel_matrix = sympy.eye(dim * len(self.__library.term))
        self.__dynex = set()

    def condition_null(self, mat):
        nullmat = find_nulls(mat * self.__kernel_matrix)
        self.__kernel_matrix = self.__kernel_matrix * nullmat

    def condition_eigen(self, mat, filter=None):
        eigenmat = find_eigen(
            self.__kernel_matrix.pinv() * mat * self.__kernel_matrix, filter=filter
        )
        self.__kernel_matrix = self.__kernel_matrix * eigenmat

    def shrink(self):
        l, n = sympy.shape(self.__kernel_matrix)
        d = self.__dim
        nonnullrows = []
        for i in range(0, l):
            flag = True
            for j in range(0, d):
                for k in range(0, n):
                    # print(i,d,j,k,i*d+j,k)
                    if self.__kernel_matrix[i * d + j, k] != 0:
                        flag = False
            if not flag:
                nonnullrows.append(i)
            else:
                self.__dynex.add(self.__library.term[i])
                # print("remove",library.term[i])
        lib = Library(
            self.__library.model, [self.__library.term[i] for i in nonnullrows]
        )
        self.__library = lib

        kernnr = []
        for idx in nonnullrows:
            for i in range(0, d):
                kernnr.append(d * idx + i)

        ker = self.__kernel_matrix.extract(kernnr, range(n))
        self.__kernel_matrix = ker

    @property
    def library(self):
        return self.__library


# search vector fields ([follower]) satisfying the below conditions:
# - each component of each [follower] is expressible as the linear combiantion of the terms of the [library]
# - each [trans] works as the same matrix production to both [lhs] and [follower]
# this means the governing eq: [lhs]=sum_i a_i*[follower]_i is covariant with [trans].


def collect_follower_raw(
    lhs,
    library,
    *trans,
    lhs_library=None,
    nullspace_finder=find_nulls,
):
    if lhs_library == None:
        lhs_library = library

    # dimension of LHS
    m = len(lhs)
    a = len(library.term)
    lhsm, dicarr = lhs_library.vector_to_cmatrix(lhs)
    if not all(len(dic) == 0 for dic in dicarr):
        raise RuntimeError("field " + str(lhs) + " is not expressible in the library")

    kernel_matrix = sympy.eye(m * a)
    logger.info("initial:")

    logger.info("  kernel dim: " + str(m * a))
    logger.info("  library dim: " + str(a))

    dynex = set()

    for step in range(0, len(trans)):
        # original expression

        # X = L t L^+
        # L d = 0
        # X R = R T
        # R D = 0

        # debug
        t = trans[step]
        logger.info("stage " + str(step) + ":")
        logger.debug("  constrainer: " + str(t))

        ltm = ConstraintMatrix(lhs_library, t)
        if library == lhs_library:
            tm = ltm
        else:
            tm = ConstraintMatrix(library, t)
        logger.debug("  transformation linearized")

        # L d = 0
        ltme = ltm.exmat
        if ltme != None:
            _R, lexl = sympy.shape(ltme)
            if lhsm * ltme != sympy.zeros(m, lexl):
                print(ltm.exmat_key)
                raise RuntimeError("lhs should not have excluded terms")
        logger.debug("  excluded terms checked")

        # collect all condition
        cond = []

        # lt_pinv=(lhsm*tp.mat).pinv()
        # X=sympy.simplify(lhsm*lt_pinv)
        # if t.differential != None:
        #     X_diff_gen_l=-X*lhsm
        #     X_diff_gen_r=lt_pinv

        lhspinv = lhsm.pinv()

        X = lhsm * ltm.mat * lhspinv
        T = tm.mat
        if (not tm.is_diff) and X.rank() < m:
            raise RuntimeError("transformation should be invertible")
        # RT - XR = 0
        for ia in range(0, a):
            for nu in range(0, m):
                arr = []
                for ib in range(0, a):
                    for mu in range(0, m):
                        tmp = 0
                        if mu == nu:
                            tmp += T[ib, ia]
                        if ia == ib:
                            tmp -= X[nu, mu]
                        arr.append(tmp)
                cond.append(arr)

        D = tm.exmat
        if D != None:
            _m, b = sympy.shape(D)
            for il in range(0, b):
                for nu in range(0, m):
                    arr = []
                    for ib in range(0, a):
                        for mu in range(0, m):
                            if mu == nu:
                                arr.append(D[ib, il])
                            else:
                                arr.append(0)
                    cond.append(arr)

            logger.debug("  excluded terms:")
            for f in tm.exmat_key:
                if f not in dynex:
                    logger.info(f)

        cond = sympy.Matrix(cond)

        logger.debug("  condition collected")
        kermat = sympy.simplify(cond * kernel_matrix)
        logger.debug("  condition simplified")

        nullmat = nullspace_finder(kermat)
        kernel_matrix = kernel_matrix * nullmat

        _l, n = sympy.shape(kernel_matrix)

        # shrink library

        # print("old")
        # display(kernel_matrix)
        # display(library.term)

        l = len(library.term)
        nonnullrows = []
        for i in range(0, l):
            flag = True
            for j in range(0, m):
                for k in range(0, n):
                    # print(i,d,j,k,i*d+j,k)
                    if kernel_matrix[i * m + j, k] != 0:
                        flag = False
            if not flag:
                nonnullrows.append(i)
            else:
                dynex.add(library.term[i])
                # print("remove",library.term[i])
        _lib = Library(library.model, [library.term[i] for i in nonnullrows])
        library = _lib
        a = len(library.term)

        kernnr = []
        for idx in nonnullrows:
            for i in range(0, m):
                kernnr.append(m * idx + i)

        kernel_matrix = kernel_matrix.extract(kernnr, range(n))

        logger.debug("  library shrinked")

        logger.info("  kernel dim: " + str(n))
        logger.info("  library dim: " + str(a))

        # print("new")
        # display(kernel_matrix)
        # display(library.term)

    # display(kernel_matrix)
    rk, _p = kernel_matrix.transpose().rref()
    kernel_matrix = rk.transpose()
    # display(kernel_matrix)

    return kernel_matrix, library


def collect_follower(
    lhs,
    library,
    *trans,
    lhs_library=None,
    nullspace_finder=find_nulls,
):
    d = len(lhs)
    kernel_matrix, library = collect_follower_raw(
        lhs,
        library,
        *trans,
        lhs_library=lhs_library,
        nullspace_finder=nullspace_finder,
    )
    R = len(library.term)
    li = []
    _ld, n = sympy.shape(kernel_matrix)
    for i in range(0, n):
        v = kernel_matrix.col(i)
        m = cdvector_to_cmatrix(v, d)
        li.append(library.cmatrix_to_vector(m))
    return tuple(li)

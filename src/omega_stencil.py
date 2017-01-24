"""Omega stencil

Creates a C-language subroutine that calculates the stencil weights
corresponding to a simple finite difference approximation of the LHS
of the generalized omega equation. The resulting subroutine
can be used in Petsc, for example.
"""

from sympy import symbols, IndexedBase, Idx, apply_finite_diff, expand
from sympy.printing.ccode import CCodePrinter


def def_fd(xyz, ijk, hxhyhz, nns):
    """A wrapper function for sympy's apply_finite_diff()

    Returns a function that is suitable for generating finite difference
    stencils.

    Arguments:

    xyz     --  the list of free variables
    ijk     --  the dictionary of indices corresponding to the free variables
    hxhyhz  --  the dictionary of grid spacings in directions
                corresponding to the free variables
    nns     --  the list of points to use for the stencil

    Example:

    >>> f = IndexedBase('f')
    >>> x = Symbol('x')
    >>> i = Idx('i')
    >>> hx = Symbol('hx')
    >>> D = def_fd([x],{x:i},{x:hx},[-1,0,1])
    >>> D(f[i],x)
    -f[i - 1]/(2*hx) + f[i + 1]/(2*hx)
    >>> D(f[i],x,x)
    f[i - 1]/hx**2 + f[i + 1]/hx**2 - 2*f[i]/hx**2
    """
    def _fd(func, *ds):
        if len(ds) == 0:
            return func
        else:
            assert ds[0] in xyz
            assert all([ix in ijk.values() for ix in func.atoms(Idx)])
            order = ds.count(ds[0])
            idx = ijk[ds[0]]
            rrs, ffs = zip(*[(m * hxhyhz[ds[0]], func.replace(idx, idx + m))
                             for m in nns])
            res = apply_finite_diff(order, rrs, ffs, 0)
            assert res != 0
            return _fd(res, *[d for d in ds if d != ds[0]])
    return _fd

def get_stencil(lhs, wrt, nns):
    """get_stencil returns a list of stencil elements

    Each stencil element is a pair of stencil k,j,i - coordinates
    and the weight coefficient.

    Arguments:

    lhs     --  the expression from which to extract the stencil
    wrt     --  the name of the array which coefficients are extracted
    nns     --  the list of points to use for the stencil
    """
    lhs = expand(lhs)
    stencil_elements = list(
        ((p, n, m), lhs.coeff(wrt[k + p, j + n, i + m]))
        for m in nns for n in nns for p in nns
        if lhs.coeff(wrt[k + p, j + n, i + m]) != 0)
    newlhs = sum(list(c * wrt[[x + y for x, y in zip((k, j, i), m)]]
                      for m, c in stencil_elements))
    assert expand(lhs - newlhs) == 0
    return stencil_elements


# Symbols in the discretized equations

f, V, sigma, zeta, omega = symbols('f V sigma zeta omega', cls=IndexedBase)
x, y, z = symbols('x y z')
i, j, k = symbols('i j k', cls=Idx)
hx, hy, hz = symbols('hx hy hz')

D = def_fd([x, y, z], {x: i, y: j, z: k}, {x: hx, y: hy, z: hz}, [-1, 0, 1])


# The LHS terms in the generalized omega equation (slightly re-organized)

L1 = (D(sigma[k, j, i] * omega[k, j, i], x, x) +
      D(sigma[k, j, i] * omega[k, j, i], y, y))

L2 = (f[j] * D(zeta[k, j, i], z) * D(omega[k, j, i], z) +
      f[j] * (zeta[k, j, i] + f[j]) * D(omega[k, j, i], z, z))

L3 = (- f[j] * D(zeta[k, j, i], z, z) * omega[k, j, i]
      - f[j] * D(zeta[k, j, i], z) * D(omega[k, j, i], z))

L4 = (f[j] * D(V[k, j, i, 0], z, z) * D(omega[k, j, i], y)
      - f[j] * D(V[k, j, i, 1], z, z) * D(omega[k, j, i], x)
      + f[j] * D(V[k, j, i, 0], z) * D(omega[k, j, i], y, z)
      - f[j] * D(V[k, j, i, 1], z) * D(omega[k, j, i], x, z))

L = (L1 + L2 + L3 + L4) * hx * hy * hz

omega_stencil = get_stencil(L, omega, [-1,0,1])


# Print as a C-function

class MyCCodePrinter(CCodePrinter):
    """MyCCodePrinter

    Minor modifications to CCodePrinter.
    """
    def _print_Indexed(self, expr):
        return "%s[%s]" % (self._print(expr.base.label), ']['.join(
            [self._print(ind).replace(' ', '') for ind in expr.indices]))

    def _print_Pow(self, expr):
        if ("Pow" not in self.known_functions) and (expr.exp == 2):
            return '%s*%s' % tuple([self._print(expr.base)] * 2)
        else:
            return super(MyCCodePrinter, self)._print_Pow(expr)


def myccode(expr, assign_to=None, **settings):
    """Like ccode, but with MyCCodePrinter"""
    return MyCCodePrinter(settings).doprint(expr, assign_to)

def print_function(stencil):
    """Writes C source for Petsc
    """
    print """static inline void omega_stencil(int i,int j,int k,
    PetscScalar hx,PetscScalar hy,PetscScalar hz,PetscScalar *f,
    const PetscScalar ***sigma,const PetscScalar ***zeta,
    const PetscScalar ****V,int *n,
    MatStencil *col,PetscScalar *w)\n{"""
    print "*n = %d;" % len(stencil)
    for idx, val in enumerate(stencil):
        print ' '.join(["col[{:2}].{} = {:5} ;".format(idx, n, n + m)
                        for n, m in zip((k, j, i), val[0])])
    for idx, val in enumerate(stencil):
        print myccode(val[1], assign_to='w[{:2}]'.format(idx),
                      contract=False)
    print "}"

print_function(omega_stencil)

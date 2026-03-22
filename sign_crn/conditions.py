r"""
Conditions for CBE
==================

We consider conditions for existence and uniqueness of positive complex-balanced equilibria (CBE)
of (chemical) reaction networks.
These conditions also describe injectivity and bijectivity
of polynomial and exponential maps as seen in [MHR19]_.

In both cases, we apply the conditions to a pair of matrices.

Robustness of existence and uniqueness
--------------------------------------

We consider the following matrices::

    sage: from sign_crn.conditions import *
    sage: P = matrix([[1, 0, 1, 0], [0, 0, 0, 1]])
    sage: P
    [1 0 1 0]
    [0 0 0 1]
    sage: Pt = matrix([[1, 0, 1, 1], [0, 1, 0, -1]])
    sage: Pt
    [ 1  0  1  1]
    [ 0  1  0 -1]

To study robustness,
we consider a condition involving maximal minors::

    sage: closure_condition(P, Pt)
    True

There is also an equivalent condition using sign vectors::

    sage: from sign_crn.conditions import closure_condition_sign_vectors
    sage: closure_condition_sign_vectors(P, Pt)
    True

Now, we consider an example involving parameters::

    sage: var("a, b, c")
    (a, b, c)
    sage: P = matrix([[c, 1, c]])
    sage: P
    [c 1 c]
    sage: Pt = matrix([[a, b, -1]])
    sage: Pt
    [ a  b -1]

We obtain the following condition on the variables::

    sage: closure_condition(P, Pt) # random
    [{-b > 0, c == 0},
     {-b < 0, c == 0},
     {-b > 0, c > 0, -a*c > 0},
     {-b < 0, c < 0, -a*c < 0}]

Thus, there are four possibilities to set the variables:
From the first two sets of conditions, we see that the closure condition is satisfied
if :math:`c` is zero and :math:`b` is nonzero.
The closure condition is also satisfied if :math:`a` and :math:`b` are negative and :math:`c` is positive
or if :math:`a` and :math:`b` are positive and :math:`c` is negative.

Uniqueness
----------

We define the following matrices::

    sage: P = matrix([[1, 1, 1]])
    sage: Pt = matrix([[1, 0, 1]])

The package uses maximal minors to study uniqueness::

    sage: uniqueness_condition(P, Pt)
    True

Instead, we can consider the oriented matroids determined by these matrices::

    sage: from sign_crn.conditions import uniqueness_condition_sign_vectors
    sage: uniqueness_condition_sign_vectors(P, Pt)
    True

Now, we consider another example::

    sage: P = matrix([[1, 1, 1]])
    sage: Pt = matrix([[1, -1, 1]])

The condition is violated::

    sage: uniqueness_condition(P, Pt)
    False

Finally, we consider an example with parameters :math:`a, b \in \mathbb{R}`::

    sage: var("a, b")
    (a, b)
    sage: P = matrix([[1, 1, 1]])
    sage: P
    [1 1 1]
    sage: Pt = matrix([[a, b, -1]])
    sage: Pt
    [ a  b -1]

Here, the function returns a system of inequalities::

    sage: uniqueness_condition(P, Pt) # random order
    [{-a >= 0, -b >= 0}]

Hence, the condition holds if and only if :math:`a, b \leq 0`.

Note that we cannot apply :func:`~uniqueness_condition_sign_vectors` because of the parameters.

Existence and uniqueness
------------------------

We consider the following matrices::

    sage: P = matrix([[1, 0, 1, 0], [0, 1, 0, 1]])
    sage: P
    [1 0 1 0]
    [0 1 0 1]
    sage: Pt = matrix([[1, 0, 0, -1], [0, 1, 1, 1]])
    sage: Pt
    [ 1  0  0 -1]
    [ 0  1  1  1]

The uniqueness condition is satisfied::

    sage: uniqueness_condition(P, Pt)
    True

We check the face condition::

    sage: face_condition(P, Pt)
    True

Finally, the nondegeneracy condition delivers::

    sage: nondegeneracy_condition(P, Pt)
    True

Let us consider another example by swapping the matrices.
The uniqueness condition is symmetric in the two matrices and thus holds again::

    sage: uniqueness_condition(Pt, P)
    True

However, the face condition is violated::

    sage: face_condition(Pt, P)
    False

Now, we consider a map involving a parameter.
(see Example 20 of [MHR19]_)::

    sage: var("a")
    a
    sage: assume(a > 0)
    sage: P = matrix([[0, 0, 1, 1, -1, 0], [1, -1, 0, 0, 0, -1], [0, 0, 1, -1, 0, 0]])
    sage: P
    [ 0  0  1  1 -1  0]
    [ 1 -1  0  0  0 -1]
    [ 0  0  1 -1  0  0]
    sage: Pt = matrix([[1, 1, 0, 0, -1, a], [1, -1, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0]])
    sage: Pt
    [ 1  1  0  0 -1  a]
    [ 1 -1  0  0  0  0]
    [ 0  0  1 -1  0  0]

The first two conditions depend on the sign vectors corresponding
to the rows of these matrices which are independent of the specific value for :math:`a`::

    sage: uniqueness_condition(P, Pt)
    True
    sage: face_condition(P, Pt)
    True

The nondegeneracy condition is satisfied
for :math:`a \in (0, 1) \cup (1, 2)`.
We demonstrate this for some values::

    sage: nondegeneracy_condition(P, Pt(a=1/2))
    True
    sage: nondegeneracy_condition(P, Pt(a=3/2))
    True

On the other hand, this condition does not hold if
:math:`a \in \{1\} \cup [2, \infty)`::

    sage: nondegeneracy_condition(P, Pt(a=1))
    False
    sage: nondegeneracy_condition(P, Pt(a=2))
    False
    sage: nondegeneracy_condition(P, Pt(a=3))
    False

To certify the result, we call::

    sage: nondegeneracy_condition(P, Pt(a=1), certify=True)
    (False, (1, 1, 0, 0, -1, 1))

Hence, the positive support of the vector ``v = (1, 1, 0, 0, -1, 1)`` in ``row(Pt)``
can be covered by a sign vector ``(++000+)`` corresponding to ``ker(P)``.
"""

#############################################################################
#  Copyright (C) 2026                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from copy import copy

from sage.combinat.combination import Combinations
from sage.matrix.constructor import Matrix
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject

from sign_vectors import SignVector, OrientedMatroid
from elementary_vectors import circuit_kernel_matrix
from elementary_vectors.utility import is_constant
from certlin import Intervals, LinearInequalitySystem
from .utility import (
    closure_minors_utility,
    intervals_to_sign_vectors,
    sign_vector_to_intervals,
    non_negative_circuits_from_matrix,
    non_negative_cocircuits_from_matrix,
    equal_entries_lists,
    vector_from_sign_vector
)


def closure_condition(matrix1: Matrix, matrix2: Matrix) -> bool | list[set]:
    r"""
    Closure condition using maximal minors.

    OUTPUT:
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in (at least) one of these sets are satisfied.

    .. SEEALSO::

        :func:`closure_condition_sign_vectors`

    .. NOTE::

        The matrices need to have maximal rank and the same dimensions.
        Otherwise, a ``ValueError`` is raised.
    """
    positive_found = False
    negative_found = False
    symbolic_pairs = set()
    for indices in Combinations(matrix1.ncols(), matrix1.nrows()):
        minor1 = matrix1.matrix_from_columns(indices).det()
        if not minor1:
            continue
        minor2 = matrix2.matrix_from_columns(indices).det()
        if not minor2:
            return False
        product = minor1 * minor2
        if not is_constant(product):
            symbolic_pairs.add((minor1, product))
            continue
        if product > 0:
            positive_found = True
        elif product < 0:
            negative_found = True
        if positive_found and negative_found:
            return False

    return closure_minors_utility(symbolic_pairs, positive_found, negative_found)


def closure_condition_sign_vectors(matrix1: Matrix, matrix2: Matrix) -> bool:
    r"""
    Closure condition using sign vectors.

    This condition is about
    whether a set of sign vectors is contained in the closure of another set of sign vectors,
    or equivalently, whether covectors are covered by topes.

    .. SEEALSO::

        :func:`closure_condition`

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~closure_condition`.
    """
    topes = OrientedMatroid(matrix2).topes()
    for covector in OrientedMatroid(matrix1).topes():
        if not any(covector <= tope for tope in topes):
            return False
    return True


def uniqueness_condition(matrix1: Matrix, matrix2: Matrix) -> bool | list[set]:
    r"""
    Uniqueness condition using maximal minors.

    OUTPUT:
    Return whether there exists at most one equilibrium.
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in exactly one of these sets are satisfied.

    .. SEEALSO::

        :func:`uniqueness_condition_sign_vectors`

    .. NOTE::

        The matrices need to have maximal rank and the same dimensions.
        Otherwise, a ``ValueError`` is raised.

    TESTS::

        sage: from sign_crn import *
        sage: var("a, b")
        (a, b)
        sage: P = matrix([[1, 0, -1], [0, 1, -1]])
        sage: P
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Pt = matrix([[1, 0, a], [0, 1, b]])
        sage: Pt
        [1 0 a]
        [0 1 b]
        sage: uniqueness_condition(P, Pt) # random order
        [{-a >= 0, -b >= 0}]
        sage: conditions = uniqueness_condition(P, Pt)[0]
        sage: conditions # random order
        sage: (-a >= 0) in conditions and (-b >= 0) in conditions
        True
        sage: P = matrix([[a, 0, 1, 0], [0, 1, -1, 0], [0, 0, 0, 1]])
        sage: Pt = matrix([[1, 0, 0, -1], [0, b, 1, 1], [0, 0, a, 1]])
        sage: uniqueness_condition(P, Pt) # random order
        [{(a - 1)*a >= 0, a*b >= 0}, {(a - 1)*a <= 0, a*b <= 0}]
        sage: len(_), len(_[0])
        (2, 2)
    """
    positive_product_found = False
    negative_product_found = False
    symbolic_expressions = set()

    for indices in Combinations(matrix1.ncols(), matrix1.nrows()):
        minor1 = matrix1.matrix_from_columns(indices).det()
        if not minor1:
            continue
        product = (
            minor1 * matrix2.matrix_from_columns(indices).det()
        )
        if not is_constant(product):
            symbolic_expressions.add(product)
        elif product > 0:
            positive_product_found = True
        elif product < 0:
            negative_product_found = True
        if positive_product_found and negative_product_found:
            return False
    if positive_product_found:
        if symbolic_expressions:
            return [set(expression >= 0 for expression in symbolic_expressions)]
        return True
    if negative_product_found:
        if symbolic_expressions:
            return [set(expression <= 0 for expression in symbolic_expressions)]
        return True
    if symbolic_expressions:
        return [
            set(expression >= 0 for expression in symbolic_expressions),
            set(expression <= 0 for expression in symbolic_expressions),
        ]
    return False


def uniqueness_condition_sign_vectors(matrix1: Matrix, matrix2: Matrix) -> bool:
    r"""
    Uniqueness condition using sign vectors.

    .. SEEALSO::

        :func:`uniqueness_condition`

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~uniqueness_condition`.

    EXAMPLES::

        sage: from sign_crn.conditions import uniqueness_condition_sign_vectors
        sage: P = matrix([[1, 1, 1]])
        sage: P
        [1 1 1]
        sage: Pt = matrix([[1, 0, 1]])
        sage: Pt
        [1 0 1]
        sage: uniqueness_condition_sign_vectors(P, Pt)
        True
        sage: P = matrix([[1, 0, -1], [0, 1, -1]])
        sage: P
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Pt = matrix([[1, 0, -1], [0, 1, 1]])
        sage: Pt
        [ 1  0 -1]
        [ 0  1  1]
        sage: uniqueness_condition_sign_vectors(P, Pt)
        False

    TESTS::

        sage: from sign_crn.conditions import uniqueness_condition_sign_vectors
        sage: A = identity_matrix(3)
        sage: uniqueness_condition_sign_vectors(A, A)
        True
    """
    covectors = OrientedMatroid(matrix1).covectors()
    counter = 0
    for covector in OrientedMatroid(matrix2).vectors():
        if covector in covectors:
            counter += 1
            if counter > 1:
                return False
    return True


def face_condition(kernel_matrix1: Matrix, kernel_matrix2: Matrix) -> bool:
    r"""
    Face condition using nonnegative cocircuits.

    This condition is about nonnegative cocircuits covering other nonnegative cocircuits.
    """
    non_negative_cocircuits = non_negative_cocircuits_from_matrix(kernel_matrix1)

    for cocircuit1 in non_negative_cocircuits_from_matrix(kernel_matrix2):
        if not any(cocircuit2 <= cocircuit1 for cocircuit2 in non_negative_cocircuits):
            return False
    return True


def nondegeneracy_condition(kernel_matrix1: Matrix, kernel_matrix2: Matrix, certify: bool = False) -> bool:
    r"""
    Nondegeneracy condition.

    This condition is about whether all positive equal components of a vector
    can be covered by nonnegative cocircuits.

    If ``certify`` is true, a list is returned to certify the result.
    (see the examples)

    EXAMPLES:

    We consider the following matrices::

        sage: from sign_crn import *
        sage: P = matrix([[1, 0, 1, 0], [0, 1, 0, 1]])
        sage: Pt = matrix([[1, 0, 0, -1], [0, 1, 0, -1]])
        sage: nondegeneracy_condition(P, Pt)
        True

    Next, we certify the result.
    The corresponding subspaces are trivially nondegenerate
    since there are no nonnegative covectors in the kernel of ``P``::

        sage: nondegeneracy_condition(P, Pt, certify=True)
        (True, 'no nonnegative covectors')

    Now, we consider an example of degenerate subspaces::

        sage: P = matrix([[0, 0, 1]])
        sage: Pt = matrix([[1, 1, 0]])
        sage: nondegeneracy_condition(P, Pt, certify=True)
        (False, (1, 1, 0))

    The resulting vector lies in the row space of ``Pt``.
    The nonnegative covector ``(++0)`` in the kernel of ``P`` covers the first two equal components.

    In the next example, there exists a partial cover::

        sage: P = matrix([[1, 1, 0, 0], [0, 0, 1, -1]])
        sage: Pt = matrix([[1, 1, 0, -1], [0, 0, 1, 0]])
        sage: nondegeneracy_condition(P, Pt, certify=True)
        (True, ([], [[[2, 3]]], [[[[2, 3]], [(--++)]]]))

    In fact, a vector in the row space of ``Pt`` with equal positive components on ``[2, 3]``
    corresponding to ``(--++)`` can be fully covered by covectors.
    However, this vector would not satisfy the condition on the support.
    """
    non_negative_cocircuits = non_negative_circuits_from_matrix(kernel_matrix1)

    if not non_negative_cocircuits:
        if certify:
            return True, "no nonnegative covectors"
        return True

    non_negative_cocircuits = sorted(non_negative_cocircuits, key=lambda covector: len(covector.support()))
    length = kernel_matrix2.ncols()
    degenerate = False

    lower_bounds = [-Infinity] * length
    upper_bounds = [0] * length
    upper_bounds_inf = [Infinity] * length

    kernel_matrix = circuit_kernel_matrix(kernel_matrix2)
    covectors_support_condition = non_negative_cocircuits_from_matrix(kernel_matrix1)

    if certify:
        certificate = []
        certificates_zero_equal_components = []
        certificates_partial_cover = []
        certificate_support_condition = []

    def recursive_degenerate(
        non_negative_cocircuits: set[SignVector],
        matrix_old: Matrix,
        positive_components: list[int],
        lower_bounds: list[int],
        upper_bounds: list[int]
    ):
        r"""
        Recursive function.

        INPUT:

        - ``non_negative_cocircuits`` -- a list of nonnegative sign vectors
        - ``lower_bounds`` -- a list of values ``-Infinity`` and ``1``
        - ``upper_bounds`` -- a list of values ``0`` and ``Infinity``
        """
        nonlocal degenerate
        nonlocal certificate

        while non_negative_cocircuits:
            cocircuit = non_negative_cocircuits.pop()
            lower_bounds_new = copy(lower_bounds)
            upper_bounds_new = copy(upper_bounds)
            for i in cocircuit.support():
                lower_bounds_new[i] = 1
                upper_bounds_new[i] = Infinity

            intervals = Intervals.from_bounds(lower_bounds_new, upper_bounds_new)
            indices_new = positive_components + [cocircuit.support()]
            matrix_new = Matrix(
                matrix_old.rows() + equal_entries_lists(length, cocircuit.support())
            ).echelon_form()
            # TODO don't use kernel matrix? consider evs in row space
            system = LinearInequalitySystem(matrix_new.right_kernel_matrix().T, intervals)

            if system.certify()[0]:
                if certify:
                    covectors_certificate_support_condition = []
                for sv in intervals_to_sign_vectors(intervals):
                    if not system.with_intervals(sign_vector_to_intervals(sv)).certify()[0]:
                        continue
                    if not any(
                        set(cocircuit.support()).issubset(sv.support())
                        for cocircuit in covectors_support_condition
                    ):
                        degenerate = True
                        if certify:
                            certificate = vector_from_sign_vector(
                                system._evs_generator(kernel=False),
                                sv
                            )
                        return
                    if certify:
                        covectors_certificate_support_condition.append(sv)
                if certify:
                    certificate_support_condition.append(
                        [indices_new, covectors_certificate_support_condition]
                    )

            if system.with_intervals(Intervals.from_bounds(lower_bounds_new, upper_bounds_inf)).certify()[0]:
                if certify:
                    certificates_partial_cover.append(indices_new)
                recursive_degenerate(
                    copy(non_negative_cocircuits),
                    matrix_new,
                    indices_new,
                    lower_bounds_new,
                    upper_bounds_new,
                )
            elif certify:
                certificates_zero_equal_components.append(indices_new)

            if degenerate:
                return
        return

    recursive_degenerate(
        non_negative_cocircuits, kernel_matrix, [], lower_bounds, upper_bounds
    )

    if certify:
        if degenerate:
            return not degenerate, certificate
        return not degenerate, (
            certificates_zero_equal_components,
            certificates_partial_cover,
            certificate_support_condition,
        )
    return not degenerate


class ConditionsCRN(SageObject):
    r"""
    Class for conditions for chemical reaction networks.
    """
    def __init__(self):
        self._reduced_stoichiometric_matrix = None
        self._reduced_kinetic_order_matrix = None
        self._stoichiometric_kernel_matrix = None
        self._kinetic_order_kernel_matrix = None

    @property
    def reduced_stoichiometric_matrix(self) -> Matrix:
        r"""
        Return the reduced stoichiometric matrix.
        """
        if self._reduced_stoichiometric_matrix is None:
            if self._stoichiometric_kernel_matrix is None:
                raise ValueError("Neither kernel matrices nor reduced matrices are defined.")
            self._reduced_stoichiometric_matrix = circuit_kernel_matrix(self._stoichiometric_kernel_matrix)
        return self._reduced_stoichiometric_matrix

    @property
    def reduced_kinetic_order_matrix(self) -> Matrix:
        r"""
        Return the reduced kinetic order matrix.
        """
        if self._reduced_kinetic_order_matrix is None:
            if self._kinetic_order_kernel_matrix is None:
                raise ValueError("Neither kernel matrices nor reduced matrices are defined.")
            self._reduced_kinetic_order_matrix = circuit_kernel_matrix(self._kinetic_order_kernel_matrix)
        return self._reduced_kinetic_order_matrix

    @property
    def stoichiometric_kernel_matrix(self) -> Matrix:
        r"""
        Return the kernel matrix of the stoichiometric matrix.
        """
        if self._stoichiometric_kernel_matrix is None:
            if self._reduced_stoichiometric_matrix is None:
                raise ValueError("Neither kernel matrices nor reduced matrices are defined.")
            self._stoichiometric_kernel_matrix = circuit_kernel_matrix(self._reduced_stoichiometric_matrix)
        return self._stoichiometric_kernel_matrix

    @property
    def kinetic_order_kernel_matrix(self) -> Matrix:
        r"""
        Return the kernel matrix of the kinetic order matrix.
        """
        if self._kinetic_order_kernel_matrix is None:
            if self._reduced_kinetic_order_matrix is None:
                raise ValueError("Neither kernel matrices nor reduced matrices are defined.")
            self._kinetic_order_kernel_matrix = circuit_kernel_matrix(self._reduced_kinetic_order_matrix)
        return self._kinetic_order_kernel_matrix

    def are_kernel_matrices_defined(self) -> bool:
        r"""
        Return whether the kernel matrices are defined.
        """
        return self._stoichiometric_kernel_matrix is not None and self._kinetic_order_kernel_matrix is not None

    def are_reduced_matrices_defined(self) -> bool:
        r"""
        Return whether the reduced matrices are defined.
        """
        return self._reduced_stoichiometric_matrix is not None and self._reduced_kinetic_order_matrix is not None

    def closure_condition(self) -> bool | list[set]:
        r"""
        Closure condition using maximal minors.

        OUTPUT:
        If the result depends on variables, a list of sets is returned.
        The condition holds if the inequalities in (at least) one of these sets are satisfied.
        """
        if self.are_kernel_matrices_defined():
            matrix1 = self.stoichiometric_kernel_matrix
            matrix2 = self.kinetic_order_kernel_matrix
        elif self.are_reduced_matrices_defined():
            matrix1 = self.reduced_stoichiometric_matrix
            matrix2 = self.reduced_kinetic_order_matrix
        else:
            raise ValueError("Neither kernel matrices nor reduced matrices are defined.")

        return closure_condition(matrix1, matrix2)

    def uniqueness_condition(self) -> bool | list[set]:
        r"""
        Uniqueness condition using maximal minors.

        OUTPUT:
        If the result depends on variables, a list of sets is returned.
        The condition holds if the inequalities in exactly one of these sets are satisfied.
        """
        if self.are_kernel_matrices_defined():
            matrix1 = self.stoichiometric_kernel_matrix
            matrix2 = self.kinetic_order_kernel_matrix
        elif self.are_reduced_matrices_defined():
            matrix1 = self.reduced_stoichiometric_matrix
            matrix2 = self.reduced_kinetic_order_matrix
        else:
            raise ValueError("Neither kernel matrices nor reduced matrices are defined.")

        return uniqueness_condition(matrix1, matrix2)

    def face_condition(self) -> bool:
        r"""
        Face condition using nonnegative cocircuits.

        This condition is about nonnegative cocircuits covering other nonnegative cocircuits.
        """
        return face_condition(self.stoichiometric_kernel_matrix, self.kinetic_order_kernel_matrix)

    def nondegeneracy_condition(self, certify: bool = False) -> bool:
        r"""
        Nondegeneracy condition.

        This condition is about whether all positive equal components of a vector
        can be covered by nonnegative cocircuits.

        If ``certify`` is true, a list is returned to certify the result.
        (see the examples)
        """
        return nondegeneracy_condition(self.stoichiometric_kernel_matrix, self.kinetic_order_kernel_matrix, certify=certify)

    @staticmethod
    def from_kernel_matrices(stoichiometric_kernel_matrix: Matrix, kinetic_order_kernel_matrix: Matrix) -> "ConditionsCRN":
        r"""
        Construct the object from kernel matrices.

        The kernel matrices are the right kernel matrices of the stoichiometric matrix and the kinetic order matrix.
        """
        obj = ConditionsCRN()
        obj._stoichiometric_kernel_matrix = stoichiometric_kernel_matrix
        obj._kinetic_order_kernel_matrix = kinetic_order_kernel_matrix
        return obj

    @staticmethod
    def from_reduced_matrices(reduced_stoichiometric_matrix: Matrix, reduced_kinetic_order_matrix: Matrix) -> "ConditionsCRN":
        r"""
        Construct the object from the stoichiometric matrix and the kinetic order matrix.
        """
        obj = ConditionsCRN()
        obj._reduced_stoichiometric_matrix = reduced_stoichiometric_matrix.matrix_from_rows(reduced_stoichiometric_matrix.pivot_rows())
        obj._reduced_kinetic_order_matrix = reduced_kinetic_order_matrix.matrix_from_rows(reduced_kinetic_order_matrix.pivot_rows())
        return obj



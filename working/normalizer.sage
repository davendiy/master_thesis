#!/usr/bin/env sage

from itertools import permutations
from collections import defaultdict
import os

a = b = c = d = ...


# noinspection PyShadowingNames
def normalize_expressions(exps, allowed=None):
    """Replace all variables except allowed (default: a b c d) with {x0,x1,x2,x3...}"""
    extra_args = set()

    if allowed is None:
        allowed = var('a b c d')

    for exp in exps:
        extra_args |= {el for el in exp.args() if el not in allowed}

    res = []
    for exp in exps:
        for i, el in enumerate(extra_args):
            exp = exp.subs({var(el): var(f'x{i}')})
        res.append(exp)
    return tuple(res)


def prepare_gap_env():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(dir_path, 'el2word.g')
    gap(f'LoadPackage("Cryst");;Read("{file_path}")')


# noinspection PyShadowingNames
def el2word(G, el, gens=None, verbose=False):
    """Represent, given an element of a finite group G,
    as the word over generators."""

    prepare_gap_env()
    gap.set('G', G)
    gap.set('g', el)

    if gens is None:
        gens = G.GeneratorsOfGroup()

    gens_dict = {f'g_{i}': el for i, el in enumerate(list(gens))}

    gens = '"' + '", "'.join(sorted(gens_dict.keys())) + '"'
    expression = f'FreeGroup({gens})'
    if verbose:
        print("Free group to be created:", expression)
        print("Translation of generators:", gens_dict)

    gap.set("F", expression)

    res = gap('ElementAsWordGeneratorsPointGroup(G, g, F)')
    return str(res).split("*"), gens_dict


# noinspection PyShadowingNames
def extend_permutation(G, el, gens, maps_to):
    """Extend given permutation on generators onto the entire group."""
    word, trans = el2word(G, el, gens=gens)

    if el == G.One():
        return el

    res = G.One()
    for g in word:
        power = 1
        if "^" in g:
            g, power = g.split('^')
            power = int(power)
        index = gens.index(trans[g])
        g_to = maps_to[index]
        #         index_to = perm[index]
        #         g_to = gens[index_to] ** power
        res *= g_to ** power
    return res


alphabet = 'abcdefghijklmnopqrstuvwxyz'


# noinspection PyShadowingNames
def create_symbolic_matrix(dim, use_alphabet=False):
    """Create symbolic matrix of given dimention.

    :param dim: integer
    :param use_alphabet: if True, then symbolic matrices are created like
                                [a, b]
                                [c, d]
                         if False, then
                                [a_00, a_01]
                                [a_10, a_11]
    """
    if use_alphabet and dim * dim > len(alphabet):
        raise ValueError(f"Can't use alphabet for matrix {dim}x{dim} due to lack of letters.")
    A = []
    args = []
    for i in range(dim):
        row = []
        for j in range(dim):
            if use_alphabet:
                row.append(var(alphabet[i * dim + j]))
            else:
                row.append(var(f'a_{i}{j}'))
            args.append(row[-1])
        A.append(row)
    A = matrix(A)
    return A, args


def sol2matrix(solution, dim=2, use_alphabet=False):
    """Transform solution into symbolic matrix."""
    if use_alphabet and dim * dim > len(alphabet):
        raise ValueError(f"Can't use alphabet for matrix {dim}x{dim} due to lack of letters.")
    A = []
    sol_dict = {sol.left(): sol.right() for sol in solution}

    for i in range(dim):
        row = []
        for j in range(dim):
            if use_alphabet:
                row.append(sol_dict[var(alphabet[i * dim + j])])
            else:
                row.append(sol_dict[var(f'a_{i}{j}')])

        A.append(row)
    A = matrix(A)
    return A


# noinspection PyShadowingNames
def normalizers(n, dim=2, use_alphabet=False, verbose=False, normalize_exp=True,
                to_matrix=True, ignore_trivial=True):
    """Find normalizer of the PointGroup in GL(n, QQ).

    Tries to find normalizer as a solution of
    system of linear equations N*A_i = A_j*N, where A_i and A_j are
    two elements of the PointGroup with respect to a permutation. Since
    Normalizer is a group N, that satisfies condition NA = AN,
    we can check all the permutations of A and solve the respective system
    of linear equations.

    Baseline version, slow.
    ----------------------------------------------------------------------------
    :param n: index of the Crystallographic group in the Gap CrystCat package
    :param dim: dimension of crystallographic group
    :param use_alphabet: if True, then symbolic matrices are created like
                                [a, b]
                                [c, d]
                         if False, then
                                [a_00, a_01]
                                [a_10, a_11]
    :param normalize_exp: if True, then solution of the linear system will be
                          normalized, i.e. all the independed variables will be
                          renamed into x_0, x_1, x_2 ...
    :param to_matrix: if True, then the result solutions will be transformed
                      into symbolic matrices instead of tuple of Expressions
    :param verbose: True to see the results on the fly.
    :param ignore_trivial: if True then solutions with zero determinant will be
                           ignored
    """

    G = gap(f'SpaceGroupIT({dim}, {n})')
    S = G.PointGroup()
    s_elements = [matrix(QQ, el) for el in S.AsList()]

    if verbose:
        print('\n====================================================')
        print(n, 'point group:', S)

        print('group elements:')
        print(*s_elements, sep='\n')
        print('\n----------------normalizers-------------------------')

    A, args = create_symbolic_matrix(dim=dim, use_alphabet=use_alphabet)

    found_solutions = set()
    for perm in permutations(range(len(s_elements))):

        # check whether permutation maps elements with same order
        check = True
        for i in range(len(s_elements)):
            if S.AsList()[i + 1].Order() != S.AsList()[perm[i] + 1].Order():
                check = False
                break
        if not check:
            continue

        # build conditions AX = YA for every X -> Y due to chosen permutation
        cond_perm = set()
        for i in range(len(s_elements)):
            cond_i = A * s_elements[i] - s_elements[perm[i]] * A
            for el in cond_i:
                cond_perm = cond_perm.union(set(el))

        eq = [cond == 0 for cond in cond_perm]
        for res in solve(eq, *args):
            if not res:
                if verbose: print(res)
                continue
            res = tuple(res)
            if normalize_exp:
                res = normalize_expressions(res, allowed=args)

            if res not in found_solutions and verbose:
                if to_matrix:
                    mtx = sol2matrix(res, dim=dim, use_alphabet=use_alphabet)
                    if ignore_trivial and mtx.det() == 0:
                        continue
                    print(mtx, end='\n\n')
                else:
                    print(res)

            found_solutions.add(res)

    if to_matrix:
        return [sol2matrix(solution, dim=dim, use_alphabet=use_alphabet) for solution in found_solutions]
    else:
        return found_solutions


def els_by_order(group_elements):
    res = defaultdict(list)

    for el in group_elements:
        res[el.order()].append(el)
    return res


def gens_mappings(gens, group_elements):
    orders = els_by_order(group_elements)
    res = []
    for gen in gens:
        res.append(orders[gen.order()])
    return res


def carthesian_wo_duplicates(*spaces):
    for res, _ in _carthesian_wo_duplicates(*spaces):
        yield res


def _carthesian_wo_duplicates(*spaces):
    if not spaces:
        yield [], set()
        return

    for res, used in _carthesian_wo_duplicates(*spaces[:-1]):
        for el in spaces[-1]:
            if str(el) in used:
                continue

            used.add(str(el))
            yield res + [el], used
            used.remove(str(el))


def get_point_group_gens(cryst_num, dim=3):
    group = gap(f"SpaceGroupIT({dim}, {cryst_num})")
    gens = group.PointGroup().MinimalGeneratingSet()
    return [matrix(QQ, el) for el in gens]


# fixme: ensure all works the same way as normalizers_v1.
#   Make tests for all the parameters
#   refactor maybe
# noinspection PyShadowingNames
def normalizers_v2(n, dim=2, verbose=False, use_alphabet=False,
                   normalize_exp=True, to_matrix=True, ignore_trivial=True):
    """Find normalizer of the PointGroup in GL(n, QQ).

    Tries to find normalizer as a solution of
    system of linear equations N*A_i = A_j*N, where A_i and A_j are
    two elements of the PointGroup with respect to a permutation. Since
    Normalizer is a group N, that satisfies condition NA = AN,
    we can check all the permutations of A and solve the respective system
    of linear equations.

    Second version, fast. Exploits the fact, that instead of full permutation of
    A we can just declare where to map generators of A.
    ----------------------------------------------------------------------------
    :param n: index of the Crystallographic group in the Gap CrystCat package
    :param dim: dimension of crystallographic group
    :param use_alphabet: if True, then symbolic matrices are created like
                                [a, b]
                                [c, d]
                         if False, then
                                [a_00, a_01]
                                [a_10, a_11]
    :param normalize_exp: if True, then solution of the linear system will be
                          normalized, i.e. all the independed variables will be
                          renamed into x_0, x_1, x_2 ...
    :param to_matrix: if True, then the result solutions will be transformed
                      into symbolic matrices instead of tuple of Expressions
    :param verbose: True to see the results on the fly.
    :param ignore_trivial: if True then solutions with zero determinant will be
                           ignored
    """
    S = MatrixGroup(get_point_group_gens(n, dim=dim))
    gens = S.gens()
    s_elements = list(S)
    possible_mappings = gens_mappings(gens, s_elements)

    if verbose:
        print('\n====================================================')
        print(n, 'point group:', S)

        print('group elements:')
        print(*s_elements, sep='\n')
        print('\n----------------normalizers-------------------------')

    A, args = create_symbolic_matrix(dim, use_alphabet=use_alphabet)
    found_solutions = set()

    for maps_to in carthesian_wo_duplicates(*possible_mappings):
        # if len(set(str(el) for el in maps_to)) < len(maps_to):
        #     continue

        # build conditions AX = YA for every X -> Y due to chosen permutation
        cond_perm = set()
        try:
            homm = S.hom(maps_to)
        except ValueError:
            continue
        S_cur = MatrixGroup(maps_to)
#         if not S.is_isomorphic(S_cur):
#             continue

        for X in s_elements:
            Y = homm(X)

            X, Y = matrix(QQ, X), matrix(QQ, Y)
            cond_i = A * X - Y * A
            for el in cond_i:
                cond_perm = cond_perm.union(set(el))

        eq = [cond == 0 for cond in cond_perm]
        for res in solve(eq, *args):
            if not res:
                if verbose: print(res)
                continue

            res = tuple(res)
            if normalize_exp:
                res = normalize_expressions(res, allowed=args)

            if res not in found_solutions and verbose:
                if to_matrix:
                    mtx = sol2matrix(res, dim=dim, use_alphabet=use_alphabet)
                    if ignore_trivial and mtx.det() == 0:
                        continue
                    print(mtx, end='\n\n')
                else:
                    print(res)

            found_solutions.add(res)

    if to_matrix:
        return [sol2matrix(solution, dim=dim, use_alphabet=use_alphabet) for solution in found_solutions]
    else:
        return found_solutions


# check whether matrix is simply factorizable
def check_div(A):
    x = var('x')

    n = len(list(A))
    test1 = A - x * matrix.identity(n)
    test2 = test1.T()

    # check rows
    for row in test1:
        if list(row).count(0) >= n - 1:
            return True

    # check columns
    for row in test2:
        if list(row).count(0) >= n - 1:
            return True
    return False


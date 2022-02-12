#!/usr/bin/env sage

from sage.all_notebook import *
from itertools import permutations

a = b = c = d = ...


def normalize_expressions(exps, allowed=None):
    """Replace all variables except allowed (default: a b c d) with {x0,x1,x2,x3...}"""
    extra_args = set()

    if allowed is None:
        allowed = var('a b c d')

    for exp in exps:
        print(exp)
        extra_args |= {el for el in exp.args() if el not in allowed}

    res = []
    for exp in exps:
        for i, el in enumerate(extra_args):
            exp = exp.subs({var(el): var(f'x{i}')})
        res.append(exp)
    return tuple(res)


def prepare_gap_env():
    gap('LoadPackage("Cryst");;Read("./el2word.g")')


def el2word(G, el, verbose=False):
    """Represent, given an element of a finite group G,
    as the word over generators."""

    prepare_gap_env()
    gap.set('G', G)
    gap.set('g', el)

    gens_dict = {f'g_{i}': el for i, el in enumerate(list(G.GeneratorsOfGroup()))}

    gens = '"' + '", "'.join(sorted(gens_dict.keys())) + '"'
    expression = f'FreeGroup({gens})'
    if verbose:
        print("Free group to be created:", expression)
        print("Translation of generators:", gens_dict)

    gap.set("F", expression)

    res = gap('ElementAsWordGeneratorsPointGroup(G, g, F)')
    return str(res).split("*"), gens_dict


def extend_permutation(G, el, gens, perm):
    """Extend given permutation on generators onto the entire group."""
    word, trans = el2word(G, el)

    if el == G.One():
        return el

    res = G.One()
    for g in word:
        power = 1
        if "^" in g:
            g, power = g.split('^')
            power = int(power)
        index = gens.index(trans[g])
        index_to = perm[index]

        g_to = gens[index_to] ** power
        res *= g_to
    return res


alphabet = 'abcdefghijklmnopqrstuvwxyz'


def create_symbolic_matrix(dim, use_alphabet=False):
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
    return A, set(args)


def normalizers(n, verbose=False):
    G = gap(f'SpaceGroupIT(2, {n})')
    S = G.PointGroup()
    s_elements = [matrix(QQ, el) for el in S.AsList()]

    if verbose:
        print('\n====================================================')
        print(n, 'point group:', S)

        print('group elements:')
        print(*s_elements, sep='\n')
        print('\n----------------normalizers-------------------------')

    var('a, b, c, d')
    A = matrix([[a, b], [c, d]])

    conditions = []
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
        for res in solve(eq, a, b, c, d):
            if not res:
                if verbose: print(res)
                continue

            res = normalize_expressions(res)
            if res not in found_solutions and verbose:
                print(res)

            found_solutions.add(res)

    return found_solutions

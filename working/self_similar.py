

# This file was *autogenerated* from the file ../working/self_similar.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_4 = Integer(4); _sage_const_2 = Integer(2)
from itertools import product

alphabet = 'abcdefghijklmnopqrstuvwxyz'

# ruff: noqa: F821
def to_L_basis(n, dim=_sage_const_3 ):
    """ Change basis of the crystallographic group to transform the lattice
    L to the Z^dim

    Parameters
    ----------
    n : int
        number of the crystallographic group
    dim : int
        dimension

    Returns
    -------
    MatrixGroup with new generators
    """
    G = gap(f'SpaceGroupOnLeftIT({dim}, {n})')
    gens = [matrix(QQ, el) for el in G.GeneratorsOfGroup()]

    v = matrix(QQ, G.TranslationBasis()).T
    trans = matrix(QQ, [_sage_const_0  for _ in range(dim)]).T
    conj = block_matrix(QQ, [[v, trans], [_sage_const_0 , _sage_const_1 ]])
    new_gens = [conj.inverse() * el * conj for el in gens]
    return gap.AffineCrystGroupOnLeft(new_gens)


def extend_names(names, gens, deep=_sage_const_4 , include_inverse=True, skip_trans=True):
    n = len(list(gens[_sage_const_0 ])) - _sage_const_1 

    res_names = names.copy()

    if include_inverse:
        res_gens = gens.copy()
        for el in gens:
            el_inv = el.inverse()
            if el_inv == el:
                continue
            name = names[str(el)]
            res_names[str(el_inv)] = f'{name}^(-1)'
            res_gens.append(el_inv)
        gens = res_gens

    res_rows = []
    for i in range(n):
        row = list(matrix.identity(n)[i]) + [_sage_const_0 ]
        res_rows.append(row)
    res_rows.append([_sage_const_0  for _ in range(n)] + [_sage_const_1 ])
    ident = matrix(QQ, res_rows)

    res_names[str(ident)] = 'e'

    if skip_trans:
        res_gens = []
        for el in gens:
            if el[:n, :n] == matrix.identity(n):
                continue
            res_gens.append(el)
        gens = res_gens

    for prod in product(gens, repeat=deep):
        res = prod[_sage_const_0 ]
        res_name = res_names[str(res)]

        for el in prod[_sage_const_1 :]:
            res *= el
            el_name = res_names[str(el)]
            res_name += el_name

        if str(res) not in res_names:
            res_names[str(res)] = res_name

    return res_names


def get_name(el, names):
    reduced_el = copy(el)
    n = len(list(el)) - _sage_const_1 
    res = ''
    for i in range(n):
        reduced_el[i, -_sage_const_1 ] = el[i, -_sage_const_1 ] - int(el[i, -_sage_const_1 ])

        cur_trans = [_sage_const_0  for _ in range(n)]
        cur_trans[i] = _sage_const_1 
        res_rows = []
        for i in range(n):
            row = list(matrix.identity(n)[i]) + [cur_trans[i]]
            res_rows.append(row)
        res_rows.append([_sage_const_0  for _ in range(n)] + [_sage_const_1 ])
        cur_trans = matrix(QQ, res_rows)
        cur_trans_name = names[str(cur_trans)]

        power = int(el[i, -_sage_const_1 ])
        if power > _sage_const_0 :
            power = str(power)
        elif power < _sage_const_0 :
            power = f'({power})'
        else:
            continue
        res += f'{cur_trans_name}^{power}'
    if str(reduced_el) in names:
        name = (names[str(reduced_el)] + res)
        if name == 'e':
            return 'e'
        else:
            return name.replace('e', '').replace('^1', '')
    else:
        return str(el)


def self_similar(n, T, dim=_sage_const_2 , verbose=False,
                 gen_alphabet=False, safe=True, change_basis=False, deep=_sage_const_4 ):
    """ Construct self-similar action for a crystallographic group
    given element of affine group that is conjugation for virtual
    endomorphism construction.

    Parameters
    ----------
    n : int
        a number of crystallographic group from the Gap package
    T : matrix of dimension `dim` + 1
        an element of the affine group A(dim) that can be represented as
        (M + t), where M is a `dim`x`dim` matrix and t is a vector that represents
        translation. Should be given in a matrix form, i.e. (M + t) is a
        block matrix:
                         _           _
                        ||     |     ||
                        ||  M  |  t  ||
                        ||_____|_____||
                        ||  0  |  1  ||
                         -           -
        dim : int
        dimension of euclidean space, where we consider a crystallographic
        group
    verbose : bool
        True to show auxiliary messages
    gen_alphabet : bool
        True to use alphabet for generators instead of a_i
    safe : bool
        if True, then function raises error if `T` doesn't generate
        virtual endomorphism.
    deep: int
        specifies maximal length of the precomputed words of generators to
        compute explicit formula of the self-similar action.
    Returns
    -------
    dict { (a, i): [j, b] }
        a self-similar action
    """
    if change_basis:
        G = to_L_basis(n, dim)
    else:
        G = gap(f'SpaceGroupOnLeftIT({dim}, {n})')

    gens_G = G.GeneratorsOfGroup()
    gens_G = [matrix(QQ, el) for el in gens_G]
    if verbose:
        print("=====================================================================")

    # check whether there exists virtual endomorphism
    conj_els = []
    for el in gens_G:

        # conjugate in the opposite direction in order to get subgroup of
        # finite index, i.e. apply \phi^{-1} (el)
        conj = T.inverse() * el * T
        if verbose:
            print("\nconjugate el:")
            print(conj)
            print("conj in G:", conj in G)
        if conj not in G and safe:
            raise ValueError("Bad matrix T, there is no virtual endomorphism")
        elif conj not in G:
            print("Bad matrix T, there is no virtual endomorphism")
            return -_sage_const_1 

        conj_els.append(conj)

    # create subgroup as image of virtual endomorphism
    H = G.Subgroup(conj_els)

    if verbose:
        print("----------------------------------------------------")
        print("Index of subgroup H:", G.Index(H))

    trans = G.RightTransversal(H)
    trans_els = [matrix(QQ, el) for el in trans.AsList()]
    if verbose:
        print("Transversal:")
        print(*trans_els, sep='\n\n')

    if gen_alphabet:
        names_G = {str(el): letter for el, letter in zip(gens_G, alphabet)}
    else:
        names_G = {str(el): f"a_{i}" for i, el in enumerate(gens_G, _sage_const_1 )}

    extended_names = extend_names(names_G, gens_G, deep=deep)

    # create self-similar action
    res_map = {}
    for a in gens_G:

        # NOTE: ENUMERATION STARTS FROM 1
        for i, d_i in enumerate(trans_els, _sage_const_1 ):
            adi = a * d_i

            if verbose:
                print(f"{names_G[str(a)]}d_{i}:")
                print(adi, end='\n\n')

            # (a * d_i)^{-1} * (a * d_i) = e \in H
            # ==> d_j^{-1} = (a * d_i)^{-1}
            #
            # firstly, find coset for (a * d_i)^{-1}
            d_j_index = trans.PositionCanonical(adi.inverse())

            # then get d_j^{-1}^{-1}
            d_j = trans_els[int(d_j_index) - _sage_const_1 ].inverse()

            if d_j.inverse() * a * d_i not in H:
                raise ValueError(f"This shouln't happen--wrong coset: {(a, d_i, d_j)}")

            # conjugation in the right direction, i.e. apply \phi(d_j^{-1} a d_i)
            tmp_res = T * (d_j.inverse() * a * d_i) * T.inverse()
            res_map[(names_G[str(a)], i)] = (int(d_j_index),
                                             get_name(tmp_res, extended_names))
    return res_map, extended_names


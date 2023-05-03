
from itertools import product

alphabet = 'abcdefghijklmnopqrstuvwxyz'


def to_L_basis(n, dim=3):
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
    trans = matrix(QQ, [0 for _ in range(dim)]).T
    conj = block_matrix(QQ, [[v, trans], [0, 1]])
    new_gens = [conj.inverse() * el * conj for el in gens]
    return gap.AffineCrystGroupOnLeft(new_gens)


def extend_names(names, gens, deep=4, include_inverse=True, skip_trans=True):
    n = len(list(gens[0])) - 1

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
        row = list(matrix.identity(n)[i]) + [0]
        res_rows.append(row)
    res_rows.append([0 for _ in range(n)] + [1])
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
        res = prod[0]
        res_name = res_names[str(res)]

        for el in prod[1:]:
            res *= el
            el_name = res_names[str(el)]
            res_name += el_name

        if str(res) not in res_names:
            res_names[str(res)] = res_name

    return res_names


def get_name(el, names):
    reduced_el = copy(el)
    n = len(list(el)) - 1
    res = ''
    for i in range(n):
        reduced_el[i, -1] = el[i, -1] - int(el[i, -1])

        cur_trans = [0 for _ in range(n)]
        cur_trans[i] = 1
        res_rows = []
        for i in range(n):
            row = list(matrix.identity(n)[i]) + [cur_trans[i]]
            res_rows.append(row)
        res_rows.append([0 for _ in range(n)] + [1])
        cur_trans = matrix(QQ, res_rows)
        cur_trans_name = names[str(cur_trans)]

        power = int(el[i, -1])
        if power > 0:
            power = str(power)
        elif power < 0:
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


def self_similar(n, T, dim=2, verbose=False,
                 gen_alphabet=False, safe=True, change_basis=False, deep=4):
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
            return -1

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
        names_G = {str(el): f"a_{i}" for i, el in enumerate(gens_G, 1)}

    extended_names = extend_names(names_G, gens_G, deep=deep)

    # create self-similar action
    res_map = {}
    for a in gens_G:

        # NOTE: ENUMERATION STARTS FROM 1
        for i, d_i in enumerate(trans_els, 1):
            adi = a * d_i

            # (a * d_i)^{-1} * (a * d_i) = e \in H
            # ==> d_j^{-1} = (a * d_i)^{-1}
            #
            # firstly, find coset for (a * d_i)^{-1}
            d_j_index = trans.PositionCanonical(adi.inverse())

            # then get d_j^{-1}^{-1}
            d_j = trans_els[int(d_j_index) - 1].inverse()

            if d_j.inverse() * a * d_i not in H:
                raise ValueError(f"This shouln't happen--wrong coset: {(a, d_i, d_j)}")

            # conjugation in the right direction, i.e. apply \phi(d_j^{-1} a d_i)
            tmp_res = T * (d_j.inverse() * a * d_i) * T.inverse()
            res_map[(names_G[str(a)], i)] = (int(d_j_index),
                                             get_name(tmp_res, extended_names))
    return res_map, extended_names

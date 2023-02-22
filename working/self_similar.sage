

alphabet = 'abcdefghijklmnopqrstuvwxyz'


def self_similar(n, T, dim=2, verbose=False, gen_alphabet=False, safe=True):
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

    Returns
    -------
    dict { (a, i): [j, b] }
        a self-similar action
    """
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
            res_map[(names_G[str(a)], i)] = (int(d_j_index), tmp_res)
    return res_map



alphabet = 'abcdefghijklmnopqrstuvwxyz'


def self_similar(n, T, dim=2, verbose=False, gen_alphabet=False, safe=True):
    G = gap(f'SpaceGroupIT({dim}, {n})')
    gens_G = G.GeneratorsOfGroup()
    gens_G = [matrix(QQ, el) for el in gens_G]
    if verbose:
        print("=====================================================================")

    # check whether there exists virtual endomorphism
    conj_els = []
    for el in gens_G:
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
    gens_H = [matrix(QQ, el) for el in H.GeneratorsOfGroup()]
    trans_els = [matrix(QQ, el) for el in trans.AsList()]
    trans_inv_els = [el.inverse() for el in trans_els]
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
        for i, d_i in enumerate(trans_inv_els, 1):
            adi = a * d_i
            d_j_index = trans.PositionCanonical(adi.inverse())
            d_j = trans_inv_els[int(d_j_index) - 1]

            if d_j.inverse() * a * d_i not in H:
                print("Error, wrong coset")
                print(a, d_i, d_j)
                return -1

            tmp_res = T * (d_j.inverse() * a * d_i) * T.inverse()
            res_map[(names_G[str(a)], i)] = (int(d_j_index), tmp_res)
    return res_map
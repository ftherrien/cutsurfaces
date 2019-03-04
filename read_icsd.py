def icsd_cif_a(filename, make_primitive=True):
    """ Reads lattice from the ICSD \*cif files.

        It will not work in the case of other \*cif.
        It is likely to produce wrong output if the site occupations are fractional.
        If the occupation is > 0.5 it will treat it as 1 and
        in the case occupation < 0.5 it will treat it as 0 and
        it will accept all occupation = 0.5 as 1 and create a mess!
    """
    import re
    from copy import deepcopy
    from os.path import basename
    from numpy.linalg import norm
    from numpy import array, transpose
    from numpy import pi, sin, cos, sqrt, dot
    from sys import version_info

    if version_info[0] >= 3:
        lines = open(filename, 'r',encoding='latin1').readlines()
    else:
        lines = open(filename, 'r').readlines()

    sym_big = 0
    sym_end = 0
    pos_big = 0
    pos_end = 0

    for l in lines:
        x = l.split()
        if len(x) > 0:
                    # CELL
            if x[0] == '_cell_length_a':
                if '(' in x[-1]:
                    index = x[-1].index('(')
                else:
                    index = len(x[-1])
                a = float(x[-1][:index])

            if x[0] == '_cell_length_b':
                if '(' in x[-1]:
                    index = x[-1].index('(')
                else:
                    index = len(x[-1])
                b = float(x[-1][:index])

            if x[0] == '_cell_length_c':
                if '(' in x[-1]:
                    index = x[-1].index('(')
                else:
                    index = len(x[-1])
                c = float(x[-1][:index])

            if x[0] == '_cell_angle_alpha':
                if '(' in x[-1]:
                    index = x[-1].index('(')
                else:
                    index = len(x[-1])
                alpha = float(x[-1][:index])

            if x[0] == '_cell_angle_beta':
                if '(' in x[-1]:
                    index = x[-1].index('(')
                else:
                    index = len(x[-1])
                beta = float(x[-1][:index])

            if x[0] == '_cell_angle_gamma':
                if '(' in x[-1]:
                    index = x[-1].index('(')
                else:
                    index = len(x[-1])
                gamma = float(x[-1][:index])

        if len(x) > 0 and x[0] == '_symmetry_Int_Tables_number':
            spg = int(x[1])
                
                
        # SYMMETRY OPERATIONS

        if len(x) > 0 and x[0] == '_symmetry_equiv_pos_as_xyz':
            sym_big = lines.index(l)

        if len(x) > 0 and x[0] == '_atom_type_symbol':
            sym_end = lines.index(l)

        ## FT: reads oxydation states
        # OXYDATION STATE

        if len(x) > 0 and x[0] == '_atom_site_label':
            ox_end = lines.index(l)
            
        # WYCKOFF POSITIONS

        if len(x) > 0 and x[0] == '_atom_site_attached_hydrogens':
            pos_big = lines.index(l)

        if len(x) > 0 and x[0] == '_atom_site_B_iso_or_equiv':
            pos_big = lines.index(l)

        if len(x) > 0 and x[0] == '_atom_site_U_iso_or_equiv':
            pos_big = lines.index(l)

        if len(x) > 0 and x[0] == '_atom_site_0_iso_or_equiv':
            pos_big = lines.index(l)

        # if pos_end == 0 and l in ['\n', '\r\n'] and lines.index(l) > pos_big:
        if pos_end == 0 and pos_big > 0 \
                and (l in ['\n', '\r\n'] or l.startswith('#')) \
                and lines.index(l) > pos_big:
            pos_end = lines.index(l)

    # _symmetry_equiv_pos_* lines are like:
    #     1     'x, x-y, -z+1/2'

    symm_ops = ['(' + x.split()[1][1:] + x.split()[2] + x.split()[3][:-1] + ')'
                for x in lines[sym_big + 1:sym_end - 1]]
    
    # ['(x,x-y,-z+1/2)', '(-x+y,y,-z+1/2)', ...]

    # Insert decimal points after integers
    symm_ops = [re.sub(r'(\d+)', r'\1.', x) for x in symm_ops]
    # ['(x,x-y,-z+1./2.)', '(-x+y,y,-z+1./2.)', ...]

    # _atom_site_* lines are like:
    #   Mo1 Mo4+ 2 c 0.3333 0.6667 0.25 1. 0
    ## FT: replaced [0] by [1] to take the ion name instead (ex: Mo4+ instead of Mo1)
    wyckoff = [[x.split()[1], [x.split()[4], x.split()[5], x.split()[6]], x.split()[7]]
               for x in lines[pos_big + 1:pos_end]]
    # [['Mo1', ['0.3333', '0.6667', '0.25'], '1.'], ['S1', ['0.3333', '0.6667', '0.621(4)'], '1.']]

    wyckoff = [w for w in wyckoff if int(float(w[-1][:4]) + 0.5) != 0]
    # [['Mo1', ['0.3333', '0.6667', '0.25'], '1.'], ['S1', ['0.3333', '0.6667', '0.621(4)'], '1.']]

    ##FT: reading the proper oxidation states
    oxidation = [[x.split()[0], int(round(float(x.split()[1]),0))] for x in lines[sym_end+2:ox_end-1]]

    # Setting up a good wyckoff list

    for w in wyckoff:

        ## FT: Not stripping anymore to keep different oxydation states different
        # Strip trailing numerals from w[0] == 'Mo1'
        # pom = 0
        # for i in range(len(w[0])):
        #     try:
        #         int(w[0][i])
        #         if pom == 0:
        #             pom = i
        #     except:
        #         pass

        # w[0] = w[0][:pom]

        # Strip trailing standard uncertainty, if any, from w[1], ..., w[3]
        for i in range(3):
            if '(' in w[1][i]:
                index = w[1][i].index('(')
            else:
                index = len(w[1][i])
            w[1][i] = float(w[1][i][:index])

        # Delete w[4]
        del w[-1]
    ##########################################

    # List of unique symbols ["Mo", "S"]
    symbols = list({w[0] for w in wyckoff})

    # List of position vectors for each symbol
    positions = [[] for i in range(len(symbols))]

    for w in wyckoff:
        symbol = w[0]
        x, y, z = w[1][0], w[1][1], w[1][2]
        for i in range(len(symm_ops)):
            # Set pom = new position based on symmetry transform
            pom = list(eval(symm_ops[i]))
            # [0.3333, -0.3334, 0.25]

            # Move positions to range [0,1]:
            for j in range(len(pom)):
                if pom[j] < 0.:
                    pom[j] = pom[j] + 1.
                if pom[j] >= 0.999:
                    pom[j] = pom[j] - 1.
            # [0.3333, 0.6666, 0.25]

            # If pom is not in positions[symbol], append pom
            if not any(norm(array(u) - array(pom)) < 0.01 for u in positions[symbols.index(symbol)]):
                ix = symbols.index(symbol)
                positions[ix].append(pom)

    ################ CELL ####################

    a1 = a * array([1., 0., 0.])
    a2 = b * array([cos(gamma * pi / 180.), sin(gamma * pi / 180.), 0.])
    c1 = c * cos(beta * pi / 180.)
    c2 = c / sin(gamma * pi / 180.) * (-cos(beta * pi / 180.) *
                                       cos(gamma * pi / 180.) + cos(alpha * pi / 180.))
    a3 = array([c1, c2, sqrt(c**2 - (c1**2 + c2**2))])
    cell = array([a1, a2, a3])
    ##########################################

    from pylada.crystal import Structure, primitive

    structure = Structure(
        transpose(cell),
        scale=1,
        name=basename(filename), group = spg)

    for i in range(len(symbols)):
        # crystal/read: i:  0   symbol:  Mo   len position:  2

        for j in range(len(positions[i])):
            atpos = dot(transpose(cell), positions[i][j])
            #  j:  0   pos:  [0.3333, 0.6666000000000001, 0.25]
            #  atpos:  [  6.32378655e-16   1.81847148e+00   3.07500000e+00]

            ## FT: Finds the corresponding oxidation
            for o in oxidation:
                if o[0] == symbols[i]:
                    ox = o[1]
                    break

            structure.add_atom(atpos[0], atpos[1], atpos[2], symbols[i], ox = ox)


    if make_primitive:
        prim = primitive(structure)
    else:
        prim = deepcopy(structure)

    return prim

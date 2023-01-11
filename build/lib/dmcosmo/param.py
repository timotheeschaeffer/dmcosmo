"""
External Parameters
"""



class Bunch(object):
    """
    translates dic['name'] into dic.name 
    """

    def __init__(self, data):
        self.__dict__.update(data)


def PS_par():
    par = {
        "filter": 'tophat',           # tophat, sharpk or smoothk
        "c": 1,                    # scale to halo mass relation (1 for tophat, 2.5 for sharp-k, 3 for smooth-k)
        "q" : 1,                     # q for f(nu) [0.707,1,1] for [ST,smoothk or sharpk,PS]
        "p" : 0.3,                   # p for f(nu) [0.3,0.3,0] for [ST,smoothk or sharpk,PS]
        "delta_c" : 1.686,
        "A" : 0.322,                    # A = 0.322 except 0.5 for PS Spherical collapse (to double check)
        "beta" : 3,             # smooth-k parameter
        }
    return Bunch(par)

def code_par():
    par = {
        "m_min" :1e4,
        "m_max" : 1e16,
        "Mbin" : 300,
        'z' : 0,  #output z value
        }
    return Bunch(par)


def cosmo_par():
    par = {
    'Om' : 0.321,
    'Ob' : 0.04945,
    'Ol' : 0.679,
    'rho_c' : 2.775e11,
    'h' : 0.6688,
    'ps': './../files/PCDM.dat',
    }
    return Bunch(par)




def par():
    par = Bunch({
        "PS": PS_par(),
        "code": code_par(),
        "cosmo" : cosmo_par(),
        })
    return par

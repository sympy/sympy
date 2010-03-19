"""Configuration utilities for polynomial manipulation algorithms. """

_default_config = {
    'USE_COLLINS_RESULTANT'  : False,
    'USE_SIMPLIFY_GCD'       : True,
    'USE_HEU_GCD'            : True,

    'USE_CYCLOTOMIC_FACTOR'  : True,

    'EEZ_RESTART_IF_NEEDED'  : True,
    'EEZ_NUMBER_OF_CONFIGS'  : 3,
    'EEZ_NUMBER_OF_TRIES'    : 5,
    'EEZ_MODULUS_STEP'       : 2,

    'GF_IRRED_METHOD'        : 'rabin',
    'GF_FACTOR_METHOD'       : 'zassenhaus',
}

_current_config = {}

def setup(key, value=None):
    """Assign a value to (or reset) a configuration item. """
    if value is not None:
        _current_config[key] = value
    else:
        _current_config[key] = _default_config[key]

def query(key):
    """Ask for a value of the given configuration item. """
    return _current_config.get(key, None)

def configure():
    """Initialized configuration of polys module. """
    from os import getenv

    for key, default in _default_config.items():
        value = getenv('SYMPY_' + key)

        if value is not None:
            try:
                _current_config[key] = eval(value)
            except NameError:
                _current_config[key] = value
        else:
            _current_config[key] = default

configure()

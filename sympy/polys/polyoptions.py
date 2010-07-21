"""Options manager for :class:`Poly` and public API functions. """

from sympy.core import S, Basic, sympify

from sympy.polys.polyerrors import (
    PolynomialError,
    GeneratorsError,
    OptionError,
    FlagError,
)

from sympy.polys.monomialtools import monomial_key

from sympy.polys.domains import FF, GF, ZZ, QQ, RR, EX

from sympy.ntheory import isprime

import re

class Option(object):
    """ """

    option = None

    is_Flag = False

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        """ """
        return None

    @classmethod
    def preprocess(cls, option):
        """ """
        return None

    @classmethod
    def postprocess(cls, options):
        """ """
        pass

class Flag(Option):
    """ """

    is_Flag = True

class BooleanOption(Option):
    """ """

    @classmethod
    def preprocess(cls, value):
        if value is True or value is False or value is 1 or value is 0:
            return bool(value)
        else:
            raise OptionError("'%s' must have a boolean value assigned, got %s" % (cls.option, value))

class OptionType(type):
    """ """

    def __init__(cls, *args, **kwargs):
        """ """
        @property
        def getter(self):
            try:
                return self[cls.option]
            except KeyError:
                return cls.default()

        setattr(Options, cls.option, getter)
        Options.__options__[cls.option] = cls

class Options(dict):
    """ """

    __options__ = {}

    def __init__(self, gens, args, flags=None, strict=False):
        """ """
        dict.__init__(self)

        if gens and args.get('gens', ()):
            raise OptionError("both '*gens' and keyword argument 'gens' supplied")
        elif gens:
            args = dict(args)
            args['gens'] = gens

        for option, value in args.iteritems():
            try:
                cls = self.__options__[option]
            except KeyError:
                raise OptionError("'%s' is not a valid option" % option)

            if issubclass(cls, Flag):
                if flags is None or option not in flags:
                    if strict:
                        raise OptionError("'%s' flag is not allowed in this context" % option)

            if value is not None:
                self[option] = cls.preprocess(value)

        for option in self.keys():
            cls = self.__options__[option]

            for require_option in cls.requires:
                if self.get(require_option) is None:
                    raise OptionError("'%s' option is only allowed together with '%s'" % (option, require_option))

            for exclude_option in cls.excludes:
                if self.get(exclude_option) is not None:
                    raise OptionError("'%s' option is not allowed together with '%s'" % (option, exclude_option))

        for option in self.keys():
            self.__options__[option].postprocess(self)

    @property
    def args(self):
        args = {}

        for option, value in self.iteritems():
            if value is not None and option != 'gens':
                cls = self.__options__[option]

                if not issubclass(cls, Flag):
                    args[option] = value

        return args

    @property
    def options(self):
        options = {}

        for option, cls in self.__options__.iteritems():
            if not issubclass(cls, Flag):
                options[option] = getattr(self, option)

        return options

    @property
    def flags(self):
        flags = {}

        for option, cls in self.__options__.iteritems():
            if issubclass(cls, Flag):
                flags[option] = getattr(self, option)

        return flags

class Expand(BooleanOption):
    """ """

    __metaclass__ = OptionType

    option = 'expand'

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        return True

class Gens(Option):
    """ """

    __metaclass__ = OptionType

    option = 'gens'

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        return ()

    @classmethod
    def preprocess(cls, gens):
        if isinstance(gens, Basic):
            gens = (gens,)
        elif len(gens) == 1 and hasattr(gens[0], '__iter__'):
            gens = gens[0]

        if len(set(gens)) != len(gens):
            raise GeneratorsError("duplicated generators: %s" % str(gens))

        return tuple(gens)

class Wrt(Option):
    """ """

    __metaclass__ = OptionType

    option = 'wrt'

    requires = []
    excludes = []

    _re_split = re.compile(" |,")

    @classmethod
    def preprocess(cls, wrt):
        if isinstance(wrt, Basic):
            return [str(wrt)]
        elif isinstance(wrt, str):
            return [ gen.strip() for gen in cls._re_split.split(wrt) ]
        elif hasattr(wrt, '__getitem__'):
            return list(map(str, wrt))
        else:
            raise OptionError("invalid argument for 'wrt' option")

class Sort(Option):
    """ """

    __metaclass__ = OptionType

    option = 'sort'

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        return []

    @classmethod
    def preprocess(cls, sort):
        if isinstance(sort, str):
            return [ gen.strip() for gen in sort.split('>') ]
        elif hasattr(sort, '__getitem__'):
            return list(map(str, sort))
        else:
            raise OptionError("invalid argument for 'sort' option")

class Order(Option):
    """ """

    __metaclass__ = OptionType

    option = 'order'

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        return monomial_key('lex')

    @classmethod
    def preprocess(cls, order):
        return monomial_key(order)

class Field(BooleanOption):
    """ """

    __metaclass__ = OptionType

    option = 'field'

    requires = []
    excludes = ['domain', 'split', 'gaussian', 'extension']

class Greedy(BooleanOption):
    """ """

    __metaclass__ = OptionType

    option = 'greedy'

    requires = []
    excludes = ['domain', 'split', 'gaussian', 'extension', 'modulus', 'symmetric']

class Domain(Option):
    """ """

    __metaclass__ = OptionType

    option = 'domain'

    requires = []
    excludes = ['field', 'greedy', 'split', 'gaussian', 'extension']

    _re_finitefield = re.compile("^(FF|GF)\((\d+)\)$")
    _re_polynomial  = re.compile("^(Z|ZZ|Q|QQ)\[(.+)\]$")
    _re_fraction    = re.compile("^(Z|ZZ|Q|QQ)\((.+)\)$")
    _re_algebraic   = re.compile("^(Q|QQ)\<(.+)\>$")

    @classmethod
    def preprocess(cls, domain):
        if not isinstance(domain, str):
            return domain
        else:
            if domain in ['Z', 'ZZ']:
                return ZZ

            if domain in ['Q', 'QQ']:
                return QQ

            if domain in ['R', 'RR']:
                return RR

            if domain == 'EX':
                return EX

            r = cls._re_finitefield.match(domain)

            if r is not None:
                _, modulus = r.groups()
                return FF(int(modulus))

            r = cls._re_polynomial.match(domain)

            if r is not None:
                ground, gens = r.groups()

                gens = map(sympify, gens.split(','))

                if ground in ['Z', 'ZZ']:
                    return ZZ.poly_ring(*gens)
                else:
                    return QQ.poly_ring(*gens)

            r = cls._re_fraction.match(domain)

            if r is not None:
                ground, gens = r.groups()

                gens = map(sympify, gens.split(','))

                if ground in ['Z', 'ZZ']:
                    return ZZ.frac_field(*gens)
                else:
                    return QQ.frac_field(*gens)

            r = cls._re_algebraic.match(domain)

            if r is not None:
                gens = map(sympify, r.groups()[1].split(','))
                return QQ.algebraic_field(*gens)

            raise OptionError('expected a valid domain specification, got %s' % domain)

    @classmethod
    def postprocess(cls, options):
        if 'gens' in options and options['domain'].is_Composite:
            if set(options['domain'].gens) & set(options['gens']):
                raise GeneratorsError("ground domain and generators interferes together")

class Split(BooleanOption):
    """ """

    __metaclass__ = OptionType

    option = 'split'

    requires = []
    excludes = ['field', 'greedy', 'domain', 'gaussian', 'extension', 'modulus', 'symmetric']

    @classmethod
    def postprocess(cls, options):
        raise NotImplementedError("'split' option is not implemented yet")

class Gaussian(BooleanOption):
    """ """

    __metaclass__ = OptionType

    option = 'gaussian'

    requires = []
    excludes = ['field', 'greedy', 'domain', 'split', 'extension', 'modulus', 'symmetric']

    @classmethod
    def postprocess(cls, options):
        if options['gaussian'] is True:
            options['extension'] = set([S.ImaginaryUnit])
            Extension.postprocess(options)

class Extension(Option):
    """ """

    __metaclass__ = OptionType

    option = 'extension'

    requires = []
    excludes = ['field', 'greedy', 'domain', 'split', 'gaussian', 'modulus', 'symmetric']

    @classmethod
    def preprocess(cls, extension):
        if extension is True or extension is 1:
            return bool(extension)
        elif extension is False or extension is 0:
            raise OptionError("'False' is an invalid argument for 'extension'")
        else:
            if not hasattr(extension, '__iter__'):
                extension = set([extension])
            else:
                if not extension:
                    extension = None
                else:
                    extension = set(extension)

            return extension

    @classmethod
    def postprocess(cls, options):
        if options['extension'] is not True:
            options['domain'] = QQ.algebraic_field(*options['extension'])

class Modulus(Option):
    """ """

    __metaclass__ = OptionType

    option = 'modulus'

    requires = []
    excludes = ['greedy', 'split', 'domain', 'gaussian', 'extension']

    @classmethod
    def preprocess(cls, modulus):
        modulus = sympify(modulus)

        if modulus.is_Integer and modulus > 0:
            return int(modulus)
        else:
            raise OptionError("'modulus' must a positive integer, got %s" % modulus)

    @classmethod
    def postprocess(cls, options):
        options['domain'] = FF(options['modulus'], options.get('symmetric', True))

class Symmetric(BooleanOption):
    """ """

    __metaclass__ = OptionType

    option = 'symmetric'

    requires = ['modulus']
    excludes = ['greedy', 'domain', 'split', 'gaussian', 'extension']

class Strict(BooleanOption):
    """ """

    __metaclass__ = OptionType

    option = 'strict'

    @classmethod
    def default(cls):
        return True

class Auto(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'auto'

    @classmethod
    def default(cls):
        return True

class Frac(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'frac'

    @classmethod
    def default(cls):
        return False

class Formal(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'formal'

    @classmethod
    def default(cls):
        return False

class Polys(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'polys'

class Include(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'include'

    @classmethod
    def default(cls):
        return False

class Monic(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'monic'

    @classmethod
    def default(cls):
        return True

class All(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'all'

    @classmethod
    def default(cls):
        return False

class Gen(Flag):
    """ """

    __metaclass__ = OptionType

    option = 'gen'

    requires = []
    excludes = []

    @classmethod
    def default(cls):
        return 0

    @classmethod
    def preprocess(cls, gen):
        if isinstance(gen, (Basic, int)):
            return gen
        else:
            raise OptionError("invalid argument for 'gen' option")

def build_options(gens, args=None):
    """Construct options from keyword arguments or ... options. """
    if args is None:
        gens, args = (), gens

    if len(args) != 1 or 'opt' not in args or gens:
        return Options(gens, args)
    else:
        return args['opt']

def allowed_flags(args, flags):
    """ """
    flags = set(flags)

    for arg in args.iterkeys():
        try:
            if Options.__options__[arg].is_Flag and not arg in flags:
                raise FlagError("'%s' flag is not allowed in this context" % arg)
        except KeyError:
            raise OptionError("'%s' is not a valid option" % arg)

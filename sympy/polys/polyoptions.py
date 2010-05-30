""" """

from sympy.core import S, Basic, sympify

from sympy.polys.polyerrors import (
    GeneratorsError,
    OptionError,
)

from sympy.polys.monomialtools import monomial_key

from sympy.ntheory import isprime

import re

from sympy.polys.algebratools import Algebra, ZZ, QQ, RR, EX

class Option(object):
    """ """

    option = None

    requires = []
    excludes = []

    default = None

    @classmethod
    def preprocess(cls, option):
        """ """
        return None

    @classmethod
    def postprocess(cls, options):
        """ """
        pass

class BooleanOption(Option):
    """ """

    @classmethod
    def preprocess(cls, value):
        if value in (True, False, 1, 0):
            return value
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
                return cls.default

        setattr(Options, cls.option, getter)
        Options.__options__[cls.option] = cls

class Flag(object):
    """ """

class Options(dict):
    """ """

    __options__ = {}

    def __init__(self, gens, args, flags=None):
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

class Expand(Option):
    __metaclass__ = OptionType

    option = 'expand'

    requires = []
    excludes = []

    default = None

    @classmethod
    def preprocess(cls, expand):
        if expand in (True, False, 1, 0):
            return expand
        elif isinstance(expand, dict):
            raise NotImplementedError
        else:
            raise OptionError("'expand' ...")

    @classmethod
    def postprocess(cls, options):
        pass

class Gens(Option):
    __metaclass__ = OptionType

    option = 'gens'

    requires = []
    excludes = []

    default = ()

    @classmethod
    def preprocess(cls, gens):
        if len(gens) == 1 and hasattr(gens[0], '__iter__'):
            gens = gens[0]

        if len(set(gens)) != len(gens):
            raise GeneratorsError("duplicated generators: %s" % str(gens))

        return tuple(gens)

    @classmethod
    def postprocess(cls, options):
        pass

class Wrt(Option):
    __metaclass__ = OptionType

    option = 'wrt'

    requires = []
    excludes = []

    default = None

    @classmethod
    def preprocess(cls, wrt):
        if isinstance(wrt, Basic):
            return [str(wrt)]
        elif isinstance(wrt, str):
            return [ gen.strip() for gen in wrt.split(",") ]
        elif hasattr(wrt, '__getitem__'):
            return list(map(str, wrt))
        else:
            raise OptionError("invalid argument for 'wrt' option")

    @classmethod
    def postprocess(cls, options):
        pass

class Sort(Option):
    __metaclass__ = OptionType

    option = 'sort'

    requires = []
    excludes = []

    default = None

    @classmethod
    def preprocess(cls, sort):
        if isinstance(sort, str):
            return sort
        else:
            raise OptionError("invalid argument for 'sort' option")

    @classmethod
    def postprocess(cls, options):
        pass

class Order(Option):
    __metaclass__ = OptionType

    option = 'order'

    requires = []
    excludes = []

    default = None

    @classmethod
    def preprocessder(cls, order):
        return monomial_key(order)

    @classmethod
    def postprocess(cls, options):
        raise NotImplementedError("'order' keyword is not implemented yet")

class Field(Option):
    __metaclass__ = OptionType

    option = 'field'

    requires = []
    excludes = ['domain', 'split', 'gaussian', 'extension', 'modulus', 'symmetric']

    default = None

    @classmethod
    def preprocess(cls, field):
        if isinstance(field, bool):
            return field
        else:
            raise OptionError("invalid argument for 'field' option")

    @classmethod
    def postprocess(cls, options):
        pass

class Greedy(Option):
    __metaclass__ = OptionType

    option = 'greedy'

    requires = []
    excludes = ['domain', 'split', 'gaussian', 'extension', 'modulus', 'symmetric']

    default = None

    @classmethod
    def preprocess(cls, greedy):
        if isinstance(greedy, bool):
            return greedy
        else:
            raise OptionError("invalid argument for 'greedy' option")

    @classmethod
    def postprocess(cls, options):
        pass

class Domain(Option):
    __metaclass__ = OptionType

    option = 'domain'

    requires = []
    excludes = ['field', 'greedy', 'split', 'gaussian', 'extension']

    default = None

    _re_polynomial = re.compile("^(Z|ZZ|Q|QQ)\[(.+)\]$")
    _re_fraction   = re.compile("^(Z|ZZ|Q|QQ)\((.+)\)$")
    _re_algebraic  = re.compile("^(Q|QQ)\<(.+)\>$")

    @classmethod
    def preprocess(cls, domain):
        if isinstance(domain, Algebra):
            return domain
        elif isinstance(domain, str):
            if domain in ['Z', 'ZZ']:
                return ZZ

            if domain in ['Q', 'QQ']:
                return QQ

            if domain in ['R', 'RR']:
                return RR

            if domain == 'EX':
                return EX

            r = re.match(cls._re_polynomial, domain)

            if r is not None:
                ground, gens = r.groups()

                gens = map(sympify, gens.split(','))

                if ground in ['Z', 'ZZ']:
                    return ZZ.poly_ring(*gens)
                else:
                    return QQ.poly_ring(*gens)

            r = re.match(cls._re_fraction, domain)

            if r is not None:
                ground, gens = r.groups()

                gens = map(sympify, gens.split(','))

                if ground in ['Z', 'ZZ']:
                    return ZZ.frac_field(*gens)
                else:
                    return QQ.frac_field(*gens)

            r = re.match(cls._re_algebraic, domain)

            if r is not None:
                gens = map(sympify, r.groups()[1].split(','))
                return QQ.algebraic_field(*gens)

        raise OptionError('expected a valid domain specification, got %s' % domain)

    @classmethod
    def postprocess(cls, options):
        if options['domain'].is_Composite and set(options['domain'].gens) & set(options['gens']):
            raise PolynomialError("ground domain and generators interferes together")

class Split(Option):
    __metaclass__ = OptionType

    option = 'split'

    requires = []
    excludes = ['field', 'greedy', 'domain', 'gaussian', 'extension', 'modulus', 'symmetric']

    default = None

    @classmethod
    def preprocess(cls, split):
        if isinstance(split, bool):
            return split
        else:
            raise OptionError("invalid argument for 'split' option")

    @classmethod
    def postprocess(cls, options):
        raise NotImplementedError("'split' option is not implemented yet")

class Gaussian(BooleanOption):
    __metaclass__ = OptionType

    option = 'gaussian'

    requires = []
    excludes = ['field', 'greedy', 'domain', 'split', 'extension', 'modulus', 'symmetric']

    default = None

    @classmethod
    def postprocess(cls, options):
        if options['gaussian'] is True:
            options['extension'] = set([S.ImaginaryUnit])

class Extension(Option):
    __metaclass__ = OptionType

    option = 'extension'

    requires = []
    excludes = ['field', 'greedy', 'domain', 'split', 'gaussian', 'modulus', 'symmetric']

    default = None

    @classmethod
    def preprocess(cls, extension):
        if isinstance(extension, bool):
            if extension is False:
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
    __metaclass__ = OptionType

    option = 'modulus'

    requires = []
    excludes = ['field', 'greedy', 'split', 'gaussian', 'extension']

    default = None

    @classmethod
    def preprocess(cls, modulus):
        modulus = sympify(modulus)

        if modulus.is_Integer and isprime(modulus):
            return int(modulus)
        else:
            raise OptionError("'modulus' must be a prime integer, got %s" % modulus)

    @classmethod
    def postprocess(cls, options):
        if 'domain' in options:
            if not options['domain'].is_ZZ:
                raise PolynomialError("'domain' must have ZZ category when give with 'modulus'")
        else:
            options['domain'] = ZZ

class Symmetric(BooleanOption):
    __metaclass__ = OptionType

    option = 'symmetric'

    requires = ['modulus']
    excludes = ['field', 'greedy', 'split', 'gaussian', 'extension']

    default = None

class Strict(BooleanOption):
    __metaclass__ = OptionType

    option = 'strict'
    default = True

class Auto(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'auto'
    default = True

class Frac(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'frac'
    default = False

class Formal(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'formal'
    default = False

class Polys(BooleanOption, Flag):
    """ """

    __metaclass__ = OptionType

    option = 'polys'
    default = None

class Gen(Option, Flag):
    __metaclass__ = OptionType

    option = 'gen'

    requires = []
    excludes = []

    default = None

    @classmethod
    def preprocess(cls, gen):
        if isinstance(gen, Basic):
            return gen
        else:
            raise OptionError("invalid argument for 'gen' option")

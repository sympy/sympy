""" Caching facility for SymPy """

# TODO: refactor CACHE & friends into class?

# global cache registry:
CACHE = []  # [] of
            #    (item, {} or tuple of {})

def print_cache():
    """print cache content"""

    for item, cache in CACHE:
        item = str(item)
        head = '='*len(item)

        print head
        print item
        print head

        if not isinstance(cache, tuple):
            cache = (cache,)
            shown = False
        else:
            shown = True

        for i, kv in enumerate(cache):
            if shown:
                print '\n*** %i ***\n' % i

            for k, v in kv.iteritems():
                print '  %s :\t%s' % (k, v)

def clear_cache():
    """clear cache content"""
    for item, cache in CACHE:
        if not isinstance(cache, tuple):
            cache = (cache,)

        for kv in cache:
            kv.clear()

########################################

def __cacheit_nocache(func):
    return func

def __cacheit(func):
    """caching decorator.

       important: the result of cached function must be *immutable*


       Example
       -------

       @cacheit
       def f(a,b):
           return a+b


       @cacheit
       def f(a,b):
           return [a,b] # <-- WRONG, returns mutable object


       to force cacheit to check returned results mutability and consistency,
       set environment variable SYMPY_USE_CACHE to 'debug'
    """

    func._cache_it_cache = func_cache_it_cache = {}
    CACHE.append((func, func_cache_it_cache))

    def wrapper(*args, **kw_args):
        if kw_args:
            keys = kw_args.keys()
            keys.sort()
            items = [(k+'=',kw_args[k]) for k in keys]
            k = args + tuple(items)
        else:
            k = args
        try:
            return func_cache_it_cache[k]
        except KeyError:
            pass
        func_cache_it_cache[k] = r = func(*args, **kw_args)
        return r

    wrapper.__doc__ = func.__doc__
    wrapper.__name__ = func.__name__

    return wrapper

def __cacheit_debug(func):
    """cacheit + code to check cache consistency"""
    cfunc = __cacheit(func)

    def wrapper(*args, **kw_args):
        # always call function itself and compare it with cached version
        r1 = func (*args, **kw_args)
        r2 = cfunc(*args, **kw_args)

        # try to see if the result is immutable
        #
        # this works because:
        #
        # hash([1,2,3])         -> raise TypeError
        # hash({'a':1, 'b':2})  -> raise TypeError
        # hash((1,[2,3]))       -> raise TypeError
        #
        # hash((1,2,3))         -> just computes the hash
        hash(r1), hash(r2)

        # also see if returned values are the same
        assert r1 == r2

        return r1

    wrapper.__doc__ = func.__doc__
    wrapper.__name__ = func.__name__

    return wrapper

def __cacheit_nondummy(func):
    func._cache_it_cache = func_cache_it_cache = {}
    CACHE.append((func, func_cache_it_cache))

    def wrapper(*args, **kw_args):
        if kw_args:
            try:
                dummy = kw_args['dummy']
            except KeyError:
                dummy = None
            if dummy:
                return func(*args, **kw_args)
            keys = kw_args.keys()
            keys.sort()
            items = [(k+'=',kw_args[k]) for k in keys]
            k = args + tuple(items)
        else:
            k = args
        try:
            return func_cache_it_cache[k]
        except KeyError:
            pass
        func_cache_it_cache[k] = r = func(*args, **kw_args)
        return r

    wrapper.__doc__ = func.__doc__
    wrapper.__name__ = func.__name__

    return wrapper

class MemoizerArg:
    """ See Memoizer.
    """

    def __init__(self, allowed_types, converter = None, name = None):
        self._allowed_types = allowed_types
        self.converter = converter
        self.name = name

    def fix_allowed_types(self, have_been_here={}):
        from basic import C
        i = id(self)
        if have_been_here.get(i): return
        allowed_types = self._allowed_types
        if isinstance(allowed_types, str):
            self.allowed_types = getattr(C, allowed_types)
        elif isinstance(allowed_types, (tuple, list)):
            new_allowed_types = []
            for t in allowed_types:
                if isinstance(t, str):
                    t = getattr(C, t)
                new_allowed_types.append(t)
            self.allowed_types = tuple(new_allowed_types)
        else:
            self.allowed_types = allowed_types
        have_been_here[i] = True
        return

    def process(self, obj, func, index = None):
        if isinstance(obj, self.allowed_types):
            if self.converter is not None:
                obj = self.converter(obj)
            return obj
        func_src = '%s:%s:function %s' % (func.func_code.co_filename, func.func_code.co_firstlineno, func.func_name)
        if index is None:
            raise ValueError('%s return value must be of type %r but got %r' % (func_src, self.allowed_types, obj))
        if isinstance(index, (int,long)):
            raise ValueError('%s %s-th argument must be of type %r but got %r' % (func_src, index, self.allowed_types, obj))
        if isinstance(index, str):
            raise ValueError('%s %r keyword argument must be of type %r but got %r' % (func_src, index, self.allowed_types, obj))
        raise NotImplementedError(`index,type(index)`)

class Memoizer:
    """ Memoizer function decorator generator.

    Features:
      - checks that function arguments have allowed types
      - optionally apply converters to arguments
      - cache the results of function calls
      - optionally apply converter to function values

    Usage:

      @Memoizer(<allowed types for argument 0>,
                MemoizerArg(<allowed types for argument 1>),
                MemoizerArg(<allowed types for argument 2>, <convert argument before function call>),
                MemoizerArg(<allowed types for argument 3>, <convert argument before function call>, name=<kw argument name>),
                ...
                return_value_converter = <None or converter function, usually makes a copy>
                )
      def function(<arguments>, <kw_arguments>):
          ...

    Details:
      - if allowed type is string object then there `C` must have attribute
        with the string name that is used as the allowed type --- this is needed
        for applying Memoizer decorator to Basic methods when Basic definition
        is not defined.

    Restrictions:
      - arguments must be immutable
      - when function values are mutable then one must use return_value_converter to
        deep copy the returned values

    Ref: http://en.wikipedia.org/wiki/Memoization
    """

    def __init__(self, *arg_templates, **kw_arg_templates):
        new_arg_templates = []
        for t in arg_templates:
            if not isinstance(t, MemoizerArg):
                t = MemoizerArg(t)
            new_arg_templates.append(t)
        self.arg_templates = tuple(new_arg_templates)
        return_value_converter = kw_arg_templates.pop('return_value_converter', None)
        self.kw_arg_templates = kw_arg_templates.copy()
        for template in self.arg_templates:
            if template.name is not None:
                self.kw_arg_templates[template.name] = template
        if return_value_converter is None:
            self.return_value_converter = lambda obj: obj
        else:
            self.return_value_converter = return_value_converter

    def fix_allowed_types(self, have_been_here={}):
        i = id(self)
        if have_been_here.get(i): return
        for t in self.arg_templates:
            t.fix_allowed_types()
        for k,t in self.kw_arg_templates.items():
            t.fix_allowed_types()
        have_been_here[i] = True

    def __call__(self, func):
        cache = {}
        value_cache = {}
        CACHE.append((func, (cache, value_cache)))

        def wrapper(*args, **kw_args):
            kw_items = tuple(kw_args.items())
            try:
                return self.return_value_converter(cache[args,kw_items])
            except KeyError:
                pass
            self.fix_allowed_types()
            new_args = tuple([template.process(a,func,i) for (a, template, i) in zip(args, self.arg_templates, range(len(args)))])
            assert len(args)==len(new_args)
            new_kw_args = {}
            for k, v in kw_items:
                template = self.kw_arg_templates[k]
                v = template.process(v, func, k)
                new_kw_args[k] = v
            new_kw_items = tuple(new_kw_args.items())
            try:
                return self.return_value_converter(cache[new_args, new_kw_items])
            except KeyError:
                r = func(*new_args, **new_kw_args)
                try:
                    try:
                        r = value_cache[r]
                    except KeyError:
                        value_cache[r] = r
                except TypeError:
                    pass
                cache[new_args, new_kw_items] = cache[args, kw_items] = r
                return self.return_value_converter(r)
        return wrapper


class Memoizer_nocache(Memoizer):

    def __call__(self, func):
        # XXX I would be happy just to return func, but we need to provide
        # argument convertion, and it is really needed for e.g. Real("0.5")
        def wrapper(*args, **kw_args):
            kw_items = tuple(kw_args.items())
            self.fix_allowed_types()
            new_args = tuple([template.process(a,func,i) for (a, template, i) in zip(args, self.arg_templates, range(len(args)))])
            assert len(args)==len(new_args)
            new_kw_args = {}
            for k, v in kw_items:
                template = self.kw_arg_templates[k]
                v = template.process(v, func, k)
                new_kw_args[k] = v

            r = func(*new_args, **new_kw_args)
            return self.return_value_converter(r)

        return wrapper



# SYMPY_USE_CACHE=yes/no/debug
def __usecache():
    import os
    return os.getenv('SYMPY_USE_CACHE', 'yes').lower()
usecache = __usecache()

if usecache=='no':
    Memoizer            = Memoizer_nocache
    cacheit             = __cacheit_nocache
elif usecache=='yes':
    cacheit = __cacheit
elif usecache=='debug':
    cacheit = __cacheit_debug   # a lot slower
else:
    raise RuntimeError('unknown argument in SYMPY_USE_CACHE: %s' % usecache)

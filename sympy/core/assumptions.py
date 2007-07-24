
class AssumeMeths(object):
    """ Define default assumption methods.
    
    AssumeMeths should be used to derive Basic class only.

    All symbolic objects have assumption attributes that
    can be accessed via .is_<assumption name> attribute.
    Assumptions determine certain properties of symbolic
    objects. Assumptions can have 3 possible values: True, False, None.
    None is returned when it is impossible to say something
    about the property. For example, ImaginaryUnit() is
    not positive neither negative. By default, all symbolic
    values are in the largest set in the given context without
    specifying the property. For example, a symbol that
    has a property being integer, is also real, complex, etc.

    Here follows a list of possible assumption names:

        - commutative   - object commutes with any other object with
                          respect to multiplication operation.
        - real          - object can have only values from the set
                          of real numbers
        - integer       - object can have only values from the set
                          of integers
        - bounded       - object absolute value is bounded
        - dummy         - used for marking dummy symbols
        - positive      - object can have only positive values
        - negative      - object can have only negative values
        - nonpositive      - object can have only nonpositive values
        - nonnegative      - object can have only nonnegative values
        - comparable    - object.evalf() returns Number object.
        - irrational    - object value cannot be represented exactly by Rational
        - unbounded     - object value is arbitrarily large
        - infinitesimal - object value is infinitesimal
        - order         - expression is not contained in Order(order).

    Example rules:

      positive=True|False -> nonpositive=not positive, real=True
      positive=None -> negative=None
      unbounded=False|True -> bounded=not unbounded
      irrational=True -> real=True

    Exceptions:
      positive=negative=False for Zero() instance

    Implementation note: assumption values are stored in
    ._assumption dictionary or are returned by getter methods (with
    property decorators) or are attributes of objects/classes.

    Examples:
    
        - True, when we are sure about a property. For example, when we are
        working only with real numbers:
        >>> from sympy import *
        >>> Symbol('x', real = True)
        x
        
        - False
        
        - None (if you don't know if the property is True or false)

    

    """

    _assume_aliases = {} # aliases, "a key means values"
    _assume_aliases['nni'] = ('integer','nonnegative')
    _assume_aliases['npi'] = ('integer','nonpositive')
    _assume_aliases['pi'] = ('integer','positive')
    _assume_aliases['ni'] = ('integer','negative')

    _properties = ['dummy','order'] # todo: rm is_order

    # inclusion relations (subset, superset)
    _assume_rels = (('prime', 'integer'),
                    ('odd', 'integer'),
                    ('integer', 'rational'),
                    ('rational', 'real'),
                    ('real', 'complex'),
                    ('positive', 'real'),
                    ('negative', 'real'),
                    ('finite', 'bounded'),
                    )

    # implications (property -> super property)
    _assume_impl = (('zero','infinitesimal'),
                    #('negative', 'nonpositive'),
                    #('positive', 'nonnegative'),
                    ('finite', 'nonzero'),
                    ('zero', 'even'),
                    ('complex', 'commutative'),
                    )

    # (property, negative property) mapping
    _assume_negs = {'bounded': 'unbounded',
                    'commutative': 'noncommutative',
                    'complex': 'noncomplex',
                    'finite': 'infinitesimal',
                    'integer': 'noninteger',
                    'negative': 'nonnegative',
                    'odd': 'even',
                    'positive': 'nonpositive',
                    'prime': 'composite',
                    'rational': 'irrational',
                    'real': 'imaginary',
                    'zero': 'nonzero',
                    'homogeneous':'inhomogeneous'}

    _assume_inegs = {}
    for k,v in _assume_negs.items(): _assume_inegs[v] = k

    _assume_defined = ('integer','rational','real','complex','noninteger','irrational',
                       'imaginary','noncomplex',
                       'even','odd','prime','composite','zero','nonzero',
                       'negative','nonnegative','positive','nonpositive',
                       'finite','infinitesimal','bounded','unbounded',
                       'commutative','noncommutative',
                       'homogeneous','inhomogeneous',
                       'comparable',
                       'dummy','order',
                       'nni','pi',
                       'evaluate')

    def _change_assumption(self, d, name, value, extra_msg = ''):
        default_assumptions = self.__class__.default_assumptions
        fixedvalue = default_assumptions.get(name, None)
        if value is None:
            oldvalue = d.pop(name, None)
        else:
            oldvalue = d.get(name, fixedvalue)
            if oldvalue is not None and  oldvalue != value and fixedvalue is not None:
                raise TypeError('%s: cannot change fixed assumption item from %s=%s to %s=%s%s'\
                                % (self.__class__.__name__, name, oldvalue, name, value, extra_msg))
            d[name] = value

    def assume(self, **assumptions):
        """ Modify object assumptions in-situ.

        Usage examples:
          obj.assume(commutative = True,  # obj is in commutative center
                     real = None          # assumption that obj is real will be removed
                     )
          obj.is_commutative              # check if object is commutative

        User is responsible for setting reasonable assumptions.
        """
        default_assumptions = self.__class__.default_assumptions
        for k,v in default_assumptions.items():
            if assumptions.has_key(k):
                nv = assumptions[k]
                if nv!=v:
                    raise TypeError("%s: assumption %r is fixed to %r for this class." \
                                    % (self.__class__.__name__,k,v))
            assumptions[k] = v
        d = self._assumptions = getattr(self, '_assumptions', {})

        ###
        if "negative" in assumptions:
            if "positive" not in assumptions:
                self._change_assumption(d, "positive", not assumptions["negative"])
        elif "positive" in assumptions:
            self._change_assumption(d, "negative", not assumptions["positive"])
        ###

        processed = {}
        aliases = self._assume_aliases 
        negs = self._assume_negs
        inegs = self._assume_inegs
        rels = self._assume_rels
        impl = self._assume_impl
        defined = self._assume_defined
        while assumptions:
            k, v = assumptions.popitem()
            # obsolete, to be removed:
            if k.startswith('is_'):
                k = k[3:]
                assumptions[k] = v
                continue
            #
            if aliases.has_key(k):
                for a in aliases[k]:
                    if a not in processed:
                        assumptions[a] = v
                    else:
                        old_v = processed[a]
                        assert old_v==v,`k,a,old_v,v`
                processed[k] = v
                continue

            if k not in defined:
                raise ValueError('unrecognized assumption item (%r:%r)' % (k,v))

            if k=='order':
                v = self.Order(v)
            elif v is not None: v = bool(v)

            if inegs.has_key(k):
                a = inegs[k]
                if v is None:
                    if a not in processed:
                        if a not in assumptions:
                            assumptions[a] = v
                        elif assumptions[a] != v:
                            raise ValueError('%s: detected inconsistency between %s=%s and %s=%s in negation' \
                                            % (self.__class__,k, v, a, assumptions[a]))
                    elif processed[a] != v:
                        raise ValueError('%s:detected inconsistency between %s=%s and %s=%s in processed negation'\
                                        % (self.__class__,k, v, a, processed[a]))
                else:
                    if a not in processed:
                        if a not in assumptions:
                            assumptions[a] = not v
                        elif assumptions[a] != (not v):
                            raise ValueError('%s: detected inconsistency between %s=%s and %s=%s in negation' \
                                            % (self.__class__,k, v, a, assumptions[a]))
                    elif processed[a] != (not v):
                        raise ValueError('%s: detected inconsistency between %s=%s and %s=%s in processed negation' \
                                        % (self.__class__,k, v, a, processed[a]))
                processed[k] = v
                continue

            self._change_assumption(d, k, v)
            if negs.has_key(k):
                a = negs[k]
                if v is None:
                    self._change_assumption(d, a, v)
                else:
                    self._change_assumption(d, a, not v)

            assert not processed.has_key(k),`processed, k, a, v`
            processed[k] = v

            for p1,p2 in rels:
                if k==p1:
                    # k is a subset of p2
                    if v is not None:
                        if p2 in processed:
                            assert processed[p2]==True, `k,v,p2,processed[p2]`
                        else:
                            assumptions[p2] = True
                if k==p2:
                    # k contains p1
                    if not v:
                        if p1 in processed:
                            assert processed[p1]==None, `self.__class__, k,v,p1,processed[p1]`
                        else:
                            assumptions[p1] = None
            for p1,p2 in impl:
                if k==p1:
                    if v:
                        if p2 in processed:
                            assert processed[p2]==True, `k,v,p2,processed[p2]`
                        else:
                            assumptions[p2] = True
            
        return

    def _assume_hashable_content(self):
        d = self._assumptions
        keys = d.keys()
        keys.sort()
        return tuple([(k+'=', d[k]) for k in keys])

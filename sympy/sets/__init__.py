from .sets import (Set, Interval, Union, EmptySet, FiniteSet, ProductSet,
        Intersection, imageset, Complement, SymmetricDifference)
from .fancysets import ImageSet, Range, ComplexRegion
from .contains import Contains
from .conditionset import ConditionSet
from ..core.singleton import S

Reals = S.Reals
del S

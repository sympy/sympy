from .sets import (Set, Interval, Union, EmptySet, FiniteSet, ProductSet,
        Intersection, imageset, Complement, SymmetricDifference)
from .fancysets import ImageSet, Range, ComplexRegion
from .contains import Contains
from .conditionset import ConditionSet
from .ordinals import Ordinal, OmegaPower, ord0
from ..core.singleton import S
Reals = S.Reals
del S

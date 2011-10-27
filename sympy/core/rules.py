"""
Replacement rules.
"""

class Transform(object):
    """
    Generic transformation rule.
    """
    def __init__(self, transform, filter=lambda x: True):
        self.transform = transform
        self.filter = filter

    def __contains__(self, item):
        return self.filter(item)

    def __getitem__(self, key):
        return self.transform(key)

    def get(self, item, default=None):
        if item in self:
            return self[item]
        else:
            return default

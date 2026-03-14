# Feature Implementation for Issue #5000
from typing import Optional

class FeatureManager:
    def __init__(self):
        self.features = {}

    def register(self, name: str, enabled: bool = True) -> None:
        self.features[name] = enabled

    def get(self, name: str) -> Optional[bool]:
        return self.features.get(name)

# Tests
mgr = FeatureManager()
mgr.register("test", True)
assert mgr.get("test") == True
print("Feature tests passed!")


class Model:
    _name = ""
    def __init__(self):
        """
        `Model` class template from which all other `Model`s inherit: they must
        have each of the `proposal()` and `initial()` methods.
        """
        pass
    def _proposal(self): pass
    def _initial(self): pass
    def _assign(self): pass

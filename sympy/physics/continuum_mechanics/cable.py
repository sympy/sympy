"""
This module can be used to solve problems related
to 2D Cables.
"""

from sympy.core.symbol import Symbol
from sympy.core.sympify import sympify

class Cable:
    """
    Cables are structures in engineering that support
    the applied transverse loads through the tesnsile
    resistance developed in its members.
    
    Cables are widely used in suspension bridges, tension 
    leg offshore platforms, transmission lines, and find 
    use in several other engineering applications.
    """
    def __init__(self):
        """
        Initializes the class.
        """
        self._left_support = []
        self._right_support = []
        self._supports = {}
        self._support_labels = []
        self._loads = {}
        self._loads_position = {}
        self._length = 0
        self._reaction_loads = {}

    @property
    def supports(self):
        """
        Returns the supports of the cable along with their 
        positions.
        """
        return self._supports

    @property
    def left_support(self):
        """
        Returns the position of the left support.
        """
        return self._left_support
    
    @property
    def right_support(self):
        """
        Returns the position of the right support.
        """
        return self._right_support
    
    @property
    def loads(self):
        """
        Returns the magnitude and direction of the loads
        acting on the cable.
        """
        return self._loads
    
    @property
    def loads_position(self):
        """
        Returns the position of the loads acting on the 
        cable. 
        """
        return self._loads_position

    @property
    def length(self):
        """
        Returns the length of the cable.
        """
        return self._length
    
    @property
    def reaction_loads(self):
        """
        Returns the reaction forces at the supports.
        """
        return self._reaction_loads

    def add_supports(self, *args):
        """
        This method adds the right and left supports to the Cable object.
        
        Parameters
        ==========
        Pass two tuples of the form (label, x, y) as inputs, where the first
        tuple represents the left support, while the later represents the
        right support. Each element in these tuples are:
        
        label: String or a symbol
            The label for the support.
            
        x: Sympifyable
            The x-coordinate of the position of the support.
            
        y: Sympifyable
            The y-coordinate of the position of the support.
            
        Examples
        ========
        
        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable()
        >>> c.add_supports(('A', 3, 4), ('B', 4, 5))
        >>> c.supports
        {'A': [3, 4], 'B': [4, 5]}
        >>> c.left_support
        [3, 4]
        >>> c.right_support
        [4, 5]
        """
        if len(self._supports) != 0:
            raise ValueError("Supports already exist")

        elif len(args) != 2:
            raise ValueError("Pass only two arguments: left support and right support")

        elif args[0][0] == args[1][0]:
            raise ValueError("Supports can not have the same label")

        elif args[0][1] == args[1][1] and args[0][2] == args[1][2]:
            raise ValueError("Supports can not be at the same position")

        elif args[0][1] > args[1][1]:
            raise ValueError("The x coordinate of left support should be less than its right support counterpart") 
        
        else:
            for i in args:
                label = i[0]
                x = sympify(i[1])
                y = sympify(i[2])
            
                self._supports[label] = [x, y]

                if len(self._left_support) == 0:
                    self._left_support.append(x)
                    self._left_support.append(y)

                else:
                    self._right_support.append(x)
                    self._right_support.append(y)

                self._support_labels.append(label)
            
    def remove_supports(self):
        """
        This method removes the left and right supports, along with 
        any existing load(s).
        
        Examples
        ========
        
        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable()
        >>> c.add_supports(('A', 3, 4), ('B', 4, 5))
        >>> c.supports
        {'A': [3, 4], 'B': [4, 5]}
        >>> c.remove_supports()
        >>> c.supports
        {}
        """
        if len(self._supports) == 0:
            raise ValueError("There are no supports to remove")
        
        else:
            self._supports.clear()
            self._left_support.clear()
            self._right_support.clear()
            self._loads.clear()
            self._loads_position.clear()
            self._length = 0
    
    def change_support(self, label, new_support):
        """
        This method changes the mentioned support with a new support.
        
        Parameters
        ==========
        label: String or symbol
            The label of the support to be changed
        
        new_support: Tuple of the form (new_label, x, y)
            new_label: String or symbol
                The label of the new support
            
            x: Sympifyable
                The x-coordinate of the position of the new support.
                
            y: Sympifyable
                The y-coordinate of the position of the new support.
        
        Examples
        ========
        
        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable()
        >>> c.add_supports(('A', 3, 4), ('B', 4, 5))
        >>> c.supports
        {'A': [3, 4], 'B': [4, 5]}
        >>> c.change_support('B',('C',5,6))
        >>> c.supports
        {'A': [3, 4], 'C': [5, 6]}
        """
        if label not in self._supports:
            raise ValueError("No support exists with the given label")
        
        else:
            if label == self._support_labels[0]:
                if new_support[1] > self._right_support[0]:
                    raise ValueError("x coordinate of left support should be less than its right support counterpart")
                else:
                    self._supports.pop(label)
                    self._left_support.clear()
                    self._support_labels.remove(label)

                    x = sympify(new_support[1])
                    y = sympify(new_support[2])
                    self._left_support = [x, y]
                    self._supports[new_support[0]] = [x, y]
                    self._support_labels.insert(0, new_support[0])
            
            else:
                if new_support[1] < self._left_support[0]:
                    raise ValueError("x coordinate of right support should be greater than its left support counterpart")
                
                else:
                    self._supports.pop(label)
                    self._right_support.clear()
                    self._support_labels.remove(label)

                    x = sympify(new_support[1])
                    y = sympify(new_support[2])
                    self._right_support = [x, y]
                    self._supports[new_support[0]] = [x, y]
                    self._support_labels.append(new_support[0])

    def add_loads(self, *args):
        """
        This method adds load(s) to the cable.
        
        Parameters
        ==========
        This method takes input as tuple(s) of the form
        (label, x, y, magnitude, direction).
        
        label: String or symbol
            The label of the load
        
        x: Sympifyable
            The x-coordinate of the position of the load.
        
        y: Sympifyable
            The y-coordinate of the position of the load.
        
        magnitude: Sympifyable
            Magnitude of the load applied. It must always be positive and any changes in
            the direction of the load are not reflected here.
        
        direction: Sympifyable
            The angle, in degrees, that the load vector makes with the horizontal
            in the counter-clockwise direction. It takes the values 0 to 360,
            inclusive.
            
        Examples
        ========
        
        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable()
        >>> c.add_supports(('A', 3, 4), ('B', 6, 5))
        >>> c.add_loads(('P', 4, 2, 6, 30), ('Q', 5, 3, 12, 0))
        >>> c.loads
        {'P': [6, 30], 'Q': [12, 0]}
        >>> c.loads_position
        {'P': [4, 2], 'Q': [5, 3]}
        """
        if len(self._supports) == 0:
            raise ValueError("Add supports before adding any load(s)")
        
        else:
            for i in args:
                label = i[0]
                x = sympify(i[1])
                y = sympify(i[2])
                magnitude = sympify(i[3])
                direction = sympify(i[4])

                if label in self._loads:
                    raise ValueError("Error adding load " + label + ": label already exists")
                
                else:
                    if x < self._left_support[0] or x > self._right_support[0]:
                        raise ValueError("Error adding load " + label + ": load should be between the supports")
                    
                    else:
                        self._loads[label] = [magnitude, direction]
                        self._loads_position[label] = [x, y]
    
    def remove_loads(self, *args):
        """
        This methods removes the specified loads.
        
        Parameters
        ==========
        This input takes multiple label(s) as input
        label(s): String or symbol
            The label(s) of the loads to be removed.
            
        Examples
        ========
        
        >>> from sympy.physics.continuum_mechanics.cable import Cable
        >>> c = Cable()
        >>> c.add_supports(('A', 3, 4), ('B', 6, 5))
        >>> c.add_loads(('P', 4, 2, 6, 30), ('Q', 5, 3, 12, 0), ('R', 5, 4, 10, 60))
        >>> c.loads
        {'P': [6, 30], 'Q': [12, 0], 'R': [10, 60]}
        >>> c.loads_position
        {'P': [4, 2], 'Q': [5, 3], 'R': [5, 4]}
        >>> c.remove_loads('P', 'R')
        >>> c.loads
        {'Q': [12, 0]}
        >>> c.loads_position
        {'Q': [5, 3]}
        """
        for i in args:
            if i not in self._loads:
                raise ValueError("Error removing load " + i + ": no such load exists")
            
            else:
                self._loads.pop(i)
                self._loads_position.pop(i)
from types import MethodType

class Parser():
    groups: = {} # type: Dict[str,Tuple(str)]

    def __init__(self,*args,**kwargs):
        self.namestack=[]
        for cls in self.__class__.__mro__[::-1]:
            for parser in cls.__dict__:
                if "do" == parser[:2]:
                    setattr(self,parser,self.decorate_parser(getattr(self,parser)))
        super().__init__(*args,**kwargs)

    def decorate_parser(self,parser):
        def run_parser(self,*args,**kwargs):
            name=parser.__name__[2:]

            #setup parser attributes
            self.setup_parser(name)

            #update the stack so that we can nest different parsers
            #(but not the same)
            self.namestack.append(name)
            #perform parsing
            parsed=parser(*args,**kwargs)
            #pop the stack to indicate we are done with parsing
            self.namestack.pop()
            return parsed
        run_parser.__doc__=parser.__doc__
        decorated=MethodType(run_parser,self)
        return decorated

    def setup_parser(self,name):
        #check if this parser is already being used
        #if name in self.namestack:
        #    error='You are not allowed to run {0} inside itself'
        #    raise RecursionError(error.format(name))
        level="_"+name+"_level"
        if not hasattr(self,level):
            setattr(self,"_"+name+"_level",0)
        #setup the name of the method to use for parsing of a class
        #that is being parsed has it.
        method=name+"method"
        if not hasattr(self,method):
            setattr(self,method,"_"+method)
        #setup the method name for when no parser is found
        empty="_"+name+"empty"
        if not hasattr(self,empty):
            setattr(self,empty,"emptyParser")


    def __getattr__(self,attr):
        if attr[:2]=='do':
            self.setup_parser(attr[2:])
            self.namestack.append(attr[2:])
            return self.generic_do
        if attr[1:] in self.namestack:
            return self.generic_parser
        #this can happen if Parser._<parse> or Parser._<parse>_<Class>
        #is used directly, but a Parser._<parse> is not explicitly defined:
        #in that case, we replace the first Parser._<parse> with a
        #generic do to set it up properly:
        if len(attr[1:].split('_'))==1:
            self.setup_parser(attr[1:])
            self.namestack.append(attr[1:])
            return self.generic_do
        raise AttributeError(attr)

    def emptyParser(self,obj):
        return obj


    def generic_do(self,*args,**kwargs):
        name=kwargs.get('name',None)
        if name is None:
            name=self.namestack[-1]
        else:
            self.namestack.append(name)
        result= getattr(self,'_'+name)(*args,**kwargs)
        self.namestack.pop()
        return result

    def generic_parser(self,expr,name=None,**kwargs):
        if name is None:
            name=self.namestack[-1]
        """Internal dispatcher
        Tries the following concepts to parse an expression:
            1. Let the object parse itself if it knows how.
            2. Take the best fitting method defined in the parser.
            3. As fall-back use the _parse_Object method for the printer
                    default is to .
        """
        setattr(self,'_'+name+'_level',getattr(self,'_'+name+'_level')+1)
        try:
            # If the printer defines a name for a printing method
            # (Printer.printmethod) and the object knows for itself how it
            # should be printed, use that method.
            parsemethod=getattr(self,name+'method')
            if (parsemethod and hasattr(expr, parsemethod) #
                    and not 'BasicMeta' in expr.__class__.__mro__):
                return getattr(expr, parsemethod)(self, **kwargs) #
            # See if the class of expr is known, or if one of its super
            # classes is known, and use that print function
            # Exception: ignore the subclasses of Undefined, so that, e.g.,
            # Function('gamma') does not get dispatched to _print_gamma
            classes = tuple(mro.__name__ for mro in type(expr).__mro__)
            if 'AppliedUndef' in classes:
                classes = classes[classes.index('AppliedUndef'):]
            if 'UndefinedFunction' in classes:
                classes = classes[classes.index('UndefinedFunction'):]
            # Another exception: if someone subclasses a known function, e.g.,
            # gamma, and changes the name, then ignore _print_gamma
            if 'Function' in classes:
                i = classes.index('Function')
                classes = tuple(c for c in classes[:i] if \
                    c == classes[0] or \
                    c.endswith("Base")) + classes[i:]
            #If a class is inside any group, add the group Class right after
            #the class which is part of the group
            for group_class,group in self.groups.items():
                idx=0
                for i,cls in enumerate(classes):
                    if cls in group:
                        idx=i+1
                if idx:
                    classes=classes[:idx] +(group_class,) + classes[idx:]
            for cls in classes:
                parsemethodname = '_'+name+'_' + cls
                parsemethod = getattr(self, parsemethodname, None)
                if parsemethod is not None:
                    return parsemethod(expr, **kwargs)
            #use the empty parser as a last resort
            empty=getattr(self,getattr(self,"_"+name+"empty"))
            return empty(expr, **kwargs)
        finally:
            setattr(self,'_'+name+'_level',getattr(self,'_'+name+'_level')-1)

@doctest_depends_on(modules=('lxml','StringIO','os',))
def openmath2cmml(omstring,simple=False):
"""
Transforms a string in Openmath to Content MathML (simple or not)
XSL Templates from https://svn.omdoc.org/repos/omdoc/projects/emacsmode/nomdoc-mode/xsl/ although they are
"""
    
    from lxml import etree
    from lxml import objectify
    from StringIO import *
    import os
    
    if simple:
        xslt_tree = etree.XML(open(os.path.join(request.folder,'mathml/data/omtosimcmml.xsl')).read())
    else:
        xslt_tree = etree.XML(open(os.path.join(request.folder,'mathml/data/omtocmml.xsl')).read())
    
    transform = etree.XSLT(xslt_tree)
    omstring= omstring.replace(' xmlns="', ' xmlnamespace="')
    parser = etree.XMLParser(ns_clean=True,remove_pis=True,remove_comments=True)
    tree   = etree.parse(StringIO(omstring), parser)
    objectify.deannotate(tree,cleanup_namespaces=True,xsi=True,xsi_nil=True)
    cmmlstring_tree=transform(tree)
    cmmlstring=etree.tostring(cmmlstring_tree.getroot())
    return(cmmlstring)


@doctest_depends_on(modules=('lxml','StringIO',))    
def parseCMML(mmlinput):
    """
    This parses Content MathML into a Python expression and a set of variables which can be sympified afterwards. At the moment, only basic operators are taking into account.
    It returns the expression and the set of variables inside it
    """
    from lxml import etree
    from StringIO import *
    from lxml import objectify
    mmlinput= mmlinput.replace(' xmlns="', ' xmlnamespace="')
    parser = etree.XMLParser(ns_clean=True,remove_pis=True,remove_comments=True)
    tree   = etree.parse(StringIO(mmlinput), parser)
    objectify.deannotate(tree,cleanup_namespaces=True,xsi=True,xsi_nil=True)
    mmlinput=etree.tostring(tree.getroot())
    exppy="" #this is the python expression
    symvars=[]  #these are symbolic variables which can eventually take part in the expression
    events = ("start", "end")
    level = 0
    context = etree.iterparse(StringIO(mmlinput),events=events)
    for action, elem in context:
        
        if level:
            continue
                
        if (action=='start') and (elem.tag=='apply'):
            exppy+='('
            level += 1
            opelem=elem[0]
            if (opelem.tag=='divide'):
                mmlaux=etree.tostring(opelem.getnext())
                (a,b)=parseCMML(mmlaux)
                symvars.append(b)
                exppy+=a
                exppy+='/'
                mmlaux=etree.tostring(opelem.getnext().getnext())
                (a,b)=parseCMML(mmlaux)
                symvars.append(b)
                exppy+=a
            if (opelem.tag=='power'):
                mmlaux=etree.tostring(opelem.getnext())
                (a,b)=parseCMML(mmlaux)
                symvars.append(b)
                exppy+=a
                exppy+='**'
                mmlaux=etree.tostring(opelem.getnext().getnext())
                (a,b)=parseCMML(mmlaux)
                symvars.append(b)
                exppy+=a
            if (opelem.tag=='plus'):
                sib=opelem.getnext()
                while sib!= None:
                    mmlaux=etree.tostring(sib)
                    (a,b)=parseCMML(mmlaux)
                    symvars.append(b)
                    exppy+=a
                    if sib.getnext()!=None:
                        exppy+='+'
                    sib=sib.getnext()
            if (opelem.tag=='minus'):
                sib=opelem.getnext()
                if sib.getnext()!= None:
                    #binary operator
                    mmlaux=etree.tostring(sib)
                    (a,b)=parseCMML(mmlaux)
                    symvars.append(b)
                    exppy+=a
                    exppy+='-'
                    mmlaux=etree.tostring(sib.getnext())
                    (a,b)=parseCMML(mmlaux)
                    symvars.append(b)
                    exppy+=a
                else:
                    #unary operator
                    exppy+='-'
                    mmlaux=etree.tostring(sib)
                    (a,b)=parseCMML(mmlaux)
                    symvars.append(b)
                    exppy+=a
            if (opelem.tag=='times'):
                sib=opelem.getnext()
                while sib!= None:
                    mmlaux=etree.tostring(sib)
                    (a,b)=parseCMML(mmlaux)
                    symvars.append(b)
                    exppy+=a
                    if sib.getnext()!=None:
                        exppy+='*'
                    sib=sib.getnext()        
            exppy+=')'        
                
        if (action=='end') and (elem.tag=='apply'):
            level -= 1
            
        if action=='start' and elem.tag=='cn': #this is a number
            exppy+=elem.text
        if action=='start' and elem.tag=='ci': #this is a variable
            exppy+=elem.text
            symvars.append(elem.text) #we'll return the variable, so sympy can sympify it afterwards
        
    return (exppy, symvars)


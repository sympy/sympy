from itertools import chain
from sympy.core.symbol import Dummy
from sympy.printing.fcode import FCodePrinter
from sympy.codegen.fnodes import Module


def render_as_module(content, name, declarations=(), printer_settings=None):
    printer_settings = printer_settings or {'standard': 2003, 'source_format': 'free'}
    printer = FCodePrinter(printer_settings)
    dummy = Dummy()
    if isinstance(content, Module):
        raise ValueError("This function expects to construct a module on its own.")
    mod = Module(name, chain(declarations, [dummy]), content)
    fstr = printer.doprint(mod)
    module_use_str = '   %s\n' % '   \n'.join(['use %s, only: %s' % (k, ', '.join(v)) for
                                                k, v in printer.module_uses.items()])
    module_use_str += '   implicit none\n'
    module_use_str += '   private\n'
    module_use_str += '   public %s\n' % ', '.join([str(node.name) for node in content if getattr(node, 'name', None)])
    return fstr.replace(printer.doprint(dummy), module_use_str)
    return fstr

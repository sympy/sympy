import importlib
m = importlib.import_module('sympy.printing.llvmjitcode')
print('Imported module:', m.__name__)
print('Has llvm_callable:', hasattr(m, 'llvm_callable'))
print('Classes:', [name for name in dir(m) if name.endswith('JitPrinter') or name.startswith('LLVMJit') or name=='CodeSignature'])

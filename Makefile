# build or clean Cythonized modules

all:
	python build.py build_ext --inplace

clean:
	rm -rf sympy/polys/*.{pyc,c,so}
	rm -rf build/


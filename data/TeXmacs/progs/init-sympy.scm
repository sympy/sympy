(plugin-configure sympy
  (:require (url-exists-in-path? "tm_sympy"))
  (:launch "tm_sympy --texmacs")
  (:session "SymPy"))

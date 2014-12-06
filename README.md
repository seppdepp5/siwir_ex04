siwir_ex03
==========

Try to run with

    ./cg 100 100 10000 0.01
  
~~I do not know why, but~~ if you do to much iterations the residual is gaining :-O.

Think this lies in the nature of cg ~~but could also be a bug~~ .

> [...] as it produces the exact solution after a finite number of iterations, which is not larger than the size of the matrix [...]

from: http://en.wikipedia.org/wiki/Conjugate_gradient_method#Convergence_properties_of_the_conjugate_gradient_method

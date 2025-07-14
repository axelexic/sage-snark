from sage.matrix.all import Matrix
from sage.rings.integer import Integer
from sage.rings.polynomial.all import PolynomialRing, Polynomial
from sage.modules.all import vector
from sage.rings.rational_field import QQ

try:
    from utils.multivariates import *
except:
    print(f"`sage-snark` directory is not on current PYTHONPATH.")


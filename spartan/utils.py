from sage.rings.finite_rings.all import FiniteField
from sage.rings.polynomial.all import Polynomial, PolynomialRing
from sage.rings.integer import Integer
from sage.modules.all import vector
from sage.matrix.all import Matrix
from sage.misc.all import prod
from typing import Any, List, Dict

def bit_decomp_le(value : int, width : int = 0) -> List[Integer]:
    """
    Do bit decomposition of integer `value` and pad it with zeros upto `width` length. Return the result in little-endian format.
    """
    bits = Integer(value).bits()
    bit_len = len(bits)
    if bit_len < width:
        bits.extend([0]* (width - bit_len))

    return bits

def bit_decomp_be(value : int, width : int = 0) -> List[Integer]:
    """
    Do bit decomposition of integer `value` and pad it with zeros upto `width` length. Return the result in big-endian format.
    """
    return bit_decomp_le(value, width)[::-1]

def bit_decomp_dict(value : int, variables : List[Any]) -> Dict[Any, Integer]:
    """
    Do bit decomposition of integer `value` and return it as a key-value pair
    where keys are assumed to represent bits in little endian format.

    For example, if value = b0 b1 b2 b3 where b_i is bits in little-endian
    format, and variables are variables = [X, Y, Z, W], then this method returns
    {
        X : b0,
        Y : b1,
        Z : b2,
        W : b3
    }
    """
    var_len = len(variables)
    bits = bit_decomp_le(value, var_len)
    return dict(zip(variables, bits))


def eq_poly(x_vals : List[Any], y_vals: List [Any]) -> Any:
    """
    Given `x_vals` and `y_vals`, compute the eq polynomial
        Π(x_i*y_i + (1-x_i)*(1-y_i))
    """

    assert len(x_vals) == len(y_vals)
    eq_one = lambda x, y: x*y + (1 - x)*(1-y)
    return prod( map(lambda xy : eq_one(xy[0], xy[1]), zip(x_vals, y_vals)) )

def eq_bit_decomp(PolyRing : PolynomialRing, index : int, poly_gens : List[Any] = None) -> Polynomial:
    """
    Computes the bit decomposition of index making sure that number of variables
    in PolyRing matches the maximum possible index value (i.e., if the number of
    variable is n, then the index must be in the range 0 and 2^n-1 (inclusive).
    """
    gens = poly_gens or list(PolyRing.gens())
    var_len = len(gens)
    max_index = 2**var_len - 1

    if index < 0 or index > max_index:
        raise ValueError(f"Invalid index {index}. Cannot be represented in polynomial ring {PolyRing}")
    bits = bit_decomp_le(index, var_len)
    return eq_poly(gens, bits)

def hadamard_product(vec1, vec2):
    assert len(vec1) == len(vec2)
    return [a*b for (a,b) in zip(vec1, vec2)]

def hypercube_sum(poly : Polynomial, skip_variables = [], hypercube = [0,1]):
    sc_val = poly
    variables = poly.variables()
    skip_variables = skip_variables or []

    for xvar in variables:
        if xvar in skip_variables:
            continue

        this_round = 0

        for hval in hypercube:
            this_round += sc_val.subs({xvar : hval})

        sc_val = this_round

    return sc_val

if __name__ == '__main__':
    from sage.rings.finite_rings.all import GF

    Fp = GF(2**31 - 1)
    R = PolynomialRing(Fp, ["X", "Y", "Z", "W"])
    gens = R.gens()
    max_index = 2**len(gens)

    for i in range(max_index):
        p = eq_bit_decomp(R, i)

        for j in range(max_index):
            decompo = bit_decomp_dict(j, gens)
            expected = p.subs(decompo)

            if i == j:
                assert expected == 1
            else:
                assert expected == 0
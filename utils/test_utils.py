from sage.rings.finite_rings.all import FiniteField, GF
from sage.modules.all import vector, random_vector
from sage.misc.functional import log
from sage.misc.prandom import random
from sage.matrix.all import Matrix, random_matrix
from typing import List

try:
    import multivariates as utils
except:
    from . import multivariates as utils


class R1CS:
    @classmethod
    def random_element(
        klass,
        Fq : FiniteField,
        num_rows : int,
        num_columns: int = None,
        force_witness : vector = None,
        sparsity : float = None):
        """
        Generate a random R1CS instance over Finite field Fq, with num_row
        constraints. If force_witness is not None, it will create R1CS instances
        that satisfied that witness. If force_witness is None, it will generate
        a random witness with `num_coulmn` entries. `sparsity` specifies show
        Sparse the matrix supposed to be.
        """

        if (num_columns, force_witness) == (None, None):
            raise ValueError(f"Both num_columns and force_witness cannot be None.")

        num_columns = len(force_witness) if force_witness is not None else num_columns

        sparsity  = sparsity or float(log(num_columns, 2))/float(num_columns)

        sparse_row = lambda count : [ Fq.random_element() if random() < sparsity else 0 for _ in range(count) ]

        witness = force_witness or  random_vector(Fq, num_columns)
        A =  Matrix(Fq, [sparse_row(num_columns) for _ in range(num_rows) ])
        B = Matrix(Fq, [sparse_row(num_columns) for _ in range(num_rows) ])

        if witness == vector(Fq, [0 for _ in range(num_columns)]):
            C = random_matrix(Fq, num_rows, num_columns)
            return klass(A, B, C, witness)

        expected_product = utils.hadamard_product(A*witness, B*witness)
        c_temp = []

        for i in range(num_rows):
            loc = int(random() * (num_columns - 1))
            while witness[loc] == 0:
                loc = random() * num_columns

            witval = witness[loc]
            assert witval != 0

            this_row = sparse_row(num_columns)
            this_row[loc] = 0
            dot_value = vector(Fq, this_row).inner_product(witness)
            entry = (expected_product[i] - dot_value) / witval
            this_row[loc] = entry
            assert vector(this_row).inner_product(witness) == expected_product[i]
            c_temp.append(this_row)

        C = Matrix(Fq, c_temp)

        # print(f"A=\n{A}\n")
        # print(f"B=\n{B}\n")
        # print(f"C=\n{C}\n")
        # print(f"W={witness}\n")
        # print(f"Hadamard: {utils.hadamard_product(A*witness, B*witness)} == {C*witness}\n")

        return klass(A, B, C, witness)


    def __init__(self, A, B, C, w) -> None:
        self.A = A
        self.B = B
        self.C = C
        self.w = vector(w)
        assert utils.hadamard_product(self.A*self.w, self.B*self.w) == self.C*self.w

    def is_satisfiable(self) -> bool:
        return vector(utils.hadamard_product(self.A*self.w, self.B*self.w)) == self.C*self.w

    def __repr__(self) -> str:
        marker = lambda v, M : f"\n========= {v} =========\n\n{M}\n"
        return f"{marker("A", self.A)}{marker("B", self.B)}{marker("C", self.C)}{marker("W", self.w)}"



if __name__ == '__main__':

    def test_r1cs_instance_gen():
        Fq = GF(15*(2**27) + 1)
        num_rows = 7
        num_columns = 4
        instance = R1CS.random_element(Fq, num_rows, num_columns);
        assert instance.is_satisfiable()

        random_witness = random_vector(Fq, 9)
        instance = R1CS.random_element(Fq, num_rows, force_witness=random_witness, sparsity=0.2)
        print(f"{instance}")
        assert instance.is_satisfiable()


    def main():
        test_r1cs_instance_gen()

    main()
from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def swap_rows(self, row1, row2):
        temp = self[row1]
        self[row1] = self[row2]
        self[row2] = temp

    def multiply_coefficient_and_row(self, coefficient, row):
        self[row] = Plane(self[row].normal_vector.times_scalar(coefficient), self[row].constant_term * coefficient)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        added_to = self[row_to_be_added_to]
        add = self[row_to_add]
        self[row_to_be_added_to] = Plane(added_to.normal_vector.plus(add.normal_vector.times_scalar(coefficient)),
                                         added_to.constant_term + (add.constant_term * coefficient))

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)

        indices = [-1] * num_equations

        for i, p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def compute_triangular_form(self):
        system = deepcopy(self)

        current_index = 0
        dimensions = system[0].dimension
        for i in range(0, len(system), 1):
            coefficient = MyDecimal(system[i][current_index])
            if coefficient.is_near_zero():              # If this row is unsuitable for triangular form
                for j in range(i+1, len(system), 1):    # Search the rows below for a non-zero index in this position
                    if system[j][current_index] != 0:
                        system.swap_rows(i, j)
                        break
            system.clear_rows_below(current_index, i)
            current_index += 1
            if current_index == dimensions:
                break

        return system

    def clear_rows_below(self, index, row):
        for i in range(row + 1, len(self), 1):                              # For each row below the reference one
            coefficient = MyDecimal(self[i][index])
            if not coefficient.is_near_zero():                              # Has to be cleared
                reference_coefficient = self[row][index]
                coefficient_to_clear = self[i][index]
                factor = coefficient_to_clear / reference_coefficient * -1  # -1 to turn it into a substraction

                self.add_multiple_times_row_to_row(factor, row, i)

    def clear_rows_above(self, index, row):
        for i in range(row - 1, -1, -1):                              # For each row below the reference one
            coefficient = MyDecimal(self[i][index])
            if not coefficient.is_near_zero():                              # Has to be cleared
                reference_coefficient = self[row][index]
                coefficient_to_clear = self[i][index]
                factor = coefficient_to_clear / reference_coefficient * -1  # -1 to turn it into a substraction

                self.add_multiple_times_row_to_row(factor, row, i)

    def compute_rref(self):
        tf = self.compute_triangular_form()

        nonzeros = tf.indices_of_first_nonzero_terms_in_each_row()

        for i in range(len(tf) - 1, -1, -1):                    # Go from the last equation up
            if nonzeros[i] == -1:                               # Skip the equation if it doesn't have coeffiecients
                continue

            index_to_check = nonzeros[i]
            tf.multiply_coefficient_and_row(Decimal(1.0) / tf[i][index_to_check], i)
            tf.clear_rows_above(index_to_check, i)

        return tf

    def solve(self):
        rref = self.compute_rref()

        has_no_solution = False
        valid_equations = 0
        dimensions = rref[0].dimension
        solution = []
        nonzeros = rref.indices_of_first_nonzero_terms_in_each_row()

        for i in range(0, len(nonzeros), 1):
            if nonzeros[i] != -1:
                solution.insert(nonzeros[i], rref[i].constant_term)
                valid_equations += 1
            elif not MyDecimal(rref[i].constant_term).is_near_zero():
                has_no_solution = True
                break

        if has_no_solution:
            return None
        elif valid_equations < dimensions:
            return self.parametrized_solve()
        else:
            return Vector(solution)

    # Used only when system has infinite solutions, reproduces manual approach of solving the parametrization
    # This method should be modularized into smaller components, but it's late and I want to sleep :D
    def parametrized_solve(self):
        rref = self.compute_rref()
        total_dimensions = rref[0].dimension
        intermediate_vectors = []
        nonzeros = rref.indices_of_first_nonzero_terms_in_each_row()

        # First I create a set of equations equal to the original ones but "solving" for the pivot variable
        #   note: I only reproduce the right side of the equation in a Vector
        #   e.g.
        #   x + y +  z = 1   =>   x = 1 -  y -  z   => [1, -1, -1]
        #       y - 2z = 2   =>   y = 2 + 0y + 2z   => [2, 0, 2]
        for i in range(0, len(rref), 1):                        # For every equation
            nonzero_for_eq = nonzeros[i]
            if nonzero_for_eq != -1:                            # If it is a valid equation
                int_vector = [rref[i].constant_term]            # Append the constant term at the start of our vector
                for j in range(1, total_dimensions, 1):         # Skip x as x will never be a free variable
                    if j != nonzero_for_eq:                     # Subtract the corresponding coefficient unless it is
                                                                # the coefficient for the pivot variable being evaluated
                        int_vector.append(rref[i][j] * -1)
                    else:
                        int_vector.append(0)                    # If it is the pivot variable then append a 0
                intermediate_vectors.append(int_vector)

        # Now I add equations for the free variables and "missing equations" (i.e. eqs in 0 = 0 form)
        last_nonzero_index = -1
        current_equation = 0
        true_results = []
        for i in range(0, len(nonzeros), 1):
            if i < total_dimensions:
                while nonzeros[i] > last_nonzero_index + 1:                 # There might be several missing pivots
                    free_variable_vector = [0] * total_dimensions
                    free_variable_vector[i] = 1
                    true_results.append(free_variable_vector)
                    last_nonzero_index += 1

                if nonzeros[i] == -1:
                    free_variable_vector = [0] * total_dimensions
                    free_variable_vector[i] = 1
                    true_results.append(free_variable_vector)
                else:
                    true_results.append(intermediate_vectors[current_equation])

                current_equation += 1
                last_nonzero_index = nonzeros[i]

        # Create the "verticals" from the solution (grab all constants and
        # coefficients for each variable in the same vector)
        verticals = []
        for i in range(0, total_dimensions, 1):
            vertical = []
            for j in range(0, len(true_results), 1):
                value = 0
                if not MyDecimal(true_results[j][i]).is_near_zero():
                    value = true_results[j][i]
                vertical.append(Decimal(value))
            verticals.append(Vector(vertical))

        # Construct final Parametrization object with the vectors we produced
        return Parametrization(verticals[0], verticals[1:])

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1, p) for i, p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


class Parametrization(object):

    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM = (
        'The basepoint and direction vectors should all live in the same '
        'dimension')

    def __init__(self, basepoint, direction_vectors):

        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.basepoint.dimension

        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension

        except AssertionError:
            raise Exception(self.BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM)

    def __str__(self):

        output = ''
        for coord in range(self.dimension):
            output += 'x_{} = {} '.format(coord + 1,
                                          round(self.basepoint[coord], 3))
            for free_var, vector in enumerate(self.direction_vectors):
                output += '+ {} t_{}'.format(round(vector[coord], 3),
                                             free_var + 1)
            output += '\n'
        return output


p1 = Plane(normal_vector=Vector([0.786, 0.786, 0.588]), constant_term=-0.714)
p2 = Plane(normal_vector=Vector([-0.131, -0.131, 0.244]), constant_term=0.319)
s = LinearSystem([p1, p2])
r = s.parametrized_solve()
print r

p1 = Plane(normal_vector=Vector([8.631, 5.112, -1.816]), constant_term=-5.113)
p2 = Plane(normal_vector=Vector([4.315, 11.132, -5.27]), constant_term=-6.775)
p3 = Plane(normal_vector=Vector([-2.158, 3.01, -1.727]), constant_term=-0.831)
s = LinearSystem([p1, p2, p3])
r = s.parametrized_solve()
print r

p1 = Plane(normal_vector=Vector([0.935, 1.76, -9.365]), constant_term=-9.955)
p2 = Plane(normal_vector=Vector([0.187, 0.352, -1.873]), constant_term=-1.991)
p3 = Plane(normal_vector=Vector([0.374, 0.704, -3.746]), constant_term=-3.982)
p4 = Plane(normal_vector=Vector([-0.561, -1.056, 5.619]), constant_term=5.973)
s = LinearSystem([p1, p2, p3, p4])
r = s.parametrized_solve()
print r

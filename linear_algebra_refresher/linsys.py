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


p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='-1') and
        r[1] == p2):
    print 'test case 1 failed'
else:
    print 'test case 1 succeeded'

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == p1 and
        r[1] == Plane(constant_term='1')):
    print 'test case 2 failed'
else:
    print 'test case 2 succeeded'

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
s = LinearSystem([p1,p2,p3,p4])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='0') and
        r[1] == p2 and
        r[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
        r[3] == Plane()):
    print 'test case 3 failed'
else:
    print 'test case 3 succeeded'

p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
s = LinearSystem([p1,p2,p3])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term=Decimal('23')/Decimal('9')) and
        r[1] == Plane(normal_vector=Vector(['0','1','0']), constant_term=Decimal('7')/Decimal('9')) and
        r[2] == Plane(normal_vector=Vector(['0','0','1']), constant_term=Decimal('2')/Decimal('9'))):
    print 'test case 4 failed'
else:
    print 'test case 4 succeeded'

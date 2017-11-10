from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class Plane(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.basepoint = None
        self.set_basepoint()

    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Plane.first_nonzero_index(n.coordinates)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            inner_output = ''

            if coefficient < 0:
                inner_output += '-'
            if coefficient > 0 and not is_initial_term:
                inner_output += '+'

            if not is_initial_term:
                inner_output += ' '

            if abs(coefficient) != 1:
                inner_output += '{}'.format(abs(coefficient))

            return inner_output

        n = self.normal_vector

        try:
            initial_index = Plane.first_nonzero_index(n.coordinates)
            terms = [write_coefficient(n[i], is_initial_term=(i == initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output

    def __eq__(self, other):
        if self.normal_vector.is_zero():
            if not other.normal_vector.is_zero():
                return False
            else:
                diff = self.constant_term - other.constant_term
                return MyDecimal(diff).is_near_zero()
        elif other.normal_vector.is_zero():
            return False

        if not self.is_parallel_to(other):
            return False

        x = self.basepoint
        y = other.basepoint
        connection = y.minus(x)

        return connection.is_orthogonal_to(self.normal_vector)

    def __getitem__(self, index):
        return self.normal_vector[index]

    def is_parallel_to(self, other):
        n1 = self.normal_vector.normalized()
        n2 = other.normal_vector.normalized()

        return n1.is_parallel_to(n2)

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


def main():
    first = Vector([-0.412, 3.806, 0.728])
    k1 = -3.46
    second = Vector([1.03, -9.515, -1.82])
    k2 = 8.65

    a = Plane(first, k1)
    b = Plane(second, k2)
    if a == b:
        print '1st pair is equal'
    elif a.is_parallel_to(b):
        print '1st pair is parallel but not equal'
    else:
        print '1st pair is not parallel'

    first = Vector([2.611, 5.528, 0.283])
    k1 = 4.6
    second = Vector([7.715, 8.306, 5.342])
    k2 = 3.76

    a = Plane(first, k1)
    b = Plane(second, k2)
    if a == b:
        print '2nd pair is equal'
    elif a.is_parallel_to(b):
        print '2nd pair is parallel but not equal'
    else:
        print '2nd pair is not parallel'

    first = Vector([-7.926, 8.625, -7.212])
    k1 = -7.952
    second = Vector([-2.642, 2.875, -2.404])
    k2 = -2.443

    a = Plane(first, k1)
    b = Plane(second, k2)
    if a == b:
        print '3rd pair is equal'
    elif a.is_parallel_to(b):
        print '3rd pair is parallel but not equal'
    else:
        print '3rd pair is not parallel'


if __name__ == "__main__":
    main()
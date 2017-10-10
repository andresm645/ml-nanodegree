from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class Line(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2

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

            initial_index = Line.first_nonzero_index(n.coordinates)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
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
            initial_index = Line.first_nonzero_index(n.coordinates)
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

    def is_parallel_to(self, other):
        n1 = self.normal_vector.normalized()
        n2 = other.normal_vector.normalized()

        return n1.is_parallel_to(n2)

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

    def intersection_with(self, other):
        if self.is_parallel_to(other):
            if self == other:
                return self
            else:
                return None
        else:

            a = self.normal_vector[0]
            b = self.normal_vector[1]
            c = other.normal_vector[0]
            d = other.normal_vector[1]
            k1 = self.constant_term
            k2 = other.constant_term

            x = ((d * k1) - (b * k2)) / ((a * d) - (b * c))
            y = ((-1 * c * k1) + (a * k2)) / ((a * d) - (b * c))

            return Vector([x, y])

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


def main():
    first = Vector([4.046, 2.836])
    k1 = 1.21
    second = Vector([10.115, 7.09])
    k2 = 3.025

    a = Line(first, k1)
    b = Line(second, k2)
    c = a.intersection_with(b)
    print 'Intersection 1: ' + c.__str__()

    first = Vector([7.204, 3.182])
    k1 = 8.68
    second = Vector([8.172, 4.114])
    k2 = 9.883

    a = Line(first, k1)
    b = Line(second, k2)
    c = a.intersection_with(b)
    print 'Intersection 2: ' + c.__str__()

    first = Vector([1.182, 5.562])
    k1 = 6.744
    second = Vector([1.773, 8.343])
    k2 = 9.525

    a = Line(first, k1)
    b = Line(second, k2)
    c = a.intersection_with(b)
    print 'Intersection 3: ' + c.__str__()


if __name__ == "__main__":
    main()

from math import sqrt, acos, pi
from decimal import Decimal, getcontext
import pdb

getcontext().prec = 30

class Vector(object):
    CANNOT_NORMALIZE_ZERO_VECTOR = 'Cannot normalize the zero vector'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple(Decimal(x) for x in coordinates)
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def plus(self, v):
        if self.dimension != other.dimension:
            raise ValueError('Vectors must have the same number of dimensions to be added')

        new_coordinates = [x+y for x,y in zip(self.coordinates, other.coordinates)]
        return Vector(new_coordinates)

    def minus(self, other):
        if self.dimension != other.dimension:
            raise ValueError('Vectors must have the same number of dimensions to be substracted')

        new_coordinates = [x-y for x,y in zip(self.coordinates, other.coordinates)]
        return Vector(new_coordinates)

    def times_scalar(self, c):
        new_coordinates = [Decimal(c) * x for x in self.coordinates]
        return Vector(new_coordinates)

    def magnitude(self):
        coordinates_squared = [x**2 for x in self.coordinates]
        return sqrt(sum(coordinates_squared))

    def normalized(self):
        try:
            magnitude = self.magnitude()
            return self.times_scalar(Decimal(1)/Decimal(magnitude))

        except ZeroDivisionError:
            raise ValueError(CANNOT_NORMALIZE_ZERO_VECTOR)

    def dot(self, other):
        zipped_products = [x * y for (x, y) in zip(self.coordinates, other.coordinates)]
        return sum(zipped_products)

    def angle_with(self, other, in_degrees = False):
        #try:
        u1 = self.normalized()
        u2 = other.normalized()
        dot_product = min(max(u1.dot(u2), -1), 1)

        #pdb.set_trace()
        angle_in_radians = acos(dot_product)

        if in_degrees:
            degrees_per_radian = 180. / pi
            return angle_in_radians * degrees_per_radian

        else:
            return angle_in_radians

        #except Exception as e:
        #    if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR:
        #        raise ValueError('Cannot compute an angle with the zero vector')
        #    else:
        #        raise e

    def is_zero(self, tolerance=1e-10):
        return self.magnitude() < tolerance

    def is_parallel_to(self, other):
        return (self.is_zero() or 
                other.is_zero() or
                self.angle_with(other) == 0 or
                Decimal(self.angle_with(other)) == Decimal(pi))

    def is_orthogonal_to(self, other, tolerance=1e-10):
        return abs(self.dot(other)) < tolerance

    def parallel_component_to(self, other):
        normalized_base = other.normalized()
        scalar = self.dot(normalized_base)
        return normalized_base.times_scalar(scalar)

    def orthogonal_component_to(self, other):
        parallel_component = self.parallel_component_to(other)
        return self.minus(parallel_component)




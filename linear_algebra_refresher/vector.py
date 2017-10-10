from math import sqrt, acos, pi
from decimal import Decimal, getcontext

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

    def __getitem__(self, index):
        return self.coordinates[index]

    def plus(self, v):
        if self.dimension != v.dimension:
            raise ValueError('Vectors must have the same number of dimensions to be added')

        new_coordinates = [x + y for x, y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def minus(self, other):
        if self.dimension != other.dimension:
            raise ValueError('Vectors must have the same number of dimensions to be substracted')

        new_coordinates = [x - y for x, y in zip(self.coordinates, other.coordinates)]
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

    def angle_with(self, other, in_degrees=False):
        try:
            u1 = self.normalized()
            u2 = other.normalized()
            dot_product = min(max(u1.dot(u2), -1), 1)

            angle_in_radians = acos(dot_product)

            if in_degrees:
                degrees_per_radian = 180. / pi
                return angle_in_radians * degrees_per_radian

            else:
                return angle_in_radians

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR:
                raise ValueError('Cannot compute an angle with the zero vector')
            else:
                raise e

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

    def cross(self, other):
        if self.dimension != 3 or other.dimension != 3:
            raise ValueError('Vectors must have 3 dimensions to calculate cross product')

        coor_1 = self.coordinates
        coor_2 = other.coordinates
        first = coor_1[1] * coor_2[2] - coor_2[1] * coor_1[2]
        second = -1 * (coor_1[0] * coor_2[2] - coor_2[0] * coor_1[2])
        third = coor_1[0] * coor_2[1] - coor_2[0] * coor_1[1]
        return Vector([first, second, third])

    def area_of_parallelogram_with(self, other):
        return self.cross(other).magnitude()

    def area_of_triangle_with(self, other):
        return self.area_of_parallelogram_with(other) / 2.0
        return self.area_of_parallelogram_with(other) / 2.0


def main():
    a = Vector([8.462, 7.893, -8.187])
    b = Vector([6.984, -5.975, 4.778])
    c = a.cross(b)
    print 'Cross: ' + c.__str__()

    a = Vector([-8.987, -9.838, 5.031])
    b = Vector([-4.268, -1.861, -8.866])
    c = a.area_of_parallelogram_with(b)
    print 'Parallelogram: ' + c.__str__()

    a = Vector([1.5, 9.547, 3.691])
    b = Vector([-6.007, 0.124, 5.772])
    c = a.area_of_triangle_with(b)
    print 'Triangle: ' + c.__str__()


if __name__ == "__main__":
    main()

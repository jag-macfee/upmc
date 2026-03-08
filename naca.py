import numpy as np

# Resolver to return functions defining NACA airfoil mean camber lines
# Required to generate relevant points to define panel approximations along
# a non-flat plate shaped airfoil


class Airfoil:
    """Base class for an airfoil."""

    def __init__(self, c):
        self.c = c

    def camber(self, x):
        """Makes airfoil objects callable, returning the y-coordinate of the
        mean camber line at a given x-station."""
        raise NotImplementedError("Subclasses should implement this method.")


class FlatPlate(Airfoil):
    """
    Represents a flat plate airfoil.
    Fallback when no digits or unsupported digit structure is supplied to naca() function
    """

    def camber(self, x):
        return 0


class Naca4Digit(Airfoil):
    """Represents a 4-digit NACA airfoil."""

    def __init__(self, digits: str, c):
        super().__init__(c)
        self.m = c * float(digits[0]) / 100
        self.p = c * float(digits[1]) / 10

    def camber(self, x):
        # Check bounds for undefined behaviour
        if x < 0 or x > self.c:
            raise ValueError(
                f"Position {x} is undefined for airfoil with chord length {self.c}"
            )

        # Case: x is before max camber
        if x >= 0 and x <= self.p:
            return self.m / (self.p**2) * (2 * self.p * x - x**2)

        # Case: x is after max camber
        return self.m / ((1 - self.p) ** 2) * ((1 - 2 * self.p) + 2 * self.p * x - x**2)


# Given a string code (eg. "2415"), return a function
# defining the mean camber line
# Note: This function will take into account absolute (non-normalised) values
# which will be relevant when calculating time stepped solutions in unsteady panel method
#
# Reference: https://web.stanford.edu/~cantwell/AA200_Course_Material/The%20NACA%20airfoil%20series.pdf
def naca(digits: str, c) -> Airfoil:
    # 4-digit
    if len(digits) == 4:
        return Naca4Digit(digits, c)

    # Default fallback
    return FlatPlate(c)

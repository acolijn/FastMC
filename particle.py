import numpy as np
import pandas as pd

class particle:
    def __init__(self, **kwargs):
        """
        Initialize a particle
        :param kwargs:
            type = particle specimen (value = gamma, <others defined later>)
            energy = energy of the particle in keV
        """
        print("particle::initialize")

        self.type = kwargs.pop('type', None)
        self.energy = kwargs.pop('energy', 0.0)

        self.x0 = np.zeros(3)
        self.direction = np.zeros(3)

    def generate(self, cylinder):
        """
        Generate the starting point of the particle and its direction
        :param cylinder:
        :return:
        """
        # generate x0 and direction of the particle
        self.x0 = cylinder.generate_point()['x']

        #
        cost = np.random.uniform(-1, 1)
        sint = np.sqrt(1 - cost ** 2)
        phi = 2 * np.pi * np.random.uniform(0, 1)

        tx = np.cos(phi) * sint
        ty = np.sin(phi) * sint
        tz = cost
        self.direction = np.array([tx, ty, tz])

        return

    def intersect(self, cylinder):
        """
        Intersect particle track with cylinder

        :param cylinder: instant of teh cylinder class, containing the geo definition
        :return number of intersections with the cylinder:
        0=particle does not hit cylinder
        1=particle probably originates somewhere on the surface, but moving away from the cylinder
        2=particle hits the cylinder at two spots

        """
        self.intersections = pd.DataFrame()

        # routines below can only be called if self.intersections is defined above
        self.intersect_with_plane(cylinder, type='top')
        self.intersect_with_plane(cylinder, type='bot')
        self.intersect_with_side(cylinder)

        self.intersections = self.intersections.sort_values('s')
        print(self.intersections)

        return len(self.intersections)

    def intersect_with_side(self, cylinder):
        """
        Intersect the particle with the cylinder shell
        :param cylinder:
        :return:
        """

        A = self.direction[0] ** 2 + self.direction[1] ** 2
        B = 2 * (self.x0[0] * self.direction[0] + self.x0[1] * self.direction[1])
        C = self.x0[0] ** 2 + self.x0[1] ** 2 - cylinder.radius ** 2

        discriminant = B ** 2 - 4 * A * C

        if discriminant >= 0:
            for sign in (-1,1):
                s = (-B + sign*np.sqrt(discriminant)) / (2 * A)
                xint = self.x0 + s * self.direction
                # is it hitting the cylinder? or is it outside?
                if (np.abs(xint[2]) < cylinder.height/2) & (s>-1e-5): # only tracks with positive pathlength, but allow for numerical -0.0000
                    # good intersection..... add to the list
                    self.intersections = self.intersections.append({'x': xint, 's': s}, ignore_index=True)

        return

    def intersect_with_plane(self, cylinder, **kwargs):
        """
        Intersect the particle track with the top/bottom plane
        :param cylinder:
        :param kwargs:

        :return:
        """

        type = kwargs.pop('type', None)

        zint = 0
        if type == "top":  # top plane
            zint = cylinder.height / 2
        else:  # bottom plane
            zint = -cylinder.height / 2

        tz = self.direction[2]

        # calculate the path length to teh intersection point
        if tz != 0.000:
            s = (zint - self.x0[2]) / tz
        else:
            print("particle::intersect_with_plane() WARNING: tz = ", tz)
            s = 0

        # calculate the position of the intersection point from the track equation
        x = self.x0 + s * self.direction

        # calculate the radius  of teh intersect and check whether it is in range
        r = np.sqrt(x[0] ** 2 + x[1] ** 2)
        if (r < cylinder.radius) & (s>-1e-5): # only tracks with positive pathlength, but allow for numerical -0.0000
            # good intersection.... add to the list
            self.intersections = self.intersections.append({'x': x, 's': s}, ignore_index=True)

        return

    def propagate(self, phys):
        """
        Propagate the particle to the next interaction location
        :param phys: physics class object
        :return:
        """
        # cdf for path length
        mu = phys.get_att(E=self.energy)
        r = np.random.uniform(0,1)
        L =  - np.log(1-r) * mu

        return L

    def scatter(self):
        """

        :return:
        """

        return

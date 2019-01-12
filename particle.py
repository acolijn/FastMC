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
        self.phys = kwargs.pop('physics',None)
        self.cryostat = kwargs.pop('geometry',None)
        self.fiducial = kwargs.pop('ficucial',None)

        self.x0 = np.zeros(3)
        self.direction = np.zeros(3)


    def generate(self):
        """
        Generate the starting point of the particle and its direction
        :param:
        :return:
        """
        # generate x0 and direction of the particle
        self.x0 = self.cryostat.generate_point()['x']

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
        intersections = pd.DataFrame()

        # routines below can only be called if self.intersections is defined above
        intersections = intersections.append(self.intersect_with_plane(cylinder, type='top'))
        intersections = intersections.append(self.intersect_with_plane(cylinder, type='bot'))
        intersections = intersections.append(self.intersect_with_side(cylinder))

        intersections = intersections.sort_values('s')
        intersections = intersections.reset_index(drop=True)

        # print(self.intersections)

        return intersections

    def intersect_with_side(self, cylinder):
        """
        Intersect the particle with the cylinder shell
        :param cylinder:
        :return: intersections
        """

        intersections = pd.DataFrame()

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
                    intersections = intersections.append({'x': xint, 's': s}, ignore_index=True)

        return intersections

    def intersect_with_plane(self, cylinder, **kwargs):
        """
        Intersect the particle track with the top/bottom plane
        :param cylinder:
        :param kwargs:

        :return: intersections
        """
        intersections = pd.DataFrame()

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
            intersections = intersections.append({'x': x, 's': s}, ignore_index=True)

        return intersections

    def propagate(self):
        """
        Propagate the particle
        :return:
        """

        # 1. check if the particle moves through the xenon
        intersections = self.intersect(self.cryostat)
        if len(intersections) == 2: # we have tracks through the LXe
            # maximum path length
            s_max = intersections['s'][1]
            s_gen = self.generate_interaction_point()

            print("s_gen = ",s_gen,"s_max =",s_max)
            if(s_gen < s_max): # we have a hit inside the xenon
                self.scatter()

    def generate_interaction_point(self):
        """
        Propagate the particle to the next interaction location
        :param phys: physics class object
        :return:
        """
        # cdf for path length
        mu = self.phys.get_att(energy=self.energy)
        r = np.random.uniform(0,1)
        L =  - np.log(1-r) * mu

        return L

    def scatter(self):
        """

        :return:
        """

        # select the scatter process
        process = self.select_scatter_process()

        if process == 'inc':
            print('Compton scatter - continue')
            # select the scatter angle

            # recalculate the new particle direction

            # calculate teh energy deposit in the xenon
        else:
            # particle is fully absorbed
            print('PE absoprtion - end particle')

        return

    def select_scatter_process(self):
        """

        :return:
        """

        # get the gamma energy
        E = self.energy
        # if we deal with Compton we will scatter, otherwise full absorption (not completely true, but I dont want to deal with pair creation)
        sigma_total = self.phys.get_sigma(process='att',energy=E)
        sigma_inc = self.phys.get_sigma(process='inc',energy=E)
        frac = sigma_inc / sigma_total

        r = np.random.uniform(0,1)
        print('stot = ',sigma_total,' sinc = ',sigma_inc,' frac = ',frac,' r= ',r)
        if r < frac:
            process = 'inc'
        else:
            process = 'pho'

        return process


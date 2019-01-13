import numpy as np
import pandas as pd

from numba import jit

class particle:
    #@jit()

    def __init__(self, **kwargs):
        """
        Initialize a particle
        :param kwargs:
            type = particle specimen (value = gamma, <others defined later>)
            energy = energy of the particle in keV
        """
        # # # print("particle::initialize")

        self.type = kwargs.pop('type', None)
        self.energy = kwargs.pop('energy', 0.0)
        self.phys = kwargs.pop('physics',None)
        self.cryostat = kwargs.pop('geometry',None)
        self.fiducial = kwargs.pop('ficucial',None)

        self.x0 = np.zeros(3)
        self.direction = np.zeros(3)

        # generate the x0 and direction of the particle
        self.generate()

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
        stop = self.intersect_with_plane(cylinder, 'top')
        sbot = self.intersect_with_plane(cylinder, 'bot')
        ssid = self.intersect_with_side(cylinder)
        # sort the list by ascending order....
        intersections = sorted([stop,sbot,ssid[0],ssid[1]],reverse=False)
        # remove s=0 entries
        for i in range(4):
            try:
                intersections.remove(0.)
            except:
                pass

        return intersections

    def intersect_with_side(self, cylinder):
        """
        Intersect the particle with the cylinder shell
        :param R,h:
        :return: intersections
        """

        intersections = [0,0]
        n = 0

        A = self.direction[0] ** 2 + self.direction[1] ** 2
        B = 2 * (self.x0[0] * self.direction[0] + self.x0[1] * self.direction[1])
        C = self.x0[0] ** 2 + self.x0[1] ** 2 - cylinder.radius ** 2
        #C = self.x0[0] ** 2 + self.x0[1] ** 2 - R ** 2

        discriminant = B ** 2 - 4 * A * C

        if discriminant >= 0:
            for sign in (-1,1):
                s = (-B + sign*np.sqrt(discriminant)) / (2 * A)
                xint = self.x0 + s * self.direction
                # is it hitting the cylinder? or is it outside?
                if (np.abs(xint[2]) < cylinder.height/2) & (s>1e-5): # only tracks with positive pathlength, remove intersect with zero pathlength
                #if (np.abs(xint[2]) < h / 2) & (s > 1e-5):  # only tracks with positive pathlength, remove intersect with zero pathlength
                    # good intersection..... add to the list
                    #intersections = intersections.append({'x': xint, 's': s}, ignore_index=True)
                    intersections[n] = s
                    n = n+1
        return intersections

    def intersect_with_plane(self, cylinder, type):
        """
        Intersect the particle track with the top/bottom plane
        :param R,h,type:

        :return: s
        """
        sint = 0

        zint = 0
        if type == "top":  # top plane
            zint = cylinder.height / 2
            #zint = h / 2
        else:  # bottom plane
            zint = -cylinder.height / 2
            #zint = -h / 2


        tz = self.direction[2]

        # calculate the path length to the intersection point
        if tz != 0.000:
            s = (zint - self.x0[2]) / tz
        else:
            print("particle::intersect_with_plane() WARNING: tz = ", tz)
            s = 0

        # calculate the position of the intersection point from the track equation
        x = self.x0 + s * self.direction

        # calculate the radius  of teh intersect and check whether it is in range
        r = np.sqrt(x[0] ** 2 + x[1] ** 2)
        if (r < cylinder.radius) & (s>1e-5): #  only tracks with positive pathlength, remove intersect with zero pathlength
        #if (r < R) & (s > 1e-5):  # only tracks with positive pathlength, remove intersect with zero pathlength
            # good intersection.... add to the list
            sint = s

        return sint

    def propagate(self):
        """
        Propagate the particle
        :return:
        """
        cost = -1.5

        # 1. check if the particle moves through the xenon
        intersections = self.intersect(self.cryostat)

        if len(intersections) >= 1: # we have a track through the LXe cylinder.
            # 1. if one intersection this gives the exit point.
            # 2. if two intersection this gives the entry point to a new volume (needed for fiducial)
            #    volume variance reduction
            s_max = intersections[0]
            s_gen = self.generate_interaction_point()
            # print("s_gen = ",s_gen,"s_max =",s_max)
            if(s_gen < s_max): # we have a hit inside the xenon
               cost = self.scatter()
            # cost = 1
        return cost

    @jit()
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

    @jit()
    def scatter(self):
        """

        :return:
        """

        # select the scatter process
        process = self.select_scatter_process()

        cost = -1.5
        if process == 'inc':
            # # #print('Compton scatter - continue')
            # select the scatter angle

            # make the cdf from the differential cross section
            cost_range = np.arange(-1.0,+1.0,0.01)
            dsigma = self.phys.KleinNishina(self.energy,cost_range)
            cdf = dsigma.cumsum()/dsigma.sum()
            # draw a random number from the dsigma distribution
            cost = np.interp(np.random.uniform(0.,1.), cdf, cost_range)

            # calculate the rotation matrix

            # recalculate the new particle direction


            # calculate teh energy deposit in the xenon
        else:
            # particle is fully absorbed
            # # #print('PE absoprtion - end particle')
            i=0

        return cost

    @jit()
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
        # # # print('stot = ',sigma_total,' sinc = ',sigma_inc,' frac = ',frac,' r= ',r)
        if r < frac:
            process = 'inc'
        else:
            process = 'pho'

        return process


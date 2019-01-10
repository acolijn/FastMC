import numpy as np

class em_physics:
    """
    The electromagnetic physics driver for the simulation
    (routines routinely stole from https://github.com/ErikHogebirk/DMPlots)
    """

    def __init__(self, **kwargs):

        d = self.extract_nist('data/xenon_cross_section/gamma_sigma.txt')
        self.rho = 3.0 # gram/cm3
        self.e = d[:, 0] # energy range
        self.sigma_coh = d[:, 1]  # Coherent
        self.sigma_inc = d[:, 2]  # Incoherent. This is Compton scattering
        self.sigma_pho = d[:, 3]  # photoelectric absorption
        self.sigma_pp_nuc = d[:, 4]  # pair production nuclear field
        self.sigma_pp_el = d[:, 5]  # pair production electron field (order of magnitude below)
        self.sigma_pp = self.sigma_pp_nuc + self.sigma_pp_el
        self.sigma_att = d[:, 7]  # Total attenuation without coherent scattering
        self.att = 1 / (self.rho * self.sigma_att)  # Attenuation length

        return

    def extract_nist(self,fn):
        """
        Extract the cross section and attenuation from the data file (stolen from Erik Hogenbirk)

        :param fn: file containing the cross section information for xenon from NIST
        :return: cross section data
        """
        all_params = []
        with open(fn, 'r') as f:
            for i, line in enumerate(f):
                # Skip header
                if i < 6:
                    continue
                # Skip blank lines
                if len(line) < 5:
                    continue
                # Cut off the shell indication and newline character
                line = line[7:-1]
                params = [np.float(el) for el in line.split(sep=' ')]

                all_params.append(params)

        return np.array(all_params)

    def get_sigma(self, **kwargs):
        """

        :param kwargs:
        process=["att" = total without coherent,
                 "pho" = PE absorption,
                 "inc" = incoherent scatter - Compton,
                 "pp"  = pair creation
        E=energy of gamma in keV,

        :return: cross section in cm2/g
        """
        process = kwargs.pop('process','att')
        energy = kwargs.pop('E',-1.0)

        if process == "att":  # total cross section
            mu = np.interp(energy / 1e3, self.e, self.sigma_att)
        elif process == "pho": # PE absorption
            mu = np.interp(energy / 1e3,self.e,self.sigma_pho)
        elif process == "inc": # Compton
            mu = np.interp(energy / 1e3,self.e,self.sigma_inc)
        elif process == "pp": # pair creation
            mu = np.interp(energy / 1e3,self.e,self.sigma_pp)
        else:
            print('em_physics::get_att ERROR wrong process selected')

        return mu

    def get_att(self,**kwargs):
        """
        Calculate the attenuation length

        :param kwargs:
            E=energy of gamma in keV
        :return: attenuation length in cm
        """
        energy = kwargs.pop('E',-1.0)

        mu = np.interp(energy * 1e3, self.e, self.sigma_att)
        return 1/ ( mu * self.rho)


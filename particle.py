# import cylinder as cyl

class particle:
    def __init__(self, **kwargs):
        print("particle::initialize")

        self.type = kwargs.pop('type',None)
        self.energy = kwargs.pop('energy',0.0)

        self.x0 = [0,0,0]
        self.direction = [0,0,0]

    def generate(self, cylinder):
        # generate x0 and direction of the particle
        self.x0 = cylinder.generate_point()['x']

        return

    def intersect(self):
        # intersect with cylinder
        return

    def scatter(self):
        # scatter the particle and calculate the energy deposit
        return

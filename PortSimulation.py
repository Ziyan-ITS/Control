import random
import math


class Vessel(object):
    arrival_rate = 2

    def __init__(self, nums):
        self.V = nums
        self.vessels = {'vessel_{}'.format(i): {} for i in range(nums)}
        self.generate_arrival_time(self.arrival_rate)
        self.generate_length([10, 20, 30])
        self.generate_container_volume([100, 200, 300])

    def generate_arrival_time(self, arrival_rate):
        """
        :param arrival_rate: lam:average vessels arrival rate (nums/h)
        :return:
        """
        # automated generate the arrival time by bosong distribution
        arrival_distribution = [1 - math.exp(-arrival_rate * (time / 10)) for time in range(0, 60)]
        arrival_distribution[-1] = 1
        arrival_time = 0

        def find_index(pro, array):
            for index, value in enumerate(array):
                if value > pro:
                    return index

        for i in range(self.V):
            probability = random.uniform(0, 1)
            ind = find_index(probability, arrival_distribution)
            if ind < 0.5:
                ind = 0.5
            arrival_time = arrival_time + ind
            name = 'vessel_{}'.format(i)
            self.vessels[name]['arrival_time'] = arrival_time

    def generate_length(self, length_type):
        # automated generate the length of arrival vessel, length here is measured in number of berth section
        rand_index = len(length_type)
        for key in self.vessels:
            self.vessels[key]['vessel_length'] = length_type[random.randint(0, rand_index)]

    def generate_container_volume(self, container_volume):
        # automated generate the container volume of arrival vessel
        rand_index = len(container_volume)
        for key in self.vessels:
            self.vessels[key]['container_volume'] = container_volume[random.randint(0, rand_index)]


class QuayCrane(object):
    def __init__(self, nums, Quay_length, *args):
        super(QuayCrane, self).__init__(*args)
        self.G = nums
        self.qc = {'qc_{}'.format(g): {} for g in range(nums)}
        self.B = Quay_length
        self.generate_action_range()
        self.generate_process_rate([10, 20, 30])

    def generate_process_rate(self, process_rate_distribution):
        rand = len(process_rate_distribution)
        for key in self.qc:
            self.qc[key]['processing_rate'] = process_rate_distribution[random.randint(0, rand)]

    def generate_action_range(self):
        i = self.B / (2 * self.G)
        for key in self.qc:
            self.qc[key]['action_range'] = (i - self.B / (2 * self.G), i + self.B / self.G)
            i = i + self.B / self.G

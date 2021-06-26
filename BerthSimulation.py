import random
import docplex.mp.model as cpx


class Vessel():
    def __init__(self, nums):
        self.nums = nums
        self.vessels = {'vessel_{}'.format(i): {} for i in range(nums)}

    def generate_arrival_time(self, arrival_rate):
        # automated generate the arrival time by bosong distribution
        arrival_time = []
        for i in range(self.nums):
            name = 'vessel_{}'.format(i)
            self.vessels[name]['arrival_time'] = arrival_time[i]
        pass

    def generate_length(self, length_type):
        # automated generate the length of arrival vessel, length here is measured in number of berth section
        rand_index = len(length_type)
        for key in self.vessels:
            self.vessels[key]['vessel_length'] = length_type[random.randint(0, rand_index)]
        pass

    def generate_container_volume(self, container_volume):
        pass


class QuayCrane():
    def __init__(self, nums, Quay_length):
        self.num = nums
        self.qc = {'qc_{}'.format(g): {} for g in range(nums)}
        self.Quayl = Quay_length

    def generate_process_rate(self, process_rate_distribution):
        rand = len(process_rate_distribution)
        for key in self.qc:
            self.qc[key]['processing_rate'] = process_rate_distribution[random.randint(0, rand)]

    def generate_action_range(self):
        i = 0
        for key in self.qc:
            self.qc[key]['action_range'] = (i, i + self.Quayl / self.num)
            i = i + self.Quayl / self.num


class Berth():
    def __init__(self, berth_lenrth, ship_nums, time_horizon, safety_time, quay_length, crane_nums):
        """
        :param berth_lenrth: The length of berth: 1unit/setion = 10 meter
        """
        self.opt_model = cpx.Model(name='BerthAllocation')
        self.B = berth_lenrth
        self.V = ship_nums
        self.T = time_horizon
        self.F = safety_time
        self.J = quay_length
        self.G = crane_nums

    def generate_variables(self):

        x_kl = {(k, l): self.opt_model.binary_var('x_{}{}'.format(k, l)) for k in range(self.V) for l in range(self.V)}
        y_kl = {(k, l): self.opt_model.binary_var('y_{}{}'.format(k, l)) for k in range(self.V) for l in range(self.V)}
        bk = {(k,): self.opt_model.integer_var(name='b_{}'.format(k)) for k in range(self.V)}
        tk = {(k,): self.opt_model.integer_var(name='t_{}'.format(k)) for k in range(self.V)}
        ck = {(k,): self.opt_model.integer_var(name='c_{}'.format(k)) for k in range(self.V)}
        z_gk_j = {(g, k, j): self.opt_model.binary_var('z_{}{}{}'.format(g, k, j)) for g in range(self.G) for k in
                  range(self.V) for j in range(self.T)}
        return x_kl, y_kl, bk, tk, ck

    def generate_discrete_pos(self):
        pi_kn = {'pi_{}{}'.format(k, n): self.opt_model.binary_var(name='pi_{}{}'.format(k, n)) for k in range(self.V)
                 for n in range(self.B)}
        theta_kn = {'theta_{}{}'.format(k, n): self.opt_model.binary_var(name='theta_{}{}'.format(k, n)) for k in
                    range(self.V) for n in range(self.B)}
        return pi_kn, theta_kn

    def generate_discrete_time(self):
        alpah_kj = {'alpha_{}{}'.format(k, j): self.opt_model.binary_var(name='alpha_{}{}.format(k, j)') for k in
                    range(self.V) for j in range(self.T)}
        beta_kj = {'beta_{}{}'.format(k, j): self.opt_model.binary_var(name='beta_{}{}.format(k, j)') for k in
                   range(self.V) for j in range(self.T)}
        return alpah_kj, beta_kj

    def generate_spatial_constraints(self, x_kl, y_kl):
        """
        :return: 这三个约束保证了船舶停靠泊位之间不会发生直接的交叉重合
        """
        for k in range(self.V):
            for l in range(k, self.V):
                self.opt_model.add_constraint(ct=(y_kl[(l, k)] + y_kl[(k, l)] <= 1), ctname='ts_cons3{}{}'.format(k, l))
                self.opt_model.add_constraint(ct=(x_kl[(l, k)] + x_kl[(k, l)] <= 1), ctname='ts_cons2{}{}'.format(k, l))
                self.opt_model.add_constraint(ct=(x_kl[(l, k)] + x_kl[(k, l)] + y_kl[(l, k)] + y_kl[(k, l)] >= 1),
                                              ctname='ts_cons1{}{}'.format(k, l))

    def generate_discrete_position_cons(self, pi_kn, theta_kn, H, ykl):
        """
        H：各个船的长度集合，输入参数
        :return: 这些约束使得船舶可以在空间层面有序停靠
        """
        # self.optmodel.sum([n * pi_kn[k , n] for n in range(self.B)]) == bk
        for k in range(self.V):
            self.opt_model.add_constraint(self.opt_model.sum([theta_kn[k, n]] for n in range(self.B)) == H[k],
                                          ctname='s_cons2{}'.format(k))
            self.opt_model.add_constraint(self.opt_model.sum([pi_kn[k, n] for n in range(self.B)] == 1),
                                          ctname='s_cons3{}'.format(k))
            self.opt_model.add_constraint(pi_kn[k, 0] >= theta_kn[k, 0], ctname='s_cons5{}'.format((k)))
            for n in range(self.B):
                if n != 0:
                    self.opt_model.add_constraint(theta_kn[k, n] - theta_kn[k, n - 1] <= pi_kn[k, n],
                                                  ctname='s_cons4{}{}'.format(k, n))
                    self.opt_model.add_constraint(1 - theta_kn[k, n - 1] >= pi_kn[k, n],
                                                  ctname='s_cons7{}{}'.format(k, n))
                self.opt_model.add_constraint(pi_kn <= theta_kn, ctname='s_cons6{}{}'.format(k, n))
                for l in range(0, self.V):
                    if l != k:
                        # check this formulation in paper, it seems to be wrong
                        self.opt_model.add_constraint(
                            ykl[k, l] + self.opt_model.sum(pi_kn[k, m] for m in range(max(0, n - H[l] + 1), self.J)) +
                            pi_kn[l, n] <= 2)

    def generation_discrete_temporal_cons(self, alpha_kj, belta_kj, ck, z_gjk, x_kl):
        # self.optmodel.sum([j * alpha_kj[k , j] for j in range(self.T)]) == tk
        for k in range(self.V):
            self.opt_model.add_constraint(self.opt_model.sum([alpha_kj[k, j] for j in range(self.T)]) == 1)
            self.opt_model.add_constraint(alpha_kj[k, 1] >= belta_kj[k, 1])
            for j in range(self.T):
                self.opt_model.add_constraint((j + 1) * belta_kj[k, j] <= ck[k])
                self.opt_model.add_constraint(alpha_kj[k, j] <= belta_kj[k, j])
                if j != 0:
                    self.opt_model.add_constraint(alpha_kj[k, j] >= belta_kj[k, j] - belta_kj[k, j-1])
                    self.opt_model.add_constraint(alpha_kj[k, j] <= 1 - belta_kj[k, j-1])
                for g in range(self.G):
                    self.opt_model.add_constraint(z_gjk[g, j, k] <= belta_kj[k, j])
                for l in range(self.V):
                    for i in range(self.T):
                        if k != l and i >= j - 1:
                            # recheck this constraint in the paper
                            self.opt_model.add_constraint(x_kl[k, l] + belta_kj[k, j] + belta_kj[l, j] <= 2)

    def generate_quay_crane_cons(self, z_gkj, alpha_kj):
        """
        use binary_val to replace the integer_val
        :param z_gkj:
        :param alpha_kj:
        :return:
        """
        for g in range(self.G):
            for j in range(self.T):
                self.opt_model.add_constraint(self.opt_model.sum([z_gkj[g, k, j] for k in range(self.V)]) <=1)
        for g in range(self.G):
            for k in range(self.V):
                serve_time = self.opt_model.sum([j * alpha_kj[k, j] for j in range(self.T)])
                for j in range(self.J):
                    self.opt_model.add_constraint(serve_time <= j * z_gkj[g, k, j] + (1 - z_gkj[g, k, j]))

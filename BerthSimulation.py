import random
import docplex.mp.model as cpx


class Vessel(object):
    def __init__(self, nums):
        self.V = nums
        self.vessels = {'vessel_{}'.format(i): {} for i in range(nums)}

    def generate_arrival_time(self, arrival_rate):
        # automated generate the arrival time by bosong distribution
        arrival_time = []
        for i in range(self.V):
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
        # automated generate the container volume of arrival vessel
        pass


class QuayCrane(object):
    def __init__(self, nums, Quay_length, *args):
        super(QuayCrane, self).__init__(*args)
        self.G = nums
        self.qc = {'qc_{}'.format(g): {} for g in range(nums)}
        self.B = Quay_length

    def generate_process_rate(self, process_rate_distribution):
        rand = len(process_rate_distribution)
        for key in self.qc:
            self.qc[key]['processing_rate'] = process_rate_distribution[random.randint(0, rand)]

    def generate_action_range(self):
        i = self.B / (2 * self.G)
        for key in self.qc:
            self.qc[key]['action_range'] = (i - self.B / (2 * self.G), i + self.B / self.G)
            i = i + self.B / self.G


class Berth(object, QuayCrane, Vessel):
    """
    Reproduction of the paper 'MIP approaches for the integrated berth allocation and quay crane assignment and
    scheduling problem'. In this paper, two different modeling method are proposed (i.e., a discrete one and a
    continues one).
    """

    def __init__(self, berth_length, ship_nums, time_horizon, safety_time, crane_nums, pattern='discrete'):
        """
        :param berth_length: The numbers of berth section: 1unit/setion = 10 meter
        :param ship_nums: Total numbers of vessels
        :param time_horizon: Total time_horizon for the berth allocation ana quay crane assignment
        :param safety_time: The safe time interval between departure vessels
        :param crane_nums:
        :param pattern: discrete or continue.
        """
        super(Berth, self).__init__(crane_nums, berth_length, ship_nums)
        self.opt_model = cpx.Model(name='BerthAllocation')
        self.T = time_horizon
        self.F = safety_time
        self.pattern = pattern

    def generate_variables(self):
        x_kl = {(k, l): self.opt_model.binary_var('x_{}{}'.format(k, l)) for k in range(self.V) for l in range(self.V)}
        y_kl = {(k, l): self.opt_model.binary_var('y_{}{}'.format(k, l)) for k in range(self.V) for l in range(self.V)}
        ck = {(k,): self.opt_model.integer_var(name='c_{}'.format(k), lb=0) for k in range(self.V)}
        z_gk_j = {(g, k, j): self.opt_model.binary_var('z_{}{}{}'.format(g, k, j)) for g in range(self.G) for k in
                  range(self.V) for j in range(self.T)}
        if self.pattern == 'continue':
            bk = {(k,): self.opt_model.integer_var(name='b_{}'.format(k)) for k in range(self.V)}
            tk = {(k,): self.opt_model.integer_var(name='t_{}'.format(k)) for k in range(self.V)}
            return x_kl, y_kl, bk, tk, ck, z_gk_j
        else:
            alpha_kj = {(k, j): self.opt_model.binary_var(name='alpha_{}{}'.format(k, j)) for k in range(self.V) for j
                        in range(self.T)}
            beta_kj = {(k, j): self.opt_model.binary_var(name='beta_{}{}'.format(k, j)) for k in range(self.V) for j
                       in range(self.T)}
            pi_kn = {(k, n): self.opt_model.binary_var(name='pi_{}{}'.format(k, n)) for k in range(self.V) for n
                     in range(self.B)}
            theta_kn = {(k, n): self.opt_model.binary_var(name='theta_{}{}'.format(k, n)) for k in range(self.V) for n
                        in range(self.B)}
            return x_kl, y_kl, pi_kn, theta_kn, alpha_kj, beta_kj, ck, z_gk_j

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
        # self.opt_model.sum([n * pi_kn[k , n] for n in range(self.B)]) == bk
        for k in range(self.V):
            vessel = 'vessel_{}'.format(k)
            self.opt_model.add_constraint(self.opt_model.sum([theta_kn[k, n]] for n in range(self.B)) ==
                                          H[vessel]['vessel_length'], ctname='s_cons2{}'.format(k))
            self.opt_model.add_constraint(self.opt_model.sum([pi_kn[k, n] for n in range(self.B)] == 1),
                                          ctname='s_cons3{}'.format(k))
            self.opt_model.add_constraint(pi_kn[k, 0] >= theta_kn[k, 0], ctname='s_cons5{}'.format(k))
            # constraint 7
            self.opt_model.add_constraint(self.opt_model.sum([n * pi_kn[k, n] for n in range(self.B)]) <= self.B - H[k])
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
                            ykl[k, l] + self.opt_model.sum(pi_kn[k, m] for m in range(max(0, n - H[l] + 1), self.B)) +
                            pi_kn[l, n] <= 2)

    def generation_discrete_temporal_cons(self, alpha_kj, beta_kj, ck, z_gjk, x_kl, A):
        # self.opt_model.sum([j * alpha_kj[k , j] for j in range(self.T)]) == tk
        """
        :param alpha_kj:
        :param beta_kj:
        :param ck:
        :param z_gjk:
        :param x_kl:
        :param A: arrival_time of vessel k
        :return:
        """
        for k in range(self.V):
            vessel = 'vessel_{}'.format(k)
            self.opt_model.add_constraint(self.opt_model.sum([alpha_kj[k, j] for j in range(self.T)]) == 1)
            self.opt_model.add_constraint(alpha_kj[k, 1] >= beta_kj[k, 1])
            self.opt_model.add_constraint(
                self.opt_model.sum([j * alpha_kj[k, j] for j in range(self.T)]) >= A[vessel]['arrival_time'])
            for j in range(self.T):
                self.opt_model.add_constraint((j + 1) * beta_kj[k, j] <= ck[k])
                self.opt_model.add_constraint(alpha_kj[k, j] <= beta_kj[k, j])
                if j != 0:
                    self.opt_model.add_constraint(alpha_kj[k, j] >= beta_kj[k, j] - beta_kj[k, j - 1])
                    self.opt_model.add_constraint(alpha_kj[k, j] <= 1 - beta_kj[k, j - 1])
                for g in range(self.G):
                    self.opt_model.add_constraint(z_gjk[g, j, k] <= beta_kj[k, j])
                for l in range(self.V):
                    for i in range(self.T):
                        if k != l and i >= j - 1:
                            # recheck this constraint in the paper
                            self.opt_model.add_constraint(x_kl[k, l] + beta_kj[k, j] + beta_kj[l, j] <= 2)

    def generate_quay_crane_cons(self, z_gkj, alpha_kj, ck, pi_kn, y_kl, qc, vessels):
        """
        use binary_val to replace the integer_val
        P: processing rate of each quay crane
        H：length of vessel k
        E: upper action range of crane g
        S: lower action range of crane g
        :param vessels:
        :param qc:
        :param y_kl:
        :param pi_kn:
        :param ck:
        :param z_gkj:
        :param alpha_kj:
        :return:
        """
        for g in range(self.G):
            for j in range(self.T):
                self.opt_model.add_constraint(self.opt_model.sum([z_gkj[g, k, j] for k in range(self.V)]) <= 1)

        for k in range(self.V):
            vessel = 'vessel_{}'.format(k)
            self.opt_model.add_constraint(
                self.opt_model.sum(
                    [self.opt_model.sum([qc['qc_{}'.format(g)]['processing_rate'] * z_gkj[g, j, k] for g in self.G]) for
                     j in self.T])
                >= vessels[vessel]['container_volume'])

        for g in range(self.G):
            qc_index = 'qc_{}'.format(g)
            for g_hat in range(g):
                for k in range(self.V):
                    for l in range(self.V):
                        for j in range(self.T):
                            self.opt_model.add_constraint(z_gkj[g, k, j] + z_gkj[g_hat, l, j] <= 2 - y_kl[k, l])

            for k in range(self.V):
                vessel = 'vessel_{}'.format(k)
                serve_time = self.opt_model.sum([j * alpha_kj[k, j] for j in range(self.T)])
                berth_location = self.opt_model.sum([n * pi_kn[k, n] for n in range(self.B)])
                for j in range(self.T):
                    self.opt_model.add_constraint(serve_time <= j * z_gkj[g, k, j] + (1 - z_gkj[g, k, j]))
                    self.opt_model.add_constraint(ck >= (j + 1) * z_gkj[g, k, j])
                    self.opt_model.add_constraint(
                        berth_location + vessels[vessel]['vessel_length'] <= qc[qc_index]['action_range'][1] * z_gkj[
                            g, k, j] + (1 - z_gkj[g, k, j]) * self.B)
                    self.opt_model.add_constraint(berth_location >= qc[qc_index]['action_range'][0] * z_gkj[g, k, j])

    def model(self):
        if self.pattern != 'discrete':
            x_kl, y_kl, pi_kn, theta_kn, alpha_kj, beta_kj, ck, z_gk_j = self.generate_variables()
            self.generate_discrete_position_cons(pi_kn, theta_kn, self.vessels, y_kl)
            self.generate_spatial_constraints(x_kl, y_kl)
            self.generation_discrete_temporal_cons(alpha_kj, beta_kj, ck, z_gk_j, x_kl, self.vessels)
            self.generate_quay_crane_cons(z_gk_j, alpha_kj, ck, pi_kn, y_kl, self.qc, self.vessels)
        else:
            x_kl, y_kl, bk, tk, ck, z_gk_j = self.generate_variables()
        obj = self.opt_model.sum([ck[k] for k in range(self.V)])
        self.opt_model.minimize(obj)
        self.opt_model.solve()

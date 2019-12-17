import numpy as np
import random
from aminoagent import AminoAgent
from rotate import rotate

class AminoModel():

    def __init__(self, aa_lst, pos = None):
        self.num_aa = len(aa_lst)
        self.schedule = {}
        self.grid = {}
        for i in range(2 * self.num_aa):
              for j in range(2 * self.num_aa):
                  for k in range(2 * self.num_aa):
                      self.grid[(i, j, k)] = 0
        # create aa
        for i, aa in enumerate(aa_lst):
            a = AminoAgent((self.num_aa + i, self.num_aa, self.num_aa), (self.num_aa + i, self.num_aa, self.num_aa + 1), self, i, aa)
            self.grid[(self.num_aa + i, self.num_aa, self.num_aa)] = i + 1
            self.grid[(self.num_aa + i, self.num_aa, self.num_aa + 1)] = i + 1
            self.schedule[i] = a
        self.total_energy = 0
    def find_open_spots(self, aa_num):
        open_spots = []
        rotations = [(-1, 1, 0), (1, -1, 0), (-1, 0, 1), (1, 0, -1), (0, -1, 1), (0, 1, -1)]
        for rot in rotations:
            new_pos = []
            conflicts = False
            i = 0
            pivot = self.schedule[aa_num - 1].pos_ca
            # print(pivot, aa_num)
            while not conflicts and i < len(self.schedule) - aa_num:
                curr_pos1 = self.schedule[aa_num + i].pos_ca
                curr_pos2 = self.schedule[aa_num + i].pos_side
                nxt_pos1 = rotate(pivot, curr_pos1, rot)
                nxt_pos2 = rotate(pivot, curr_pos2, rot)
                # print(pivot, curr_pos, nxt_pos)
                if (self.grid[nxt_pos1] == 0 or self.grid[nxt_pos1] >= aa_num + 1) and (self.grid[nxt_pos2] == 0 or self.grid[nxt_pos2] >= aa_num + 1):
                    new_pos.append((nxt_pos1, nxt_pos2))
                else:
                    conflicts = True
                i += 1
            if len(new_pos) >= len(self.schedule) - aa_num:
                open_spots.append(new_pos)
        return open_spots
            
    def step(self, update=False):
        '''Advance the model by one step.'''
        moved = False
        while not moved:
            aa_chosen = random.randint(2, self.num_aa - 1)
            open_spots = self.find_open_spots(aa_chosen)
            # seq = [x.pos_ca for x in self.schedule.values()]
            # print(seq)
            if len(open_spots) > 0:
                chosen_rot = random.randint(0, len(open_spots)-1)
                for i, x in enumerate(open_spots[chosen_rot]):
                    newpos1, newpos2 = x
                    pos1 = self.schedule[aa_chosen + i].pos_ca
                    pos2 = self.schedule[aa_chosen + i].pos_side
                    if self.grid[pos1] == aa_chosen + i + 1:
                        self.grid[pos1] = 0
                    if self.grid[pos2] == aa_chosen + i + 1:
                        self.grid[pos2] = 0
                    self.grid[newpos1] = aa_chosen + i + 1
                    self.grid[newpos2] = aa_chosen + i + 1
                    self.schedule[aa_chosen + i].set_pos(newpos1, newpos2)
                moved = True
        if update:
            for i in self.schedule:
                self.schedule[i].update()
            total_energy = sum([self.schedule[x].energy for x in self.schedule])
            self.total_energy = total_energy
            return total_energy
        else:
            return -1


class AminoModelWeighted():

    def __init__(self, aa_lst, pos = None):
        self.num_aa = len(aa_lst)
        self.schedule = {}
        self.grid = {}
        for i in range(2 * self.num_aa):
              for j in range(2 * self.num_aa):
                  for k in range(2 * self.num_aa):
                      self.grid[(i, j, k)] = 0
        # create aa
        for i, aa in enumerate(aa_lst):
            a = AminoAgent((self.num_aa + i, self.num_aa, self.num_aa), (self.num_aa + i, self.num_aa, self.num_aa + 1), self, i, aa)
            self.grid[(self.num_aa + i, self.num_aa, self.num_aa)] = i + 1
            self.grid[(self.num_aa + i, self.num_aa, self.num_aa + 1)] = i + 1
            self.schedule[i] = a
        self.energies = [1/(self.num_aa-1) for i in range(self.num_aa-1)]
    def find_open_spots(self, aa_num):
        open_spots = []
        rotations = [(-1, 1, 0), (1, -1, 0), (-1, 0, 1), (1, 0, -1), (0, -1, 1), (0, 1, -1)]
        for rot in rotations:
            new_pos = []
            conflicts = False
            i = 0
            pivot = self.schedule[aa_num - 1].pos_ca
            # print(pivot, aa_num)
            while not conflicts and i < len(self.schedule) - aa_num:
                curr_pos1 = self.schedule[aa_num + i].pos_ca
                curr_pos2 = self.schedule[aa_num + i].pos_side
                nxt_pos1 = rotate(pivot, curr_pos1, rot)
                nxt_pos2 = rotate(pivot, curr_pos2, rot)
                # print(pivot, curr_pos, nxt_pos)
                if (self.grid[nxt_pos1] == 0 or self.grid[nxt_pos1] >= aa_num + 1) and (self.grid[nxt_pos2] == 0 or self.grid[nxt_pos2] >= aa_num + 1):
                    new_pos.append((nxt_pos1, nxt_pos2))
                else:
                    conflicts = True
                i += 1
            if len(new_pos) >= len(self.schedule) - aa_num:
                open_spots.append(new_pos)
        return open_spots
            
    def step(self, update=False):
        '''Advance the model by one step.'''
        moved = False
        while not moved:
            aa_chosen = random.choices([i + 1 for i in range(self.num_aa - 1)], weights=self.energies, k=1)[0]
            open_spots = self.find_open_spots(aa_chosen)
            # seq = [x.pos_ca for x in self.schedule.values()]
            # print(seq)
            if len(open_spots) > 0:
                chosen_rot = random.randint(0, len(open_spots)-1)
                for i, x in enumerate(open_spots[chosen_rot]):
                    newpos1, newpos2 = x
                    pos1 = self.schedule[aa_chosen + i].pos_ca
                    pos2 = self.schedule[aa_chosen + i].pos_side
                    if self.grid[pos1] == aa_chosen + i + 1:
                        self.grid[pos1] = 0
                    if self.grid[pos2] == aa_chosen + i + 1:
                        self.grid[pos2] = 0
                    self.grid[newpos1] = aa_chosen + i + 1
                    self.grid[newpos2] = aa_chosen + i + 1
                    self.schedule[aa_chosen + i].set_pos(newpos1, newpos2)
                moved = True
        if update:
            for i in self.schedule:
                self.schedule[i].update()
            energies = np.array([self.schedule[x].energy for x in self.schedule])[1:]
            total_energy = sum(energies)
            min_energy = min(energies)
            self.energies = energies - min_energy
            return total_energy
        else:
            return -1

class MultiAminoModelWeighted():

    def __init__(self, lst_of_models, opp=False):
        self.num_aa = [x.num_aa for x in lst_of_models]
        self.schedule = {}
        self.grid = {}
        self.max_aa = max(self.num_aa)
        for i in range(4 * self.max_aa):
              for j in range(4 * self.max_aa):
                  for k in range(4 * self.max_aa):
                      self.grid[(i, j, k)] = (0, 0)
        # create aa
        for k in range(len(lst_of_models)):
            aa_lst = [x for x in lst_of_models[k].schedule]
            if opp:
                for i, aa in enumerate(aa_lst):
                    x1, y1, z1 = lst_of_models[k].schedule[aa].pos_ca
                    x2, y2, z2 = lst_of_models[k].schedule[aa].pos_side
                    if i % 2 == 0:
                        move = len(aa_lst)
                    else:
                        move = 0
                        piv = (len(aa_lst), len(aa_lst), len(aa_lst))
                        pt = (x1, y1, z1)
                        x1, y1, z1 = rotate(piv, rotate(piv, rotate(piv, pt, (0, 1, -1)), (-1, 0, 1)), (-1, 1, 0))
                        pt = (x2, y2, z2)
                        x2, y2, z2 = rotate(piv, rotate(piv, rotate(piv, pt, (0, 1, -1)), (-1, 0, 1)), (-1, 1, 0))
                    pt1 = (x1 + move, y1 + move, z1 + move)
                    pt2 = (x2 + move, y2 + move, z2 + move)
                    a = AminoAgent(pt1, pt2, self, i, aa)
                    self.grid[pt1] = (i + 1, k)
                    self.grid[pt2] = (i + 1, k)
                    self.schedule[(i, k)] = a
            else:
                for i, aa in enumerate(aa_lst):
                    x1, y1, z1 = lst_of_models[k].schedule[aa].pos_ca
                    x2, y2, z2 = lst_of_models[k].schedule[aa].pos_side
                    if i % 2 == 0:
                        move = len(aa_lst)
                    else:
                        move = 0
                        piv = (len(aa_lst), len(aa_lst), len(aa_lst))
                        pt = (x1, y1, z1)
                        x1, y1, z1 = rotate(piv, rotate(piv, rotate(piv, pt, (0, 1, -1)), (-1, 0, 1)), (-1, 1, 0))
                        pt = (x2, y2, z2)
                        x2, y2, z2 = rotate(piv, rotate(piv, rotate(piv, pt, (0, 1, -1)), (-1, 0, 1)), (-1, 1, 0))
                    pt1 = (x1 + move, y1 + move, z1 + move)
                    pt2 = (x2 + move, y2 + move, z2 + move)
                    a = AminoAgent(pt1, pt2, self, i, aa)
                    self.grid[pt1] = (i + 1, k)
                    self.grid[pt2] = (i + 1, k)
                    self.schedule[(i, k)] = a
        self.energies = [[1/(x-1) for i in range(x-1)] for x in self.num_aa]
    def find_open_spots(self, aa_num, aa_id):
        open_spots = []
        rotations = [(-1, 1, 0), (1, -1, 0), (-1, 0, 1), (1, 0, -1), (0, -1, 1), (0, 1, -1)]
        for rot in rotations:
            new_pos = []
            conflicts = False
            i = 0
            pivot = self.schedule[(aa_num - 1, aa_id)].pos_ca
            # print(pivot, aa_num)
            while not conflicts and i < self.num_aa[aa_id] - aa_num:
                curr_pos1 = self.schedule[(aa_num + i, aa_id)].pos_ca
                curr_pos2 = self.schedule[(aa_num + i, aa_id)].pos_side
                nxt_pos1 = self.rotate(pivot, curr_pos1, rot)
                nxt_pos2 = self.rotate(pivot, curr_pos2, rot)
                if all([x >=0 and x < self.max_aa for x in nxt_pos1]) and all([x >=0 and x < self.max_aa for x in nxt_pos2]):
                    if (self.grid[nxt_pos1][0] == 0 or (self.grid[nxt_pos1][0] >= aa_num + 1 and self.grid[nxt_pos1][1] >= aa_id)) and (self.grid[nxt_pos2][0] == 0 or (self.grid[nxt_pos2][0] >= aa_num + 1 and self.grid[nxt_pos1][1] == aa_id)):
                        new_pos.append((nxt_pos1, nxt_pos2))
                    else:
                        conflicts = True
                else:
                    conflicts = True
                i += 1
            if len(new_pos) >= len(self.schedule) - aa_num:
                open_spots.append(new_pos)
        return open_spots
            
    def step(self, update=False):
        '''Advance the model by one step.'''
        moved = False
        while not moved:
            model_chosen = random.randint(0, len(self.num_aa) - 1)
            aa_chosen = random.choices([i + 1 for i in range(self.num_aa[model_chosen] - 1)], weights=self.energies[model_chosen], k=1)[0]
            open_spots = self.find_open_spots(aa_chosen, model_chosen)
            # seq = [x.pos_ca for x in self.schedule.values()]
            # print(seq)
            if len(open_spots) > 0:
                chosen_rot = random.randint(0, len(open_spots)-1)
                for i, x in enumerate(open_spots[chosen_rot]):
                    newpos1, newpos2 = x
                    pos1 = self.schedule[(aa_chosen + i, model_chosen)].pos_ca
                    pos2 = self.schedule[(aa_chosen + i, model_chosen)].pos_side
                    if self.grid[pos1][0] == aa_chosen + i + 1 and self.grid[pos1][1] == model_chosen:
                        self.grid[pos1] = (0, 0)
                    if self.grid[pos2][0] == aa_chosen + i + 1 and self.grid[pos1][1] == model_chosen:
                        self.grid[pos2] = (0, 0)
                    self.grid[newpos1] = (aa_chosen + i + 1, model_chosen)
                    self.grid[newpos2] = (aa_chosen + i + 1, model_chosen)
                    self.schedule[(aa_chosen + i, model_chosen)].set_pos(newpos1, newpos2)
                moved = True
        if update:
            for i in self.schedule:
                self.schedule[i].update()
            energies = [[i.energy for i in x][1:] for x in self.schedule]
            total_energy = sum(sum(energies))
            min_energy = [min(x) for x in energies]
            self.energies = [(np.array(energies[i]) - min_energy[i]) for i in range(len(self.num_aa))]
            return total_energy
        else:
            return -1
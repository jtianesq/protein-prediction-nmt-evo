from energy import freeEnergy, get_energy
from aminomodel import AminoModel, AminoModelWeighted, MultiAminoModelWeighted

class AminoAgent():
  
    def __init__(self, pos1, pos2, model, unique_id, agent_type):
        self.pos_ca = pos1
        self.pos_side = pos2
        self.a_type = agent_type
        self.energy = np.inf
        self.unique_id = unique_id
        self.model = model

    def set_pos(self, coord1, coord2):
        self.pos_ca = coord1
        self.pos_side = coord2
    
    def update(self):
        energy = 0
        for i in self.model.schedule:
            other = self.model.schedule[i]
            if other.unique_id != self.unique_id:
                energy += freeEnergy([self.pos_ca, get_energy(self.a_type), self.a_type], [other.pos_ca, get_energy(other.a_type), other.a_type])
        self.energy = energy

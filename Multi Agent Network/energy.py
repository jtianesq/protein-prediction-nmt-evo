

# Based on mccammon1984 and wallqvist1994

EnergyTable = {}

def energyInput(energyRadius):
  EnergyTable[energyRadius[0]] = (energyRadius[1], energyRadius[2])

def get_energy(a):
  if a in EnergyTable:
    return EnergyTable[a]
  else:
    raise ValueError("{} not found in energy table.".format(a))


energyInput(['G', 0, 2.18])
energyInput(['A', 1.76, 1.77])
energyInput(['V', 6.95, 2.07])
energyInput(['P', 3.08, 2.13])
energyInput(['T', 1.46, 2.05])
energyInput(['S', -0.21, 1.64])
energyInput(['N', -3.43, 2.20])
energyInput(['D', -4.39, 2.05])
energyInput(['R', -1.70, 2.17])
energyInput(['K', -1.66, 2.49])
energyInput(['E', -1.37, 2.00])
energyInput(['Q', -0.18, 2.39])
energyInput(['L', 4.85, 2.31])
energyInput(['I', 5.15, 2.81])
energyInput(['F', 5.11, 2.72])
energyInput(['Y', 2.74, 2.81])
energyInput(['W', 8.96, 2.90])
energyInput(['M', 3.52, 2.41])
energyInput(['C', 2.81, 2.28])
energyInput(['H', 0.83, 2.28])

def dist(a,b):
  x = a[0]**2 + b[0]**2
  y = a[1]**2 + b[1]**2
  z = a[2]**2 + b[2]**2
  return math.sqrt(x+y+z)

# each amino acid a is a list defined by a = [(x,y,z), (g,r)]
def freeEnergy(a,b):
  interaction = 1
  lst1 = ['R', 'K', 'D', 'E']
  lst2 = ['Q', 'N', 'H', 'S', 'T', 'Y', 'C']
  if a[2] in lst1 and b[2] in lst1:
      interaction = 0.5
  elif a[2] in lst2 and b[2] in lst2:
      interaction = 0.8
  deltaG = a[1][0] + b[1][0]
  modRad = a[1][1] + b[1][1]
  rad = dist(a[0],b[0])*modRad
  
  return interaction * deltaG / rad
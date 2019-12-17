from plots import plot_3d_model, plot_3d_model_multi
from aminomodel import AminoModel, AminoModelWeighted, MultiAminoModelWeighted
import random

  
def multifold(primary_seq, prob, k, total_steps, epoch_intervals, plot=False, anneal=False, model='unweighted', begin=None, opp=False):
    if model == 'multi':
      model = MultiAminoModelWeighted(begin, opp=opp)
    elif begin:
      model = begin
    elif model == 'unweighted':
      model = AminoModel(primary_seq)
    else:
      model = AminoModelWeighted(primary_seq)
    topk_models = []
    energies = []
    for i in range(total_steps):
        if i % epoch_intervals == 0:
          global_energy = model.step(update=True)
          energies.append(global_energy)
          topk_models.append((global_energy, model))
          if i / epoch_intervals > k:
              topk_models.sort()
              topk_models = topk_models[:k]
        else:
            global_energy = model.step()
        if i % epoch_intervals == 0 and i != 0:
            # simulate transition
            p = random.random()
            if anneal:
                prob_to_beat = np.sqrt(i/total_steps)
            else:
                prob_to_beat = prob
            if p < prob_to_beat:
                model = topk_models[0][1]
            else:
                rint = random.randint(1, len(topk_models))
                if rint != len(topk_models):
                  model = topk_models[rint][1]
    final_energy, model = topk_models[0]
    if plot:
        if model != "multi":
            plot_3d_model(model)
        else:
            plot_3d_model_multi(model)
    return energies, topk_models
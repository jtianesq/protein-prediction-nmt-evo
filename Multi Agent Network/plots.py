import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def check_valid(seq1, seq2):
    prev_points = set([])
    for i in range(len(seq1) - 1):
        pair = seq1[i: i + 2]
        if pair[1] in prev_points:
            return False
        else:
            if np.sum(np.abs(np.array(pair[1]) - np.array(pair[0]))) != 1:
                return False
            prev_points.add(pair[1])
    for i in range(len(seq1) - 1):
        pair = [seq1[i]] + [seq2[i]]
        if pair[1] in prev_points:
            return False
        else:
            if np.sum(np.abs(np.array(pair[1]) - np.array(pair[0]))) != 1:
                return False
            prev_points.add(pair[1])
    return True

def plot_3d_model(model):
    pts = np.array([model.schedule[x].pos_ca for x in model.schedule])
    pts_side = np.array([model.schedule[x].pos_side for x in model.schedule])
    labels = [model.schedule[x].a_type for x in model.schedule]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for label in np.unique(labels):
        ix = np.where(np.array(labels) == label)
        xs = pts[ix][:,0]
        ys = pts[ix][:,1]
        zs = pts[ix][:,2]
        ax.scatter3D(xs, ys, zs, label = label)
    
    for i in range(len(pts) - 1):
        pair = pts[i:i + 2]
        xs = pair[:, 0]
        ys = pair[:, 1]
        zs = pair[:, 2]
        ax.plot3D(xs, ys, zs, 'gray')
    
    for i in range(len(pts) - 1):
        xs = [pts[i][0]] + [pts_side[i][0]]
        ys = [pts[i][1]] + [pts_side[i][1]]
        zs = [pts[i][2]] + [pts_side[i][2]]
        ax.plot3D(xs, ys, zs, 'gray')

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.axis('equal')
    ax.legend(loc='center left', bbox_to_anchor=(0.9, 0.5))
    figure = plt.figure()
    plt.show()   
  
def plot_3d_model_multi(model):
    for i, lst in enumerate(model.num_aa):
        pts = np.array([model.schedule[(x, i)].pos_ca for x in lst])
        pts_side = np.array([model.schedule[(x, i)].pos_side for x in lst])
        labels = [model.schedule[(x, i)].a_type for x in lst]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for label in np.unique(labels):
            ix = np.where(np.array(labels) == label)
            xs = pts[ix][:,0]
            ys = pts[ix][:,1]
            zs = pts[ix][:,2]
            ax.scatter3D(xs, ys, zs, label = label)

        for i in range(len(pts) - 1):
            pair = pts[i:i + 2]
            xs = pair[:, 0]
            ys = pair[:, 1]
            zs = pair[:, 2]
            ax.plot3D(xs, ys, zs, 'gray')

        for i in range(len(pts) - 1):
            xs = [pts[i][0]] + [pts_side[i][0]]
            ys = [pts[i][1]] + [pts_side[i][1]]
            zs = [pts[i][2]] + [pts_side[i][2]]
            ax.plot3D(xs, ys, zs, 'gray')

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.axis('equal')
    ax.legend(loc='center left', bbox_to_anchor=(0.9, 0.5))
    figure = plt.figure()
    plt.show()
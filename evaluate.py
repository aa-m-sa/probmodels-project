
import sys
import pickle

import matplotlib.pyplot as plt
#import seaborn as sns
import math


def evaluate_arcs(dg, arcs):
    n = len(dg)
    n_total_arcs = n * n - n

    dg_arcs = set()
    for node, parents in dg:
        for parent in parents:
            dg_arcs.add((parent, node))

    n_present_arcs = len(dg_arcs)
    n_absent_arcs = n_total_arcs - n_present_arcs

    n_correct = 0
    n_incorrect = 0
    auc = 0
    xs = [0]
    ys = [0]

    for arc in arcs:
        if arc in dg_arcs:
            n_correct += 1
            xs.append(xs[-1])
            ys.append(n_correct / n_present_arcs)
        else:
            n_incorrect += 1
            auc += (n_correct / n_present_arcs) * (1 / n_absent_arcs)
            xs.append(n_incorrect / n_absent_arcs)
            ys.append(ys[-1])

    xs.append(1)
    ys.append(1)

    return auc, xs, ys


def probability_in_network(bn, row):
    dag, cpts = bn

    row_prob = 1
    for node, parents in dag:
        node_value = row[node]
        parent_values = tuple(row[p] for p in parents)
        row_prob *= cpts[node][parent_values][node_value]

    return row_prob


def evaluate_probs(bn, test_data, probs):
    def kullback_leibler_divergence(p, q):
        return sum(p[i] * math.log(p[i] / q[i]) for i in range(len(p)))

    def normalize(probs):
        normalizer = sum(probs)
        return [p / normalizer for p in probs]

    true_probs = [probability_in_network(bn, row) for row in test_data]
    true_probs = normalize(true_probs)

    return kullback_leibler_divergence(true_probs, probs), probs, true_probs


def draw_roc(auc, xs, ys):
    t = 0.03

    plt.xlim([-t, 1+t])
    plt.ylim([-t, 1+t])
    plt.gca().set_aspect('equal', adjustable='box')

    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')
    plt.xlabel("False positives", fontsize=22, labelpad=16)
    plt.ylabel("True positives", fontsize=22, labelpad=16)

    plt.tick_params(axis='both', which='major', labelsize=18, length=8, direction="out", pad=10)

    plt.plot([-t,1+t], [-t,1+t], color="#888888")

    plt.plot(xs, ys, lw=2, color="k")

    plt.title("AUC = %.3f" % auc, fontsize=24)

    plt.tight_layout(rect=[0, 0, 1, 1])

    #plt.savefig("roc.pdf", bbox_inches='tight', transparent=True)
    plt.show()
    plt.close()


def draw_scatter(kld, probs, true_probs):
    plt.title("KL = %.3f" % kld, fontsize=26)

    plt.xscale('log')
    plt.yscale('log')

    plt.xlim((min(true_probs), max(true_probs)))

    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')
    plt.xlabel("True", fontsize=22, labelpad=16)
    plt.ylabel("Predicted", fontsize=22, labelpad=16)

    plt.tick_params(axis='both', which='major', labelsize=18, length=8, direction="out", pad=10)

    plt.plot([min(true_probs), max(true_probs)], [min(true_probs), max(true_probs)], "r", lw=3)

    plt.plot(true_probs, probs, 'o', ms=5, markerfacecolor="None", markeredgecolor='k', markeredgewidth=1)

    #plt.savefig("scatter.pdf", bbox_inches='tight', transparent=True)
    plt.show()
    plt.close()

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("Usage: ./evaluate.py fulldata.pickle arcs.txt probs.txt")
        exit()


    with open(sys.argv[1], "rb") as f:
        fulldata = pickle.load(f)

    with open(sys.argv[2], "r") as f:
        arcs = [tuple(a.strip().split(" ")) for a in f.readlines()]

    with open(sys.argv[3], "r") as f:
        probs = [float(p) for p in f.read().split()]


    network = fulldata["network"]
    test_data = fulldata["test-data"]

    auc, xs, ys = evaluate_arcs(network[0], arcs)
    kld, probs, true_probs = evaluate_probs(network, test_data, probs)

    print("AUC:", auc)
    print("KLD:", kld)

    draw_roc(auc, xs, ys)
    draw_scatter(kld, probs, true_probs)


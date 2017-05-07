import sys
import pickle

import matplotlib.pyplot as plt
import seaborn as sns
import math


from evaluate import evaluate_arcs

sns.set_style("ticks")
sns.set_style("whitegrid")
sns.set_palette('Paired')
#sns.set_palette(sns.dark_palette('palegreen'))


def draw_roc_all(aucs_all):

    t = 0.03

    #plt.gca().set_aspect('equal', adjustable='box')
    f = plt.figure()
    ax = f.add_subplot(111)

    #plt.gca().xaxis.set_ticks_position('bottom')
    #plt.gca().yaxis.set_ticks_position('left')
    #plt.xlabel("False positives", fontsize=22, labelpad=16)
    #plt.ylabel("True positives", fontsize=22, labelpad=16)
    plt.xlabel("False positives", )
    plt.ylabel("True positives", )

    #plt.tick_params(axis='both', which='major', labelsize=18, length=8, direction="out", pad=10)
    #plt.tick_params(axis='both', which='major',  direction="out",)

    plt.plot([-t,1+t], [-t,1+t], color="#888888")

    for (auc, xs, ys, legend) in aucs_all:



        #plt.plot(xs, ys, lw=2, color="k")
        ax.plot(xs, ys, lw=2, label="{0}, auc: {1:.5}".format(legend, auc))


    #plt.title("AUC = %.3f" % auc, fontsize=24)
    ax.legend()

    #ax.axis('equal')
    ax.set_ylim(ymin=-t, ymax=1+t)
    ax.set_xlim(xmin=-t, xmax=1+t)
    ax.set_aspect('equal')
    #plt.tight_layout(rect=[0, 0, 1, 1])
    plt.savefig("roc5.pdf")
    #plt.show()
    plt.close()



if __name__ == "__main__":

    with open(sys.argv[1], "rb") as f:
        fulldata = pickle.load(f)


    arc_files = ['../return-2017-04-05/aaro-salosensaari_arcs.txt',
                 '../return-2017-04-12/aaro-salosensaari_arcs.txt',
                 '../return-2017-04-26/aaro-salosensaari_arcs.txt',
                 '../return-2017-05-03/aaro-salosensaari_arcs.txt']

    legends = ['w2', 'w3', 'w4', 'w5']

    network = fulldata["network"]

    arcs_all = []
    aucs_all = []

    for arc_file, l in zip(arc_files, legends):
        with open(arc_file, "r") as f:
            arcs = [tuple(a.strip().split(" ")) for a in f.readlines()]
            arcs_all.append(arcs)
            auc, xs, ys = evaluate_arcs(network[0], arcs)
            aucs_all.append((auc, xs, ys, l))


    draw_roc_all(aucs_all)

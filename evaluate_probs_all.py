import sys
import pickle

import matplotlib.pyplot as plt
import seaborn as sns
import math
import numpy as np


from evaluate import evaluate_probs

sns.set_style("ticks")
sns.set_style("whitegrid")
sns.set_palette('Paired')


def draw_scatter_all(probs_all):

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

def draw_scatter(kld, probs, true_probs, l, ylim=None):
    plt.figure()
    plt.title('{0} KL = {1:.3}'.format(l,kld))

    plt.xscale('log')
    plt.yscale('log')

    plt.xlim((min(true_probs), max(true_probs)))

    if ylim:
        plt.ylim(ylim[0], ylim[1])

    #plt.gca().xaxis.set_ticks_position('bottom')
    #plt.gca().yaxis.set_ticks_position('left')
    #plt.xlabel("True", fontsize=22, labelpad=16)
    #plt.ylabel("Predicted", fontsize=22, labelpad=16)
    plt.xlabel('True')
    plt.ylabel('Predicted')

    #plt.tick_params(axis='both', which='major', labelsize=18, length=8, direction="out", pad=10)

    plt.plot([min(true_probs), max(true_probs)], [min(true_probs), max(true_probs)], "r", lw=3)

    plt.plot(true_probs, probs, 'o', ms=5, markerfacecolor="None", markeredgecolor='k', markeredgewidth=1)

    plt.savefig("scatter{0}.pdf".format(l))
    #plt.show()
    plt.close()


def draw_scatter_ax(kld, probs, true_probs, l, ylim=None, ax=None):

    ts = '{0} KL = {1:.3}'.format(l,kld)
    ax.set_title(ts)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim((min(true_probs), max(true_probs)))

    if ylim:
        ax.set_ylim(ylim[0], ylim[1])
    #ax.set_xlabel('True')
    ax.set_ylabel('Predicted')

    ax.plot(true_probs, probs, '.')

    ax.plot([min(true_probs), max(true_probs)], [min(true_probs), max(true_probs)], "r", lw=2)

    #ax.plot(true_probs, probs, 'o', ms=5, markerfacecolor="None", markeredgecolor='k', markeredgewidth=1)


if __name__ == "__main__":



    with open(sys.argv[1], "rb") as f:
        fulldata = pickle.load(f)



    #prb_files = ['../return-2017-04-05/aaro-salosensaari_probs.txt',
    #             '../return-2017-04-12/aaro-salosensaari_probs.txt',
    #             '../return-2017-04-26/aaro-salosensaari_probs.txt',
    #             '../return-2017-05-03/aaro-salosensaari_probs.txt']

    #legends = ['w2', 'w3', 'w4', 'w5']

    prb_files = ['../code/alpha_1_test.txt',
                 '../code/alpha_10_test.txt',
                 '../code/alpha_15_test.txt',
                 '../code/alpha_100_test.txt']
    legends = ['alpha = 1,', 'alpha = 10,', 'alpha = 15,', 'alpha = 100,']


    network = fulldata["network"]
    test_data = fulldata["test-data"]

    probs_all = []

    ymin = np.inf
    ymax = -np.inf
    for pf,l in zip(prb_files, legends):
        with open(pf, "r") as f:
            probs = [float(p) for p in f.read().split()]

            kld, probs, true_probs = evaluate_probs(network, test_data, probs)
            probs_all.append((kld,probs, true_probs, l))
            if min(probs) < ymin:
                ymin = min(probs)
            if max(probs) > ymax:
                ymax = max(probs)


    print(ymin, ymax)


    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    for (kld, probs, true_probs,l), ax  in zip(probs_all,[ax1, ax2, ax3, ax4]):
        draw_scatter_ax(kld, probs, true_probs,l,ylim=(ymin,ymax),ax=ax)


    ax3.set_xlabel('True')
    ax4.set_xlabel('True')

    plt.savefig('scatter_alpha.pdf')

    #draw_roc(auc, xs, ys)
    #draw_scatter_all(probs_all)


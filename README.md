# readme

Source codes for my Project in Probabilistic Models course submissions.

Includes sources for [BEANDisco](https://www.cs.helsinki.fi/u/tzniinim/BEANDisco/) by Teppo Niinimäki, which is GPL licensed.

## Project in Probabilistic Models?

### Problem

Project was a 'hands-free' challenge: Given a training data set generated from a unknown directed Bayesian (graphical) model, learn the best model (both the structure and the parameters). Some details:

The training data set consisted of 2500 instances like this:

```
A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
1 1 1 0 1 2 2 0 1 0 2 0 0 0 0 0 0 0 2 1 1 0 1 1 2 2
0 1 0 2 1 2 2 0 1 0 1 0 2 0 2 1 0 0 0 1 1 0 1 0 2 2
1 1 0 1 0 2 1 0 1 0 2 2 2 0 2 1 2 0 1 0 1 0 0 1 1 2
```

In other words, we had 26 variables with 3 discrete values. The task was to find a structural relationships between them.

Algorithms I wrote are based on approach of trying to find a network that has the best BDeu / AIC / BIC scores w.r.t. the training data.
I also tried algorithm implemented in the BEANDisco software.


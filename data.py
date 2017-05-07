
# The file fulldata.pickle is a serialized Python object than can be deserialized with the pickle module.

import pickle

with open("fulldata.pickle", "rb") as f:
    data = pickle.load(f)

# The data object is a dictionary containing all of course data.

print(data.keys())

# The most relevant one is the network, consisting of a DAG and the probability parameters.

dag, params = data["network"]

# The dag object is a list of nodes and their respective parents in a topologically sorted order.

for node, parents in dag:
    print(node, parents)

# The params object is a dictionary that maps each variable to its local probability table.
# For instance, the table of E is

table_E = params["E"]

# The table is a dictionary that maps each instantiation of the variable's parents to a probability distribution.
# For instance, the parents of E are A and V, and the local distribution of E when A=1 and V=2 is

distribution = table_E[(1,2)]

# The distribution is a 3-tuple of probabilities.

#print(distribution)


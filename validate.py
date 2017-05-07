
import sys

args = sys.argv[1:]

if len(args) < 2:
	print("Usage: python3 validate.py arcs.txt probs.txt")
	exit()

path_arcs = args[0]
path_probs = args[1]



nodes = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

with open(path_arcs, "r") as f:
	lines = f.readlines()
	
arcs = []
for line in lines:
	line = line.strip()
	if line.count(" ") != 1:
		print("Line '%s' is not a proper arc." % line)
		exit()
	arcs.append(tuple(line.split(" ")))

all_arcs = set((u,v) for u in nodes for v in nodes if u != v)
accounted_arcs = set()

for arc in arcs:
	if arc not in all_arcs:
		print("Error: '%s %s' is not a proper arc." % arc)
		exit()
	if arc in accounted_arcs:
		print("Error: Arc '%s %s' appears multiple times." % arc)
		exit()
	accounted_arcs.add(arc)

for arc in sorted(all_arcs):
	if arc not in accounted_arcs:
		print("Error: Arc '%s %s' is missing." % arc)
		exit()



with open(path_probs, "r") as f:
	lines = [l.strip() for l in f.readlines()]

probs = []
for line in lines:
	try:
		probs.append(float(line))
	except ValueError:
		print("Error: Line '%s' cannot be parsed as a number." % line)
		exit()

if any(p <= 0 for p in probs):
	print("Error: All probabilities are not positive.")
	exit()

psum = sum(probs)
diff = abs(psum - 1)
threshold = 0.00001
if diff > threshold:
	print("Error: The probabilities sum up to %s (should sum up to 1)." % psum)
	print()
	print("Normalize the probabilities so that the difference to 1 is at most %s." % threshold)
	exit()



print()
print("Your files are of correct form!")
print()
print("Make sure that the files you return are named:")
print("  firstname-lastname_arcs.txt")
print("  firstname-lastname_probs.txt")
print("  firstname-lastname_diary.pdf")
print("  firstname-lastname_summary.txt")

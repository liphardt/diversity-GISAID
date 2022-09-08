import argparse
import sys
parser = argparse.ArgumentParser()
parser.add_argument("foo")
args=parser.parse_args()

idea_states = ["AK","AR","DE","HI","ID","KS","KY","LA","ME","MS","MT","NE","NV","NH","NM","ND","OK","PR","RI","SC","SD","VT","WV","WY"]

with open(args.foo,"r") as handle:
    data = handle.readlines()

sys.stdout.write(data[0])

for i in range(1,len(data)):
    cut = data[i].strip().split(",")
    if cut[1] in idea_states:
        sys.stdout.write(data[i])


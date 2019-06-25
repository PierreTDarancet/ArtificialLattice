with open("nontrivial.dat","r") as f: 
    for line in f: 
        col = line.split() 
        print(len(col))

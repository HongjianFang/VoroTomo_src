with open('sources.in') as f:
    for line in f:
        ls = line.split()
        if float(ls[0])>2:
            print line

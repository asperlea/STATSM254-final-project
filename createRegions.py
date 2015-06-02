__author__ = 'adriana'

import glob
from collections import defaultdict
import matplotlib.pyplot as plt

PCnum = {"DNase": 3, "H3K4me1": 2, "H3K4me3": 4, "H3K9ac": 4, "H3K9me3": 2, "H3K27ac": 3, "H3K27me3": 2, "H3K36me3": 4, "RNAseq": 2}

def getRegions(valsToBins):
    # Sort bins based on contribution to PC
    important = sorted(valsToBins.keys(), reverse=True)
    # Cutoff at the ones that contribute at least 5%. Might play with this.
    x = 0
    while (x < len(important)) and (important[x] > 0.05):
        x += 1
    important = important[:x]
    print len(important), "  bins left."

    ends = []
    for val in important:
        ends.append(int(valsToBins[val][0]))
        ends.append(int(valsToBins[val][1]))
    ends.sort()
    ends = list(sorted(set(ends)))

    return ends
    '''
    regions = []
    curint = (ends[0], ends[0])
    for i in range(2, len(ends)):
        if ends[i] == curint[1] + 500:
            curint = (curint[0], ends[i])
        else:
            regions.append(curint)
            curint = (ends[i], ends[i])
    return regions
    '''

def main():
    ofile = open("importantBins.txt", "w")

    for filename in glob.glob("./Course_data/Output/PCBins/*.txt"):
        mod = filename.split("/")[4].split(".")[0]
        print mod
        ofile.write(mod + "\n")
        PCs = PCnum[mod]

        PCA = open(filename, "r")
        header = PCA.readline()
        mod = header.split(".")[0]

        PCA.readline() # read PC1
        valsToBins = defaultdict(float)
        ends = []
        line = PCA.readline()

        donePCs = 0
        fig = plt.figure()

        while (line and (donePCs < PCs)):
            if line.find("PC") != -1:
                donePCs += 1
                total = 0.0
                regions = getRegions(valsToBins)
                if len(regions) > 0:
                    b = [x for x in range(min(regions), max(regions), 5000)]
                    print len(b)
                    newhist = fig.add_subplot(PCs, 1, donePCs)
                    (n, newbins, patches) = newhist.hist(regions, bins = b)

                    #get important bins
                    ofile.write("PC" + str(donePCs) + "\n")
                    for i in range(len(newbins) - 1):
                        if n[i] > 0:
                            ofile.write(str(newbins[i]) + "," + str(n[i]) + "\n")
                    newhist.set_title("PC" + str(donePCs))
                    newhist.set_ylabel("Frequency")
                    if (donePCs == PCs):
                        newhist.set_xlabel("Position in DNA")
                print "PC", donePCs, ":"
                #print regions
                valsToBins = {}
            else:
                line = line.replace(';', ',')
                r = line.split(",")
                valsToBins[abs(float(r[-1]))] = (float(r[1]), float(r[2]))
            line = PCA.readline()

        plt.savefig(mod + ".jpg", type="jpg")
        plt.close()
    ofile.close()
main()

# look at which regions are in multiple histones
# look at p-values in those locations, if they're high enough look into them on the epigenome roadmap
__author__ = 'adriana'

import glob
from collections import defaultdict
import matplotlib.pyplot as plt

PCnum = {"DNase": 3, "H3K4me1": 5, "H3K4me3": 3, "H3K9ac": 4, "H3K9me3": 5, "H3K27ac": 4, "H3K27me3": 3, "H3K36me3": 5, "RNAseq": 6}
#PCnum = {"DNase": 1, "H3K4me1": 1, "H3K4me3": 1, "H3K9ac": 1, "H3K9me3": 1, "H3K27ac": 1, "H3K27me3": 1, "H3K36me3": 1, "RNAseq": 1}

def getEnhancers(enhancerDir):
    enhancers = {}
    for filename in glob.glob(enhancerDir):
        cellType = filename[32:-7]
        enhancers[cellType] = {}
        enhancerFile = open(filename, "r")
        for line in enhancerFile:
            start = int(line.split()[1])
            end = int(line.split()[2])
            enhancers[cellType][(start, end)] = line.split()[3]
    return enhancers

def getRegions(valsToBins):
    # Sort bins based on contribution to PC
    important = sorted(valsToBins.keys(), reverse=True)
    # Cutoff at the ones that contribute at least 5%. Might play with this.
    x = 0
    while (x < len(important)) and (important[x] > 0.05):
        x += 1
    important = important[:100] # taking top 10. Other ideas?
    print len(important), "  bins left."

    ends = []
    for val in important:
        ends.append(int(valsToBins[val][0]))
        ends.append(int(valsToBins[val][1]))
    ends.sort()
    ends = list(sorted(set(ends)))

    return ends

def main():
    ofile = open("importantBins.txt", "w")
    enhancerDir = "./ChromHMM/*.bed"
    enhancers = getEnhancers(enhancerDir)
    acrossHistones = defaultdict(list)

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

                    ofile.write("PC" + str(donePCs) + "\n")
                    for j in range(len(n)):
                        if n[j] > 5:
                            i = newbins[j]
                            print i
                            for cellType in enhancers:
                                for enhancer in enhancers[cellType]:
                                    start = enhancer[0]
                                    end = enhancer[1]
                                    if (start <= i) and (i <= end):
                                        print cellType, i, mod, donePCs, enhancer, enhancers[cellType][enhancer]


                        ofile.write(str(i) + "\n")
                        acrossHistones[i].append((mod))

                    newhist.set_title("PC" + str(donePCs))
                    newhist.set_ylabel("Frequency")
                    if (donePCs == PCs):
                        newhist.set_xlabel("Position in DNA")

                print "PC", donePCs, ":"
                #print regions
                valsToBins = {}
            else:
                r = line.strip().split(",")
                start = r[0].split(".")[1]
                end = r[0].split(".")[2]
                valsToBins[abs(float(r[-1]))] = (float(start), float(end))
            line = PCA.readline()

        plt.savefig(mod + ".jpg", type="jpg")
        plt.close()


    for region in acrossHistones:
        if len(acrossHistones[region]) > 1:
            print region, acrossHistones[region]

    ofile.close()
main()

# look at which regions are in multiple histones
# look at p-values in those locations, if they're high enough look into them on the epigenome roadmap
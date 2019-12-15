#!/usr/local/bin/python3
# Group 2

import math


# normalization formula implemented in a function
def noteNormalization(gene1IC, gene2IC, parIC):
    if parIC == 0:
        return 0
    else:
        num = float(2 * math.log10(parIC))
        den = float(math.log10(gene1IC) + math.log10(gene2IC))
        return float(num / den)


def ssScores(gene1, gene2):
    # reading in the mf child to parent file
    inFile = open("mf_CtoP_trans.txt", "r")

    # mfCToP: key = children terms, values: parent terms
    mfCToP = dict()
    mfRoot = "GO:0003674"

    # reading in the values and storing in a dictionary
    for currLine in inFile:
        if currLine[0] != 'C':
            lineList = currLine.split()
            mfCToP[lineList[0][0:len(lineList[0]) - 1]] = set()
            for value in lineList:
                if value != lineList[0]:
                    mfCToP[lineList[0][0:len(lineList[0]) - 1]].add(value)

    # closing the input file
    inFile.close()

    # ensuring that the dictionary is truly transitive
    for child in mfCToP:
        if mfRoot != child:
            while mfRoot not in mfCToP[child]:
                orig = set()
                for parent in mfCToP[child]:
                    if parent in mfCToP:
                        for ancestor in mfCToP[parent]:
                            orig.add(ancestor)
                if orig != mfCToP[child]:
                    for thing in orig:
                        mfCToP[child].add(thing)

    # reading in the bp child to parent file
    inFile = open("bp_CtoP_trans.txt", "r")

    # bpCToP: key = children terms, values: parent terms
    bpCToP = dict()
    bpRoot = "GO:0008150"

    # reading in the values and storing in dictionary
    for currLine in inFile:
        if currLine[0] != 'C':
            lineList = currLine.split()
            bpCToP[lineList[0][0:len(lineList[0]) - 1]] = list()
            for value in lineList:
                if value != lineList[0]:
                    bpCToP[lineList[0][0:len(lineList[0]) - 1]].append(value)

    # closing the input file
    inFile.close()

    # ensuring that the dictionary is truly transitive
    for child in bpCToP:
        if child != bpRoot:
            while bpRoot not in bpCToP[child]:
                origbp = set()
                for parent in bpCToP[child]:
                    if parent in bpCToP:
                        for ancestor in bpCToP[parent]:
                            origbp.add(ancestor)
                if bpCToP[child] != origbp:
                    for thing in origbp:
                        bpCToP[child].append(thing)

    # bpGeneToTermwithIEA: key - gene name, values - GO terms annotated to that gene
    bpGeneToTermwithIEA = dict()

    # reading in the BP annotations
    inFile = open("bp_GeneToTerm.txt")
    for currLine in inFile:
        if currLine != "Gene\t\t Term\n":
            lineList = currLine.split()
            bpGeneToTermwithIEA[lineList[0][0:len(lineList[0]) - 1]] = list()
            for value in lineList:
                if value != lineList[0]:
                    bpGeneToTermwithIEA[lineList[0][0:len(lineList[0]) - 1]].append(value)
    inFile.close()

    # including transitivity into the Gene to Term file
    for gene in bpGeneToTermwithIEA:
        if gene in bpCToP and gene != bpRoot:
            for a in bpCToP[gene]:
                if a not in bpGeneToTermwithIEA[gene]:
                    bpGeneToTermwithIEA[gene].append(a)

    # mfGeneToTermwithIEA: key - gene name, values - GO terms annotated to that gene
    mfGeneToTermwithIEA = dict()

    # reading in the MF annotations
    inFile = open("mf_GeneToTerm.txt")
    for currLine in inFile:
        lineList = currLine.split()
        mfGeneToTermwithIEA[lineList[0][0:len(lineList[0]) - 1]] = set()
        for value in lineList:
            if value != lineList[0]:
                mfGeneToTermwithIEA[lineList[0][0:len(lineList[0]) - 1]].add(value)
    inFile.close()

    # including transitivity into the Gene to Term file
    for gene in mfGeneToTermwithIEA:
        if gene in mfCToP and gene != mfRoot:
            for a in mfCToP[gene]:
                if a not in mfGeneToTermwithIEA[gene]:
                    mfGeneToTermwithIEA[gene].add(a)

    # bpGeneToTermwithIEA: key - gene name, values - GO terms annotated to that gene
    bpGeneToTermwithoutIEA = dict()

    # reading in the BP annotations
    inFile = open("bp_GeneToTerm_withoutIEA.txt")
    for currLine in inFile:
        if currLine != "Gene\t\t Term\n":
            lineList = currLine.split()
            bpGeneToTermwithoutIEA[lineList[0][0:len(lineList[0]) - 1]] = list()
            for value in lineList:
                if value != lineList[0]:
                    bpGeneToTermwithoutIEA[lineList[0][0:len(lineList[0]) - 1]].append(value)
    inFile.close()

    # including transitivity into the Gene to Term file
    for gene in bpGeneToTermwithoutIEA:
        if gene in bpCToP and gene != bpRoot:
            for a in bpCToP[gene]:
                if a not in bpGeneToTermwithoutIEA[gene]:
                    bpGeneToTermwithoutIEA[gene].append(a)

    # mfGeneToTermwithIEA: key - gene name, values - GO terms annotated to that gene
    mfGeneToTermwithoutIEA = dict()

    # reading in the MF annotations
    inFile = open("mf_GeneToTerm.txt")
    for currLine in inFile:
        lineList = currLine.split()
        mfGeneToTermwithoutIEA[lineList[0][0:len(lineList[0]) - 1]] = set()
        for value in lineList:
            if value != lineList[0]:
                mfGeneToTermwithoutIEA[lineList[0][0:len(lineList[0]) - 1]].add(value)
    inFile.close()

    # including transitivity into the Gene to Term file
    for gene in mfGeneToTermwithoutIEA:
        if gene in mfCToP and gene != mfRoot:
            for a in mfCToP[gene]:
                if a not in mfGeneToTermwithoutIEA[gene]:
                    mfGeneToTermwithoutIEA[gene].append(a)

    #bpGeneToTermwithIEA: key - gene name, values - GO terms annotated to that gene
    bpTermToGenewithIEA = dict()

    # reading in the BP annotations
    inFile = open("bp_TermToGene_withIEA.txt")
    for currLine in inFile:
        if currLine != "Gene\t\t Term\n":
            lineList = currLine.split()
            bpTermToGenewithIEA[lineList[0][0:len(lineList[0])]] = list()
            for value in lineList:
                if value != lineList[0]:
                    bpTermToGenewithIEA[lineList[0][0:len(lineList[0])]].append(value)
    inFile.close()

    # bpGeneToTermwithoutIEA: key - gene name, values - GO terms annotated to that gene
    bpTermToGenewithoutIEA = dict()

    # reading in the BP annotations
    inFile = open("bp_TermToGene_withoutIEA.txt")
    for currLine in inFile:
        if currLine != "Gene\t\t Term\n":
            lineList = currLine.split()
            bpTermToGenewithoutIEA[lineList[0][0:len(lineList[0])]] = list()
            for value in lineList:
                if value != lineList[0]:
                    bpTermToGenewithoutIEA[lineList[0][0:len(lineList[0])]].append(value)
    inFile.close()

    # mfGeneToTermwithIEA: key - gene name, values - GO terms annotated to that gene
    mfTermToGenewithIEA = dict()

    # reading in the MF annotations
    inFile = open("mf_TermToGene_withIEA.txt")
    for currLine in inFile:
        lineList = currLine.split()
        mfTermToGenewithIEA[lineList[0][0:len(lineList[0])]] = set()
        for value in lineList:
            if value != lineList[0]:
                mfTermToGenewithIEA[lineList[0][0:len(lineList[0])]].add(value)
    inFile.close()

    # mfGeneToTermwithoutIEA: key - gene name, values - GO terms annotated to that gene
    mfTermToGenewithoutIEA = dict()

    # reading in the MF annotations
    inFile = open("mf_TermToGene_withoutIEA.txt")
    for currLine in inFile:
        lineList = currLine.split()
        mfTermToGenewithoutIEA[lineList[0][0:len(lineList[0])]] = set()
        for value in lineList:
            if value != lineList[0]:
                mfTermToGenewithoutIEA[lineList[0][0:len(lineList[0])]].add(value)
    inFile.close()

    # creating the sets to save the GO terms for each gene
    bpGene1GOwithIEA = set()
    bpGene2GOwithIEA = set()
    mfGene1GOwithIEA = set()
    mfGene2GOwithIEA = set()
    bpGene1GOwithoutIEA = set()
    bpGene2GOwithoutIEA = set()
    mfGene1GOwithoutIEA = set()
    mfGene2GOwithoutIEA = set()

    # values to be printed to the file
    mfSSNoIEA = -1
    mfSSIEA = -1
    bpSSNoIEA = -1
    bpSSIEA = -1

    # translating from gene to GO term
    if gene1 in bpGeneToTermwithIEA and gene2 in bpGeneToTermwithIEA:
        bpGene1GOwithIEA = bpGeneToTermwithIEA[gene1]
        bpGene2GOwithIEA = bpGeneToTermwithIEA[gene2]
    if gene1 in mfGeneToTermwithIEA and gene2 in mfGeneToTermwithIEA:
        mfGene1GOwithIEA = mfGeneToTermwithIEA[gene1]
        mfGene2GOwithIEA = mfGeneToTermwithIEA[gene2]
    if gene1 in bpGeneToTermwithoutIEA and gene2 in bpGeneToTermwithoutIEA:
        bpGene1GOwithoutIEA = bpGeneToTermwithoutIEA[gene1]
        bpGene2GOwithoutIEA = bpGeneToTermwithoutIEA[gene2]
    if gene1 in mfGeneToTermwithoutIEA and gene2 in mfGeneToTermwithoutIEA:
        mfGene1GOwithoutIEA = mfGeneToTermwithoutIEA[gene1]
        mfGene2GOwithoutIEA = mfGeneToTermwithoutIEA[gene2]

    # if the gene is in the MF w/ IEA tree
    if len(mfGene1GOwithIEA) > 0 and len(mfGene2GOwithIEA) > 0:
        # sets of all the ancestor terms
        mfAGene1 = set()
        mfAGene2 = set()

        # finding all ancestors of the terms of gene 1
        for term in mfGene1GOwithIEA:
            if term in mfCToP:
                allAncestors = mfCToP[term]
                for a in allAncestors:
                    mfAGene1.add(a)

        # finding all ancestors of the terms of gene 2
        for term in mfGene2GOwithIEA:
            if term in mfCToP:
                allAncestors = mfCToP[term]
                for a in allAncestors:
                    mfAGene2.add(a)

        # find the union of the two parents
        parU = mfAGene1.union(mfAGene2)

        # find the intersection of the two parents
        parI = mfAGene1.intersection(mfAGene2)

        # finding the information contents of gene 1
        gene1ICwIEA = float(len(mfGeneToTermwithIEA[gene1]) / len(mfGeneToTermwithIEA.values()))

        # finding the information contents of gene 2
        gene2ICwIEA = float(len(mfGeneToTermwithIEA[gene2]) / len(mfGeneToTermwithIEA.values()))

        # find the sum of the information contents for the intersection of parents
        numIEA = 0
        for parent in parI:
            if parent in mfTermToGenewithIEA:
                parICIEA = float(len(mfTermToGenewithIEA[parent]) / len(mfTermToGenewithIEA.values()))
                currIEA = noteNormalization(gene1ICwIEA, gene2ICwIEA, parICIEA)
                numIEA += currIEA

        # find the sum of the information contents for the union of parents
        denIEA = 0
        for parent in parU:
            if parent in mfTermToGenewithoutIEA:
                parICIEA = float(len(mfTermToGenewithIEA[parent]) / len(mfTermToGenewithIEA.values()))
                currIEA = noteNormalization(gene1ICwIEA, gene2ICwIEA, parICIEA)
                denIEA += currIEA

        # check if the numerator is zero before dividing
        if numIEA == 0:
            mfSSIEA = 0
        else: # divide the intersection sum by the union sum
            mfSSIEA = numIEA / denIEA

    # if the gene is in the BP w/ IEA tree
    if len(bpGene1GOwithIEA) > 0 and len(bpGene2GOwithIEA) > 0:
        # sets of all the ancestor terms
        bpAGene1 = set()
        bpAGene2 = set()

        # finding all ancestors of the terms of gene 1
        for term in bpGene1GOwithIEA:
            if term in bpCToP:
                allAncestors = bpCToP[term]
                for a in allAncestors:
                    bpAGene1.add(a)

        # finding all ancestors of the terms of gene 2
        for term in bpGene2GOwithIEA:
            if term in bpCToP:
                allAncestors = bpCToP[term]
                for a in allAncestors:
                    bpAGene2.add(a)

        # find the union of the two parents
        parU = bpAGene1.union(bpAGene2)

        # find the intersection of the two parents
        parI = bpAGene1.intersection(bpAGene2)

        # finding the information contents of gene 1
        gene1ICwIEA = float(len(bpGeneToTermwithIEA[gene1]) / len(bpGeneToTermwithIEA.values()))

        # finding the information contents of gene 2
        gene2ICwIEA = float(len(bpGeneToTermwithIEA[gene2]) / len(bpGeneToTermwithIEA.values()))

        # find the sum of the information contents for the intersection of parents
        numIEA = 0
        for parent in parI:
            if parent in bpTermToGenewithIEA:
                parICIEA = float(len(bpTermToGenewithIEA[parent]) / len(bpTermToGenewithIEA.values()))
                currIEA = noteNormalization(gene1ICwIEA, gene2ICwIEA, parICIEA)
                numIEA += currIEA

        # find the sum of the information contents for the union of parents
        denIEA = 0
        for parent in parU:
            if parent in bpTermToGenewithoutIEA:
                parICIEA = float(len(bpTermToGenewithIEA[parent]) / len(bpTermToGenewithIEA.values()))
                currIEA = noteNormalization(gene1ICwIEA, gene2ICwIEA, parICIEA)
                denIEA += currIEA

        # check if the numerator is 0 before dividing
        if numIEA == 0:
            bpSSIEA = 0
        else: # divide the intersection sum by the union sum
            bpSSIEA = numIEA / denIEA

    # if the gene is in the MF w/o IEA tree
    if len(mfGene1GOwithoutIEA) > 0 and len(mfGene2GOwithoutIEA) > 0:
        # sets of all the ancestor terms
        mfAGene1 = set()
        mfAGene2 = set()

        # finding all ancestors of the terms of gene 1
        for term in mfGene1GOwithoutIEA:
            if term in mfCToP:
                allAncestors = mfCToP[term]
                for a in allAncestors:
                    mfAGene1.add(a)

        # finding all ancestors of the terms of gene 2
        for term in mfGene2GOwithoutIEA:
            if term in mfCToP:
                allAncestors = mfCToP[term]
                for a in allAncestors:
                    mfAGene2.add(a)

        # find the union of the two parents
        parU = mfAGene1.union(mfAGene2)

        # find the intersection of the two parents
        parI = mfAGene1.intersection(mfAGene2)

        # finding the information contents of gene 1
        gene1ICwoIEA = float(len(mfGeneToTermwithoutIEA[gene1]) / len(mfGeneToTermwithoutIEA.values()))

        # finding the information contents of gene 2
        gene2ICwoIEA = float(len(mfGeneToTermwithoutIEA[gene2]) / len(mfGeneToTermwithoutIEA.values()))

        # find the sum of the information contents for the intersection of parents
        numNoIEA = 0
        for parent in parI:
            if parent in mfTermToGenewithoutIEA:
                parICwoIEA = float(len(mfTermToGenewithoutIEA[parent]) / len(mfTermToGenewithoutIEA.values()))
                currNoIEA = noteNormalization(gene1ICwoIEA, gene2ICwoIEA, parICwoIEA)
                numNoIEA += currNoIEA

        # find the sum of the information contents for the union of parents
        denNoIEA = 0
        for parent in parU:
            if parent in mfTermToGenewithoutIEA:
                parICwoIEA = float(len(mfTermToGenewithoutIEA[parent]) / len(mfTermToGenewithoutIEA.values()))
                currNoIEA = noteNormalization(gene1ICwoIEA, gene2ICwoIEA, parICwoIEA)
                denNoIEA += currNoIEA

        # check to see if the numerator is zero before dividing
        if numNoIEA == 0:
            mfSSNoIEA = 0
        else: # divide the intersection sum by the union sum
            mfSSNoIEA = numNoIEA / denNoIEA

    # if the gene is in the BP w/o IEA tree
    if len(bpGene1GOwithoutIEA) > 0 and len(bpGene2GOwithoutIEA) > 0:
        # sets of all the ancestor terms
        bpAGene1 = set()
        bpAGene2 = set()

        # finding all ancestors of the terms of gene 1
        for term in bpGene1GOwithoutIEA:
            if term in bpCToP:
                allAncestors = bpCToP[term]
                for a in allAncestors:
                    bpAGene1.add(a)

        # finding all ancestors of the terms of gene 2
        for term in bpGene2GOwithoutIEA:
            if term in bpCToP:
                allAncestors = bpCToP[term]
                for a in allAncestors:
                    bpAGene2.add(a)

        # find the union of the two parents
        parU = bpAGene1.union(bpAGene2)

        # find the intersection of the two parents
        parI = bpAGene1.intersection(bpAGene2)

        # finding the information contents of gene 1
        gene1ICwoIEA = float(len(bpGeneToTermwithoutIEA[gene1]) / len(bpGeneToTermwithoutIEA.values()))

        # finding the information contents of gene 2
        gene2ICwoIEA = float(len(bpGeneToTermwithoutIEA[gene2]) / len(bpGeneToTermwithoutIEA.values()))

        # find the sum of the information contents for the intersection of parents
        numNoIEA = 0
        for parent in parI:
            if parent in bpTermToGenewithoutIEA:
                parICwoIEA = float(len(bpTermToGenewithoutIEA[parent]) / len(bpTermToGenewithoutIEA.values()))
                currNoIEA = noteNormalization(gene1ICwoIEA, gene2ICwoIEA, parICwoIEA)
                numNoIEA += currNoIEA

        # find the sum of the information contents for the union of parents
        denNoIEA = 0
        for parent in parU:
            if parent in bpTermToGenewithoutIEA:
                parICwoIEA = float(len(bpTermToGenewithoutIEA[parent]) / len(bpTermToGenewithoutIEA.values()))
                currNoIEA = noteNormalization(gene1ICwoIEA, gene2ICwoIEA, parICwoIEA)
                denNoIEA += currNoIEA

        # divide the intersection sum by the union sum
        bpSSNoIEA = numNoIEA / denNoIEA

    return max(bpSSNoIEA, bpSSIEA, mfSSNoIEA, mfSSIEA)


#!/usr/local/bin/python3
# Group 2

import os
import json
import cgi
import math

print('Content-type: application/json;charset=utf-8')
print()

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

# bpGeneToTermwithIEA: key - gene name, values - GO terms annotated to that gene
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

# normalization formula implemented in a function
def noteNormalization(gene1IC, gene2IC, parIC):
    if parIC == 0:
        return 0
    else:
        num = float(2 * math.log10(parIC))
        den = float(math.log10(gene1IC) + math.log10(gene2IC))
        return float(num / den)


def ssScores(gene1, gene2):
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


def findShortestPath(graph, start, end):
    shortest_paths = {start: None}
    current_node = start
    visited = set()

    while current_node != end:
        visited.add(current_node)
        destinations = graph[current_node]

        for next_node in destinations:
            if next_node not in shortest_paths:
                shortest_paths[next_node] = current_node
        next_destinations = {node: shortest_paths[node] for node in shortest_paths if node not in visited}
        if not next_destinations:
            return "no connection"
        # next node is the destination with the lowest weight
        current_node = min(next_destinations, key=lambda k: next_destinations[k])

    # Work back through destinations in shortest path
    path = []
    while current_node is not None:
        path.append(current_node)
        next_node = shortest_paths[current_node]
        current_node = next_node
    # Reverse path
    path = path[::-1]
    return path


# initializing the synonym table
synonyms = dict()

# opening the the file that saves the synonym table
with open("synTable.txt") as file:
    # reading in the information from the file
    file = file.read()

    # splitting the data into lines
    data = file.split("\n")

    # parsing through the lines
    for line in data:
        # splitting the lines into individual terms
        afterSplit = line.split(" ")

        # checking to see if the current value is in the dictionary
        if afterSplit[0] not in synonyms.keys():
            synonyms[afterSplit[0]] = set()

        # parsing through the rest of the terms on the line
        for term in afterSplit:
            # checking to see that the term is not the key value
            if term != afterSplit[0] and term != "":
                # adding the term to the set
                synonyms[afterSplit[0]].add(term.upper())


# getting the two user entered genes
form = cgi.FieldStorage()
enteredGene1 = form["startGene"].value
enteredGene2 = form["endGene"].value

# initializing the synonym values
geneSyn1 = ""
geneSyn2 = ""

# finding synonym for the gene 1
if enteredGene1 not in synonyms.keys():
    # parsing the terms in the synonyms table
    for term in synonyms.keys():
        # parsing the synonyms to gene 1
        for syn in synonyms[term]:
            # check to see if the current synonym is the entered term
            if syn == enteredGene1:
                # assigning the top term if the entered term is a synonym
                geneSyn1 = term
                break
        if geneSyn1 == term:
            break
elif enteredGene1 in synonyms.keys():
    geneSyn1 = enteredGene1

# finding synonym for the gene 2
if enteredGene2 not in synonyms.keys():
    # parsing the terms in the synonyms table
    for term in synonyms.keys():
        # parsing the synonyms to gene 2
        for syn in synonyms[term]:
            # check to see if the current synonym is the entered term
            if syn == enteredGene2:
                # assigning the top term if the entered term is a synonym
                geneSyn2 = term
                break
        if geneSyn2 == term:
            break
elif enteredGene2 in synonyms.keys():
    geneSyn2 = enteredGene2

if geneSyn1 == "" and enteredGene1 not in synonyms.keys():
    print("error: gene 1 not in database")
elif geneSyn2 == "" and enteredGene2 not in synonyms.keys():
    print("error: gene 2 not in database")
elif geneSyn1 != geneSyn2:
    # creating the node and edges sets
    nodes = dict()
    edges = dict()

    # reassigning variables
    gene1 = geneSyn1
    gene2 = geneSyn2

    # adding gene 1 and gene 2 into the nodes set
    info = dict()
    info['id'] = enteredGene1
    info['process'] = "this is a test"

    info2 = dict()
    info2['id'] = enteredGene2
    info2['process'] = "this is a test2"

    syn1ToPrint = list()
    syn2ToPrint = list()

    for i in synonyms[gene1]:
        if i != enteredGene1:
            syn1ToPrint.append(i)
        else:
            syn1ToPrint.append(gene1)

    for i in synonyms[gene2]:
        if i != enteredGene2:
            syn2ToPrint.append(i)
        else:
            syn2ToPrint.append(gene2)

    info['synonyms'] = syn1ToPrint
    info2['synonyms'] = syn2ToPrint

    nodes[enteredGene1] = info
    nodes[enteredGene2] = info2

    # opening the complexes file
    with open("complexes_human_cln_flt.txt") as data:
        # reading in the data from the complexome data file
        data = data.read()

        # splitting the life into lines
        data = data.split("\n")

        # reading through each line of the file
        for line in data:
            # splitting the line into genes
            contents = line.split(" ")

            # initializing bool variables
            isGene1 = 0
            isGene2 = 0

            # reading each gene in the line
            for gene in contents:
                # check if curr gene is gene 1
                if gene == gene1:
                    isGene1 = 1

                # check if curr gene is gene 2
                if gene == gene2:
                    isGene2 = 1

                # check if both genes have been found
                if isGene1 == 1 and isGene2 == 1:
                    # iterating through each gene again
                    for neighbor in contents:
                        # checking to add nodes & edges only for non-user-entered genes
                        if neighbor != gene1 and neighbor != gene2 and neighbor != '':
                            # adding the curr neighbor to the nodes dictionary

                            if neighbor not in synonyms.keys():
                                # parsing the terms in the synonyms table
                                for term in synonyms.keys():
                                    # parsing the synonyms to gene 2
                                    for syn in synonyms[term]:
                                        # check to see if the current synonym is the entered term
                                        if syn == neighbor:
                                            # assigning the top term if the entered term is a synonym
                                            geneSyn3 = term
                                            break
                                    if geneSyn3 == term:
                                        break
                            elif neighbor in synonyms.keys():
                                geneSyn3 = neighbor

                            info3 = dict()
                            info3['id'] = neighbor
                            info3['process'] = "this is a test for 3"

                            syn3ToPrint = list()

                            for i in synonyms[geneSyn3]:
                                if i != neighbor:
                                    syn3ToPrint.append(i)
                                else:
                                    syn3ToPrint.append(geneSyn3)

                            info3['synonyms'] = syn3ToPrint

                            nodes[neighbor] = info3

                            # creating the edge list
                            currEdge1 = list()

                            # storing the starting node
                            currEdge1.append(enteredGene1)

                            # storing the ending nose
                            currEdge1.append(neighbor)

                            # creating edge ID
                            edgeID = enteredGene1 + "-" + neighbor

                            ssEdge1 = ssScores(enteredGene1, neighbor)
                            currEdge1.append(ssEdge1)

                            # storing the edge id as key in dictionary
                            edges[edgeID] = currEdge1

                            # creating edge list for edge to gene2
                            currEdge2 = list()

                            # storing the starting node
                            currEdge2.append(enteredGene2)

                            # storing the ending nose
                            currEdge2.append(neighbor)

                            # creating edge ID
                            edgeID = enteredGene2 + "-" + neighbor

                            ssEdge2 = ssScores(enteredGene2, neighbor)
                            currEdge2.append(ssEdge2)

                            # storing the edge id as key in dictionary
                            edges[edgeID] = currEdge2

                    # resetting the bool values
                    isGene1 = 0
                    isGene2 = 0

                # checking to see if the current gene value is the last of the line
                if gene == contents[len(contents) - 1]:
                    # resetting the bool values
                    isGene1 = 0
                    isGene2 = 0

    if len(edges) == 0:
        shortestPathGraph = dict()
        with open("complexes_human_cln_flt.txt") as data:
            # reading in the data from the complexome data file
            data = data.read()

            # splitting the life into lines
            data = data.split("\n")

            # reading through each line of the file
            for line in data:
                contents = line.split(" ")

                # going through each term in line and adding it to the term's graph
                for term in contents:
                    # adding in the term to dict if not there
                    if term not in shortestPathGraph:
                        # initializing the set
                        shortestPathGraph[term] = set()

                        # adding in the rest of the line into the set to show connection
                        for inTerm in contents:
                            if inTerm != "":
                                shortestPathGraph[term].add(inTerm)
                            else:
                                # adding in the rest of the line into the set to show connection
                                for inTerm in contents:
                                    if inTerm != "":
                                        shortestPathGraph[term].add(inTerm)
        for item in shortestPathGraph:
            shortestPathGraph[item] = list(shortestPathGraph[item])

        path = findShortestPath(shortestPathGraph, gene1, gene2)

        index = 0
        while index < len(path):
            # adding gene 1 and gene 2 into the nodes set
            info = dict()
            info['id'] = path[index]

            syn1ToPrint = list()

            for i in synonyms[path[index]]:
                syn1ToPrint.append(i)

            info['synonyms'] = syn1ToPrint

            nodes[path[index]] = info

            index += 1

        index = 1
        while index < len(path):
            # creating the edge list
            currEdge = list()

            # storing the starting node
            currEdge.append(path[index - 1])

            # storing the ending nose
            currEdge.append(path[index])

            # creating edge ID
            edgeID = path[index - 1] + "-" + path[index]

            ssEdge = ssScores(path[index - 1], path[index])
            currEdge.append(ssEdge)

            edges[edgeID] = currEdge
            index += 1

# creating the list for the json output
toJson = dict()

# creating keys
key1 = "edges"
key2 = "nodes"

# adding the nodes into the json output
toJson[key2] = nodes

# adding the edges to the json output
toJson[key1] = edges

# printing out the edges and the nodes in json
print(json.dumps(toJson, indent=4))





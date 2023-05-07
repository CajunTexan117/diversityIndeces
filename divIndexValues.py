#! /usr/bin/env python3

# diversityIndex.py
# Authors: Andres Herrera, Philip Porter, Matthew Stout

# This script reads in a CSV and runs that data through
# a series of functions to return diversity index values 
# and three graphs: A rank-abundance plot, Comparison of
# Species Richness and Chao value, and a comparison of
# Diversity Indeces for each location

# Modules
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import csv
import argparse
import matplotlib.pyplot as plt

# Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# shannon function takes absolute abundance value and returns shannon index
def shannon(absolute_abundances):
    relative_abundances = []
    log_relative_abundances = []
    summary = []
    for i in absolute_abundances:
        Pi = i / sum(absolute_abundances)
        relative_abundances.append(Pi)
    for y in relative_abundances:
        if y <= 0:
            Pii = 0
        else:
            Pii = np.log(y)
        log_relative_abundances.append(Pii)
    for z in range(len(log_relative_abundances)):
        summary.append((relative_abundances[z] * log_relative_abundances[z]) * -1)
        H = sum(summary)    
    return H

# effective_H returns expected Shannon value
def effective_H(absolute_abundances):
    H = shannon(absolute_abundances)
    exp_H = np.exp(H)
    return exp_H

# simpson returns Simpson index value
def simpson(absolute_abundances):
    N = sum(absolute_abundances)
    new_list = []
    for i in absolute_abundances:
            val = (i * (i - 1))
            new_list.append(val)
    S = sum(new_list)
    D = 1 - (S / (N * (N - 1)))
    return D

# richness returns species richness value
def richness(absolute_abundances):
    r = [x for x in absolute_abundances if x != 0]
    S = len(r)
    return S

# species richness for an individual row
def rowRichness(row):
    return sum([1 for i in row if i > 0])

# pielou returns pielou evenness value
def pielou(absolute_abundances):
    H = shannon(absolute_abundances)
    S = richness(absolute_abundances)
    P = H / np.log(S)
    return P

# margalef returns margalef diversity index
def margalef(absolute_abundances):
    S = len(absolute_abundances)
    N = sum(absolute_abundances)
    M = ( S - 1) / np.log(N)
    return M

 # menhenik returns menhenik value
def menhenik(absolute_abundances):
    S = len(absolute_abundances)
    N = sum(absolute_abundances)
    ME = S / np.sqrt(N)
    return ME

# chao returns Chao diversity index
def chao(absolute_abundances):
    S = len(absolute_abundances)
    n1 = absolute_abundances.count(1)
    n2 = absolute_abundances.count(2)
    C = S + ((n1 * (n1 - 1)) / (2 * (n2 + 1)))
    return C

# chao returns Chao diversity index
def rowChao(absolute_abundances, richness):
    S = richness
    n1 = absolute_abundances.count("1")
    n2 = absolute_abundances.count("2")
    #print(S, n1, n2)
    C = S + ((n1 * (n1 - 1)) / (2 * (n2 + 1)))
    return C

# Main Block
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#file input
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="the name of the input file")
parser.add_argument("diversity_index", help="Diverity index to graph in comparing sites")
args = parser.parse_args()
inputFile = args.input_file
Index = args.diversity_index

#open csv file and import as dictionary of items and values
with open(inputFile, "r", newline="") as file:
    communitiesData = csv.DictReader(file, )
    communitiesDataList = list(communitiesData)
    
    #list of species populations
    speciesPops = []
    columnNames = communitiesData.fieldnames
    
    #loop to sum species populations    
    for col in columnNames:
        # sum the numeric values in the row
        col_sum = 0
        for row in communitiesDataList:
            col_sum += int(row[col])
        speciesPops.append(col_sum)
    
    #list of community populations
    communityPops = []
    
    # loop to sum community populations
    for i, row in enumerate(communitiesData):

        # sum the numeric values in the row
        row_sum = 0
        for value in row.values():
            try:
                row_sum += int(value)
            except ValueError:
                pass  # ignore non-numeric values
        
        # add the row sum to the sums list
        communityPops.append(row_sum)
        
    #Make a list of values with the diversity index for each row in the .csv file
    rowIndex = []
    for dic in communitiesDataList:
        x = []
        for key in dic:
            num = int(dic[key])
            x.append(num)
        rowIndex.append(x)
        
    #print(x_values)
    
    #Iterate through each list in x_values and create an index based on 
    #the one inputed by the user on the command line
    div_x = []
    for li in rowIndex:
        if sum(li) == 0:
            div_x.append(0)
            continue
        if Index == 'H':
            div = shannon(li)
            div_x.append(div)
        if Index == 'exp_H':            
            div = effective_H(li)
            div_x.append(div)
        if Index == 'D':
            div = simpson(li)
            div_x.append(div)
        if Index == 's':
            div = richness(li)
            div_x.append(div)
        if Index == 'P':
            div = pielou(li)
            div_x.append(div)
        if Index == 'M':
            div = margalef(li)
            div_x.append(div)
        if Index == 'Me':
            div = menhenik(li)
            div_x.append(div)
        if Index == 'chao':
            div = chao(li)
            div_x.append(div)
             
#print off a list of all index values
H = shannon(speciesPops)
print("Shannon:", H)

exp_H = effective_H(speciesPops)
print("Exp_Shannon:", exp_H)

D = simpson(speciesPops)
print("Simpson:", D)

s = richness(speciesPops)
print("Richness:", s)  

P = pielou(speciesPops)
print("Evenness:", P)

M = margalef(speciesPops)
print("Margalef:", M)

Me = menhenik(speciesPops)
print("Menhenik:", Me)

chaoprint = chao(speciesPops)
print("Chao:", chaoprint)

# Plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot of index specified in command line argument
# Here is the data plotted in an axes
x = np.arange(1, len(communitiesDataList)+1)
y = div_x

# Creates the plot:
fig, ax = plt.subplots()
ax.errorbar(x, y, color='forestgreen', fmt='o', linewidth=2, capsize=6)
#Adds labels
ax.set_xlabel("Sampling Sites")
ax.set_ylabel(Index)
ax.set_title("Comparison of Diversity Index for each Sampling Site")
# This sets it so that every individual sampleing site is shown on the x axis
plt.xticks(x)

plt.show()

#
# calculate species richness and chao values for each row
richness_vals = []
chao_vals = []
for row in communitiesDataList:
    # count number of species with non-zero abundances
    richness = len([int(x) for x in row.values() if int(x) > 0])
    richness_vals.append(richness)
    # calculate chao value using chao function
    feedVal = list(row.values())
    #print(feedVal)
    chao_value = rowChao(feedVal, richness)
    #print(chao_value)
    chao_vals.append(chao_value)
#print(chao_vals)
# create bar graph
x = np.arange(1, len(communitiesDataList)+1)
width = 0.35
fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, richness_vals, width, label='Species Richness')
rects2 = ax.bar(x + width/2, chao_vals, width, label='Chao Value')

# add labels and title
ax.set_ylabel('Number of Species')
ax.set_xlabel('Location')
ax.set_title('Comparison of Species Richness and Chao Value for each Location')
ax.set_xticks(x)
ax.legend()

# show the graph
plt.show()

#
# Rank abundance plot

#Create a list with the libraries
data = dict(zip(columnNames,speciesPops))

#Sort the data from the hihest to lowest abundance
sortedDic = sorted(data.items(), key=lambda x:x[1], reverse=True)
converted_dict = dict(sortedDic)

#Call the dic keys and values
names = list(converted_dict.keys())
values = list(converted_dict.values())

#Plot the rank-abun
fig, ax = plt.subplots()
plt.plot(names, values, linestyle="-", marker="o", color = "red")
locs, labels = plt.xticks()
plt.setp(labels, rotation=90)
plt.yscale("log")
ax.set_ylabel('Species abundances (Log)')
ax.set_xlabel('Species')
ax.set_title('Rank-Abundance plot')
plt.show()




















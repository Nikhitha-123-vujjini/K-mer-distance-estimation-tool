#!/usr/bin/env python3


import argparse
import csv
import math
import matplotlib.pyplot as plt
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("fasta",nargs = "+", help = "Input fasta file")
parser.add_argument("-k","--kmer_size", type=int,help = "k_mer size" )
parser.add_argument("-m","--metric", default = "euclidean", help = "Distance metric")
parser.add_argument("-o","--output_dir",default = "results", help = "output directory")
args = parser.parse_args()

def load_fasta(fasta_file):
    seqs = {}

    for file in fasta_file:
        basename = file.split(".")[0]
        #parse all sequences in FASTA file
        for record in SeqIO.parse(file,"fasta"):
            header = record.id.split("|")[-1]
            short = header.split("_")[0]
            #readable label for this sequence
            label = basename + "_" + short
            #convert the  sequence object to string
            seqs[label] = str(record.seq)

    return seqs

#extracts all the overlapping kmers of user specified kmer_size from the DNA sequence 
def extract_kmers(sequence, kmer_size):
    kmers_list = []
    end = len(sequence) - kmer_size + 1

    for idx in range(end):
        k_mer = sequence[idx:idx+kmer_size]
        kmers_list.append(k_mer)

    return kmers_list

#counts the number of occurences of each kmer in a sequence  
def count_kmers(kmer_list):
    counts = {}

    for k_mer in kmer_list:
        if k_mer in counts:
            counts[k_mer] += 1
        else:
            counts[k_mer] = 1

    return counts

#creates a master list of all k-mers across all sequences.
def unique_kmers(all_counts):
    unique_set = set()

    for count_dict in all_counts.values():
        for k_mer in count_dict.keys():
            unique_set.add(k_mer)

    unique_list = sorted(unique_set)
    return unique_list


#converts a sequences k-mer counts into a normalized numeric vector.
def count_to_vectors(counts, unique_list):

    total = sum(counts.values())

    vectors = []
    for k_mer in unique_list:
        # if the kmer appears in  the sequence, the  computes the normalized frequency
        if k_mer in counts:
            frequency = counts[k_mer]/total
        else:
            frequency = 0

        vectors.append(frequency)

    return vectors


#distance metrics for the calaculation 

def euclidean(v1,v2):
    #Euclidean(v1,v2)= √i=1∑n(v1i−v2i) ** 2
    diff = 0
    for i in range(len(v1)):
        diff += (v1[i] - v2[i]) ** 2
    return math.sqrt(diff)

def manhattan(v1,v2):
    #Manhattan(v1,v2)=i=1∑n∣v1i−v2i|
    total = 0
    for i in range(len(v1)):
        a = v1[i]
        b = v2[i]
        if a > b:
            total += (a-b)
        else:
            total += (b-a)
    return total

def cosine(v1,v2):
    #Cosine_similarity(v1,v2) = √∑i=1n v1i⋅v2i / ∑i=1n (v1i) ** 2 . √∑i=1n (v2i) ** 2
    #Cosine(v1,v2) = 1−cosine_similarity(v1,v2)
    numerator =  0
    mag1 = 0
    mag2 = 0

    for i in range(len(v1)):
        a = v1[i]
        b = v2[i]

        numerator += a*b
        mag1 += a*a
        mag2 += b*b


    mag1 = math.sqrt(mag1) 
    mag2 = math.sqrt(mag2)

    cosine_similarity = numerator/(mag1 * mag2)
    cosine_distance = 1 - cosine_similarity 
    return cosine_distance

#calculates pairwise distances between all sequences
def pairwise_comparision(all_vectors,metric):
    labels = list(all_vectors.keys())
    n = len(labels)
    matrix = []

    for i in range(n):
        row = []

        for j in range(n):
            #get the two vectors we are measuring distance between
            v1 = all_vectors[labels[i]] 
            v2 = all_vectors[labels[j]]
        
            #applies chosen distance metric
            if metric == "euclidean":
                dist = euclidean(v1,v2)
            elif metric == "manhattan":
                dist = manhattan(v1,v2)
            else:
                dist = cosine(v1,v2)

            row.append(dist)

        matrix.append(row)

    return labels, matrix

def heatmap_plot(matrix, labels, outfile):
    #visualization  of distance matrix(csv) as a heatmap
    fig, ax = plt.subplots(figsize = (10,8))

    heatmap = ax.imshow(matrix, cmap = "coolwarm")
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))

    ax.set_xticklabels(labels, rotation = 25)
    ax.set_yticklabels(labels)
    
    fig.colorbar(heatmap)
    ax.set_title("k-mer distance heatmap")
    
    fig.savefig(outfile)
    plt.close()


#main pipeline for execution 
def main():

    #create the output directory 
    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        print("output directory already exists")

    sequences = load_fasta(args.fasta)
    
    #extract and count kmers for each sequence
    all_kmer_counts = {}
    for header, sequence in sequences.items():
        kmer_list = extract_kmers(sequence,args.kmer_size)
        kmer_counts = count_kmers(kmer_list)
        all_kmer_counts[header] = kmer_counts

    #build the unique kmers
    unique_kmers_list = unique_kmers(all_kmer_counts)

    #convert each sequence to  normalized vector
    all_vectors = {}
    for header, kmer_counts in all_kmer_counts.items():
        vec = count_to_vectors(kmer_counts,unique_kmers_list)
        all_vectors[header] = vec

    #computes pairwise distances 
    labels, matrix = pairwise_comparision(all_vectors, args.metric)

    #save csv
    csv_path = "{}/distance_matrix.csv".format(args.output_dir)
    with open (csv_path, "w", newline = "") as filehandle:
        writer = csv.writer(filehandle, delimiter = ",", quotechar = "|")
        writer.writerow([""]+labels)
        for i in range(len(labels)):
            writer.writerow([labels[i]] + matrix[i])
    #save heatmap
    heatmap_path = "{}/heatmap.png".format(args.output_dir)
    heatmap_plot(matrix, labels, heatmap_path)

main()

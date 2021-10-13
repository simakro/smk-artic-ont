# Reference Guided small viral GEnome de-novo Assembly from 3rdGen AMPlicon data '
# REGGEA-AMP #

import sys
import os
from collections import Counter
import subprocess
import multiprocessing as mp
import numpy as np
import pickle
import gzip
    

class Alignment:
    def __init__(self, qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, matches, alnlen, qual):
        self.qname = qname
        self.qlen = int(qlen)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.strand = strand
        self.tname = tname
        self.tlen = int(tlen)
        self.tstart = int(tstart)
        self.tend = int(tend)
        self.matches = int(matches)
        self.alnlen = int(alnlen)
        self.qual = int(qual)
        self.read = str()


class Amplicon:
    def __init__(self, amp_name, type, amp_len, amp_start, amp_end, start_primer_boundary, end_primer_boundary, clipped_amp_len, ovl_with_prev_amp):
        self.amp_name = int(amp_name)
        self.type = type
        self.amp_len = int(amp_len)
        self.amp_start = int(amp_start)
        self.amp_end = int(amp_end)
        self.start_primer_boundary = int(start_primer_boundary)
        self.end_primer_boundary = int(end_primer_boundary)
        self.clipped_amp_len = int(clipped_amp_len)
        self.ovl_with_prev_amp = int(ovl_with_prev_amp)
        self.stacked_reads = dict() #list()


def load_reads(read_fasta):
    read_dict = dict()
    with open(read_fasta, "r") as reads:
        for line in reads:
            if line.startswith(">"):
                name_long = line.strip().split(">")[1]
                name = name_long.split(" ")[0]
                read_dict[name] = next(reads).strip()
    return read_dict


def load_mm2_paf(paf, reads):
    alignments = list()
    with open(paf,"r") as paf:
        for line in paf:
            al = Alignment(*line.strip().split("\t")[:12])
            al.read = reads[al.qname]
            alignments.append(al)
    return alignments

def load_amplicons(amp_csv):
    amplicons = list()
    with open(amp_csv, "r") as amps:
        ct = 0
        for line in amps:
            ct += 1
            if ct == 1:
               pass
            else:
                amp = Amplicon(*line.strip().split(";"))
                amplicons.append(amp)
    return amplicons

def merge_alt__amps(amps):
    prim_amps = [amp for amp in amps if amp.type == "primary"]
    alt_amps = [amp for amp in amps if amp.type == "alternative"]
    for aamp in alt_amps:
        for pamp in prim_amps:
            if aamp.amp_name == pamp.amp_name:
                pamp.amp_start = min([pamp.amp_start, aamp.amp_start])
                pamp.amp_end = max([pamp.amp_end, aamp.amp_end])
                pamp.start_primer_boundary = max([pamp.start_primer_boundary, aamp.start_primer_boundary])
                pamp.end_primer_boundary = min([pamp.end_primer_boundary, aamp.end_primer_boundary])
                pamp.amp_len = pamp.amp_end - pamp.amp_start + 1
                pamp.clipped_amp_len = pamp.end_primer_boundary - pamp.start_primer_boundary + 1
                pamp.ovl_with_prev_amp = [pamp.ovl_with_prev_amp, aamp.ovl_with_prev_amp]
    #print(prim_amps)
    return prim_amps
                

def filter_paf(als, max_amp_len):
    #in order for downstream steps to work correctly it must be ensured, that there is only one single alignement in the paf file
    #So there must still written a filter to include only the best alignemnt for each read; which is not the case yet
    print("alignments before filtering", len(als))
    len_filt = [al for al in als if al.qlen < max_amp_len*1.1]
    #print("alignments after length filtering", len(len_filt))
    q_filt = [al for al in len_filt if al.qual > 59]
    #print("alignments after quality filtering", len(q_filt))
    match_filt = [al for al in q_filt if al.matches/al.alnlen > 0.80]
    print("alignments after matches filtering", len(match_filt))
    return match_filt


def revcomp(seq):
    reverse_seq = "".join(seq[::-1]).upper()
    cbd = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "R": "Y", "Y": "R"}
    revcomp = "".join([cbd[base] for base in reverse_seq])
    return revcomp

def sort_alns_to_amps(merged_amps, filt_als):
    binned = list()
    too_long = list()
    for al in filt_als:
        for amp in merged_amps:
            if al.tend <= amp.amp_end: #amp.end_primer_boundary amp.amp_end
                if al.tstart >= amp.amp_start: #amp.start_primer_boundary amp.amp_start
                    if al.qlen < amp.clipped_amp_len * 1.1:
                        if al.strand == "+":
                            amp.stacked_reads[al.qname] = (al.read)
                        elif al.strand == "-": #in order for this to work without screwups it must be ensured, that there is only one single alignement for each read in the paf file/ALignment list
                            amp.stacked_reads[al.qname] = (revcomp(al.read))
                        else:
                            print("Unexpected strand indicator")
                    else:
                        too_long.append(al)
                    binned.append(al.qname)
    unbinnable = [al for al in filt_als if al.qname not in binned]
    print("binned:", len(binned))
    print("too long for correctly clipped amp:", len(too_long))
    print("unbinnable:", len(unbinnable))
    print("sum unbinned+binned", len(binned) + len(unbinnable))
    print("input als:", len(filt_als))

    with open("pickled_amp_binned_als.lst", "wb") as paba:
        pickle.dump(merged_amps, paba)

    return merged_amps


def write_out_bins(amp_binned_als):
    bins = list()
    folder_path = "./results/binned_reads/"
    os.makedirs(folder_path, exist_ok=True)
    for amp in amp_binned_als:
        if len(amp.stacked_reads) > 0:
            bin_path = folder_path + "Amp" + str(amp.amp_name) + "_binned_reads.fasta"
            bins.append(bin_path)
            with open(bin_path, "w") as out:
                for read in amp.stacked_reads:
                    if len(amp.stacked_reads) > 0:
                        out.write(">" + read + "\n")
                        out.write(amp.stacked_reads[read] + "\n")
        else:
            print("Amplicon 64 dropout")

    return bins


def convert_multiline_fasta(ml_in, sl_out):
	converted = open(sl_out, "a")
	row=0
	with open(ml_in, "r") as multiFile:
		for line in multiFile:
			if row==0 and line.startswith(">"):
				print(line, end="", file=converted)
				row+=1
			elif row>0 and line.startswith(">"):
				print("\n" + line, end="", file=converted)
			elif len(line)==0:
				pass
			else:
				print(line.strip(), end="", file=converted),
	converted.close()


def multiprocessing_multiple_al(bin_list):
    jobs = []
    maln_lst = list()
    for bin in bin_list:
        outfile = bin + ".maln"
        maln_lst.append(outfile)
        p = mp.Process(target=mp_multiple_alignment, args=(bin, outfile))
        jobs.append(p)
        p.start()
        
    for proc in jobs:
        proc.join()
    
    with open("pickled_maln.lst", "wb") as mln:
        pickle.dump(maln_lst, mln)
            
    return maln_lst


# def multiple_alignment(bin_list): #single threaded
#     maln_lst = list()
#     for bin in bin_list:
#         outfile = bin + ".maln"
#         maln_lst.append(outfile)
#         #with open(bin, "r") as readpile:
#         subprocess.call(["muscle", "-in", bin, "-out", outfile])
#     for mal in maln_lst:
#         convert_multiline_fasta(mal, mal + ".tmp")
#     for mal in maln_lst:
#         os.remove(mal)
#         os.rename(mal + ".tmp", mal)    
#     return maln_lst


def mp_multiple_alignment(bin, outfile): #multiprocessing
    #with open(bin, "r") as readpile:
    subprocess.call(["muscle", "-in", bin, "-out", outfile])


def gen_consensus(maln_lst):
    consensus_dict = {}
    print(maln_lst)
    for multaln in maln_lst:
        print(multaln)
        aln_lst = list()
        convert_multiline_fasta(multaln, multaln + ".tmp")
        #os.remove(multaln)
        os.rename(multaln, multaln + ".ml")
        os.rename(multaln + ".tmp", multaln) 
        with open(multaln, "r") as mal:
            for line in mal:
                if line.startswith(">"):
                    aln_lst.append(list(next(mal).strip()))
        n = len(aln_lst[0])
        #print(n)
        maln_stack = np.asarray(aln_lst, "str")
        col_occur_a = np.count_nonzero(maln_stack == "A", axis = 0)
        #print("A:\n", col_occur_a)
        col_occur_t = np.count_nonzero(maln_stack == "T", axis = 0)
        col_occur_c = np.count_nonzero(maln_stack == "C", axis = 0)
        col_occur_g = np.count_nonzero(maln_stack == "G", axis = 0)
        col_occur_gap = np.count_nonzero(maln_stack == "-", axis = 0)
        #print("-:\n", col_occur_gap)
        #scoring_array = np.assaray([col_occur_a, col_occur_t, col_occur_c, col_occur_g, col_occur_gap])
        consensus = ""
        for i in range(n):
            score_dict = {"A": col_occur_a[i], "T": col_occur_t[i], "C": col_occur_c[i], "G": col_occur_g[i], "-": col_occur_gap[i]}
            #print(score_dict)
            new = max(score_dict.items(), key=lambda x: x[1])[0]
            #print(new)
            consensus = consensus + new
        print(consensus)
        cl = list(consensus)
        compressed_consensus = "".join([char for char in cl if char != "-"])
        head, tail = os.path.split(multaln)
        amp_no = int(tail.split("_")[0].split("Amp")[1])
        print(amp_no)
        consensus_dict[amp_no] = compressed_consensus #multaln.split("_")[0]
        print(compressed_consensus)
        with open("pickled_consensus.dict", "wb") as coco:
            pickle.dump(consensus_dict, coco)
    return consensus_dict

        
def analyze_als(als):
    lct = 0
    nct = set()
    names = Counter()
    for al in als:
        names[al.qname] += 1
        if al.qlen > 462:
            #print(al.qlen)
            lct += 1
    for n in names:
        if names[n] > 1:
            #print(n, names[n])
            nct.add(n)
    print("lencount=", lct)
    print("namecount=", len(nct))


def arrange_consensus_tiles(consensus_dict, amp_binned_als): #


    ovl_dict = dict()
    genome_len = int()
    for amp in amp_binned_als:
        if type(amp.ovl_with_prev_amp) == int:
            genome_len = genome_len + amp.amp_len - amp.ovl_with_prev_amp
            ovl_dict[amp.amp_name] = amp.ovl_with_prev_amp
        else:
            genome_len = genome_len + amp.amp_len - max(amp.ovl_with_prev_amp)
            ovl_dict[amp.amp_name] = max(amp.ovl_with_prev_amp)
        
    
    amp_obj_dict = {}
    for amp in amp_binned_als:
        amp_obj_dict[amp.amp_name] = amp

    print("Genome length", genome_len)

    cons_tiles = sorted(consensus_dict.items())
    stitched_genome = str()

    
    for tile in cons_tiles:
        # seqment = tile[1]
        # amp_num = tile[0]
        # ovl_idx = ovl_dict[amp_num]
        #print(ovl_idx)
        stitched_genome = stitched_genome + tile[1][ovl_dict[tile[0]]:]

    print(stitched_genome)

    with open("result_stitched_genome", "w") as out:
        out.write(stitched_genome)



if __name__ == "__main__":
    reads_fasta = "/home/simon/smk-artic-ont/smk-artic-ont/results/BC11_corr/BC11.correctedReads.fasta"
    max_amp_len = 420
    amp_csv = "./resources/ARTICv3_amplicons.csv"
    paf = sys.argv[1]
    reads = load_reads(reads_fasta)
    amps = load_amplicons(amp_csv)
    merged_amps = merge_alt__amps(amps)
    als = load_mm2_paf(paf, reads)
    filt_als = filter_paf(als, max_amp_len)
    amp_binned_als = sort_alns_to_amps(merged_amps, filt_als)
    bin_list = write_out_bins(amp_binned_als)
    maln_lst= multiprocessing_multiple_al(bin_list)
    # with open("pickled_maln.lst", "rb") as pm:
    #     maln_lst = pickle.load(pm)
    # with open("pickled_amp_binned_als.lst", "rb") as aba:
    #     amp_binned_als = pickle.load(aba)
    # with open("pickled_consensus.dict", "rb") as pcd:
    #     consensus_dict = pickle.load(pcd)
    consensus_dict = gen_consensus(maln_lst)
    stitched_genome = arrange_consensus_tiles(consensus_dict, amp_binned_als)




    #analyze_als(filt_als)
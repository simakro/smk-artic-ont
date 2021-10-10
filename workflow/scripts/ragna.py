#Reference Assembly Guided de-Novo Assembly#
import sys
import os
from collections import Counter
import subprocess
import numpy as np
    

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
    #print(read_dict)
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
    return merged_amps

def write_out_bins(amp_binned_als):
    bins = list()
    folder_path = "./results/binned_reads/"
    os.makedirs(folder_path, exist_ok=True)
    for amp in amp_binned_als:
        bin_path = folder_path + "Amp" + str(amp.amp_name) + "_binned_reads.fasta"
        bins.append(bin_path)
        with open(bin_path, "w") as out:
            for read in amp.stacked_reads:
                out.write(">" + read + "\n")
                out.write(amp.stacked_reads[read] + "\n")
    return bins


# def multiprocessing_len_analysis(read_dir, fasta_ext):
#         manager = mp.Manager()
#         return_list = manager.list()
#         return_dict = manager.dict()
#         jobs = []

#         for entry in os.scandir(read_dir):
#             file_path = str(entry.path)
#             tail = os.path.split(file_path)[1]
#             f_name, f_ext = "".join(tail.split(".")[:-1]), tail.split(".")[-1]
#             if f_ext in fasta_ext and entry.is_file():
#                 p = mp.Process(target=get_read_len, args=(str(entry.path), return_list, return_dict))
#                 jobs.append(p)
#                 p.start()
        
#         for proc in jobs:
#             proc.join()
            
#         return return_list, return_dict

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


def multiple_alignment(bin_list):
    maln_lst = list()
    for bin in bin_list:
        outfile = bin + ".maln"
        maln_lst.append(outfile)
        with open(bin, "r") as readpile:
            subprocess.call(["muscle", "-in", bin, "-out", outfile])
    for mal in maln_lst:
        convert_multiline_fasta(mal, mal + ".tmp")
    for mal in maln_lst:
        os.remove(mal)
        os.rename(mal + ".tmp", mal)    
    return maln_lst

def gen_consensus(maln_lst):
    for multaln in maln_lst: 
        aln_lst = list()
        with open(multaln, "r") as mal:
            for line in mal:
                if line.startswith(">"):
                    aln_lst.append(list(next(mal).strip()))
        #n = len(aln_lst[0])
        maln_stack = np.asarray(aln_lst, "str")
        col_occur_a = np.count_nonzero(maln_stack == "A", axis = 0)
        #print("A:\n", col_occur_a)
        col_occur_t = np.count_nonzero(maln_stack == "T", axis = 0)
        col_occur_c = np.count_nonzero(maln_stack == "C", axis = 0)
        col_occur_g = np.count_nonzero(maln_stack == "G", axis = 0)
        col_occur_gap = np.count_nonzero(maln_stack == "-", axis = 0)
        #print("-:\n", col_occur_gap)
        cons_profile = list()
        
        # for mal in aln_lst:
        #     local_profile = {"A": 0, "T": 0, "G": 0, "C": 0, "-": 0}
        #     for char in mal:

        #consensus = 
        
        # profile = { 'T':[0]*n,'G':[0]*n ,'C':[0]*n,'A':[0]*n, '-':[0]*n}
        # bestseqs = [[]]
        # for i in range(n):
        #     d = {N:profile[N][i] for N in ['T','G','C','A', '-']}
        #     m = max(d.values())
        #     l = [N for N in ['T','G','C','A'] if d[N] == m]
        #     bestseqs = [ s+[N] for N in l for s in bestseqs ]
        #     print(bestseqs)

        # for s in bestseqs:
        #     print(''.join(s))

        
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



if __name__ == "__main__":

    # from Bio import AlignIO
    # from Bio.Align import AlignInfo
    # alignment = AlignIO.read("/Users/simon/Documents/GitHub/smk-artic-ont/results/binned_reads/Amp2_binned_reads.fasta.maln", "fasta")
    # summary_align = AlignInfo.SummaryInfo(alignment)
    # summary_align.dumb_consensus(0.5)


    # reads_fasta = "/Users/simon/Documents/IKIM_Experiments/SM_0002_GridION_SARS-CoV2_ArticV3_Etablierung/SM0002_GridION_SARSCoV2_ARTICv3_20210915/SM002_barcodes_both_ends_demultiplexed/canu_corrected_and_trimmed/BC02_corr/BC02.correctedReads.fasta"
    # max_amp_len = 420
    # amp_csv = "./resources/ARTICv3_amplicons.csv"
    # paf = sys.argv[1]
    # reads = load_reads(reads_fasta)
    # amps = load_amplicons(amp_csv)
    # merged_amps = merge_alt__amps(amps)
    # als = load_mm2_paf(paf, reads)
    # filt_als = filter_paf(als, max_amp_len)
    # amp_binned_als = sort_alns_to_amps(merged_amps, filt_als)
    # bin_list = write_out_bins(amp_binned_als)
    # maln_lst= multiple_alignment(bin_list)
    gen_consensus(["/Users/simon/Documents/GitHub/smk-artic-ont/results/binned_reads/Amp1_binned_reads.fasta.maln.sl"])



    #analyze_als(filt_als)

import sys
import statistics as stat
import os


def fastq_to_fasta(fastq_in):

    head, tail = os.path.split(fastq_in)
    spt = tail.split(".")
    new_name = "".join([".".join([i for i in spt[:-1]]), ".fasta"])

    def cont_gen():
        while True:
            yield "header"
            yield "seq"
            yield "plus"
            yield "qstring"

    content = cont_gen()

    bases = ["A", "T", "G", "C", "N"]
    with open(fastq_in, "r") as fastq:
        fasta = open(os.path.join(head, new_name), "w")
        # readlen = []
        for line in fastq:
            cur_line = next(content)
            if cur_line == "header":
                if line.startswith("@"):
                    line = line.replace("@", ">")
                    fasta.write(line)
                    # fasta.write(next(fastq))
            elif cur_line == "seq":
                if line[0] in bases:
                    fasta.write(line)
                else:
                    print("Unexpected character in sequence line")
            elif cur_line == "plus":
                if line.startswith("+"):
                    pass
                else:
                    print("Unexpected character in plus line")
            else:
                pass

                # readlen.append(len(next(fastq)))
        fasta.close()
        # print(readlen)
        # print("No. of reads", len(readlen))
        # print("Average read length:", stat.mean(readlen))
        # print("longest read", max(readlen))
        # print("Total base count", sum(readlen))
    return new_name


if __name__ == "__main__":
    fastq_in = snakemake.input[0]
    fastq_to_fasta(fastq_in)



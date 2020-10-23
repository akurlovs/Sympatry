import pandas as pd
import numpy as np
import sys
import argparse
import subprocess
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plot
import matplotlib.patches as patches
from matplotlib import rcParams
from cycler import cycler
import matplotlib.pylab as lab
import glob
from itertools import cycle

rcParams['font.sans-serif'] = 'Arial'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

######################################## DESCRIPTIONS AND PARSER #####################################################

SCRIPT_DESC="""Code estimates heterozygosity across the genome for each individual in a vcf file and plots it as pdfs"""

PARSER=argparse.ArgumentParser(
    description=SCRIPT_DESC,
    formatter_class=argparse.RawDescriptionHelpFormatter)

## IMPORTANT: if only these two are included, only coverage code will be run
PARSER.add_argument("-v", "--vcf", required=True, help="what is the vcf file to analyze?")
PARSER.add_argument("-o", "--outdir", required=True, help="in which directory should we place the files (full path)?")

# these are for the matrix code. the first part asks if the matrix needs to be created
# in the matrix, 0=homo, 1=het, 2=missing call. for the rest, see below:
PARSER.add_argument("-p", "--plot", required=False, action="store_true", help="IMPORTANT:if you put anything here, only the plotting will take place")
PARSER.add_argument("-m", "--matrix", required=False, action="store_true", help="IMPORTANT:if you put anything here, the matrix code will run")
PARSER.add_argument("-k", "--skip", required=False, action="store_true", help="IMPORTANT: if you put anything here, the program will assume that the coverage file has already been made")
PARSER.add_argument("-r", "--filter_vcf", required=False, action="store_true", help="IMPORTANT: if you put anything here, the program will output a filtered vcf")
PARSER.add_argument("-c", "--covdiff", required=False, default=0.25, help="coverage deviation from average that's allowed the default is 0.25 meaning everything above 0.75*mean_coverage and everything below 1.25*mean_coverage will get a pass. fail=3 in the matrix")
PARSER.add_argument("-qds", "--qds", required=False, default=2,
                    help="Mininum variant score")
PARSER.add_argument("-sor", "--sor", required=False, default=3,
                    help="Maximum strand bias score")
PARSER.add_argument("-mq", "--mps", required=False, default=50,
                    help="Minimum mean mapping quality score")
PARSER.add_argument("-mqrs", "--mqrs", required=False, default=8,
                    help="Minimum mean quality rank sum score")
PARSER.add_argument("-rprs", "--rprs", required=False, default=8,
                    help="Minimum read pos rank sum")

# IMPORTANT: if you provide a scaffold info file, the program will assume you want heterozygosity estimates
# if you only provide sin, but do not fill matrix or skip, the program will perform the estimates only.
# if you also provide plot, however, only the plotting will be done
PARSER.add_argument("-s", "--sin", required=False, help="what is the scaffold info file?")
PARSER.add_argument("-w", "--window", required=False, default=100000, help="what is the sliding window size?")
PARSER.add_argument("-l", "--slide", required=False, default=10000, help="by how much will the window slide?")
PARSER.add_argument("-a", "--minallele", required=False, default=20, help="what is minimum number of PASS (0 or 1) alleles in window for the results to be outputted?")
PARSER.add_argument("-f", "--minfreqpass", required=False, default=0.05, help="what is the fraction of PASS (0 or 1) alleles in window for the window to be outputted?")

# This is for plotting purposes
PARSER.add_argument("-i", "--min", required=False, default=-0.02, help ="ylim on plots")
PARSER.add_argument("-x", "--max", required=False, default=1.02, help ="ymax on plots")
PARSER.add_argument("-xx", "--xx", required=False, default=2, help ="x-axis column")
PARSER.add_argument("-yy", "--yy", required=False, default=3, help ="y-axis column")
PARSER.add_argument("-n", "--minor_output", required=False, action="store_true", help ="plot and output minor allele frequency")

##########################################################################################################

ARGIES = PARSER.parse_args()
ARGDICT = {}

ARGDICT["vcf"] = ARGIES.vcf
OUTDIR = ARGIES.outdir
ARGDICT["outdir1"] = OUTDIR+"/info_files"
ARGDICT["outdir2"] = OUTDIR+"/het_files"
ARGDICT["outdir3"] = OUTDIR+"/het_plots"

ARGDICT["plotty"] = ARGIES.plot
ARGDICT["matrix"] = ARGIES.matrix
ARGDICT["skip"] = ARGIES.skip
ARGDICT["filter_vcf"] = ARGIES.filter_vcf
ARGDICT["covdiff"] = float(ARGIES.covdiff)

ARGDICT["qds"] = float(ARGIES.qds)
ARGDICT["mps"] = float(ARGIES.mps)
ARGDICT["sor"] = float(ARGIES.sor)
ARGDICT["mqrs"] = float(ARGIES.mqrs)
ARGDICT["rprs"] = float(ARGIES.rprs)

ARGDICT["sin"] = ARGIES.sin
ARGDICT["window"] = int(ARGIES.window)
ARGDICT["slide"] = int(ARGIES.slide)
ARGDICT["minallele"] = int(ARGIES.minallele)
ARGDICT["minfreqpass"] = float(ARGIES.minfreqpass)

ARGDICT["xx"] = int(ARGIES.xx)
ARGDICT["yy"]= int(ARGIES.yy)
ARGDICT["ylim_low"] = float(ARGIES.min)
ARGDICT["ylim_high"] = float(ARGIES.max)
ARGDICT["minor_output"] = ARGIES.minor_output

########################################### FUNCTIONS ##########################################################
# COVERAGE
def coverage():
    cov = {}
    subprocess.call("mkdir %s"%(ARGDICT["outdir1"]), shell=True)
    outcov = open(ARGDICT["outdir1"]+"/coverageinfo.txt", "w")
    outcov.write("%s\t%s\t%s\n"%("#STRAIN", "AVERAGE_DEPTH", "SD_DEPTH"))
    with open(ARGDICT["vcf"], "r") as openvcf:
        for line in openvcf:
            if line[0:2] == "##": #this contains the lengths of the scaffold for the original pass I don't care about it.
                pass
            elif line[0:6] == "#CHROM": #need to store all the variable names for the samples and which columns correspond to it.
                header = line.strip().split("\t")
                header_strains = header[9:]
                for strain in header_strains:
                    cov[strain] = {}
                    cov[strain]["count"] = 0.0
                    cov[strain]["covva"] = 0.0
                    cov[strain]["all_cov"] = []
            else:
                if line[0] != "#": #need average coverage for parent1, parent2, sel, and unsel
                    SNP_line = line.strip().split("\t")
                    ref_allele = SNP_line[3]
                    alt_alleles = SNP_line[4].split(",")
                    alt_snp = 1
                    for alt in alt_alleles:
                        if len(alt) > 1:
                            alt_snp = len(alt)
                    if len(ref_allele) == 1 and alt_snp == 1:
                        haplotypes = SNP_line[9:]
                        index = 0
                        for strain in header_strains:
                            choice = haplotypes[index]
                            if "./" not in choice:
                                reads = choice.split(":")[1].split(",")
                                parent_cov = sum([int(i) for i in reads])
                                cov[strain]["count"] += 1
                                cov[strain]["covva"] += parent_cov
                                cov[strain]["all_cov"].append(parent_cov)
                                index += 1
                            else:
                                index += 1
    for strain in header_strains:
        avecovvcf = str(cov[strain]["covva"]/cov[strain]["count"])
        stdcovvcf = np.std(cov[strain]["all_cov"])
        outcov.write("%s\t%s\t%s\n"%(strain, avecovvcf, stdcovvcf))
    outcov.close()

def line_parser(line):
    """Parses VCF string to get mapping quality info"""
    linedict = {}
    line = line.split(";")
    for keyval in line:
        keyval = keyval.split("=")
        if len(keyval) > 1 and "|" not in keyval[1]:
            values = keyval[1].split(",")
            if len(values) == 1:
                linedict[keyval[0]] = float(keyval[1])
            else:
                linedict[keyval[0]] = min([float(i) for i in values])
    return(linedict)

# MATRIX
def get_matrix():
    cov = {}
    with open(COVERAGEFILE, "r") as opencov:
        for line in opencov:
            if line[0] != "#":
                linestrip = (line.rstrip()).split("\t")
                st = linestrip[0]
                co = float(linestrip[1])
                cov[st] = co
    outfile = open(ARGDICT["outdir1"]+"/matrix.txt", "w")
    if ARGDICT["filter_vcf"]:
        outvcf = open(ARGDICT["outdir1"]+"/filtered.vcf", "w")
    with open(ARGDICT["vcf"], "r") as openvcf:
        for line in openvcf:
            if line[0:2] == "##": #this contains the lengths of the scaffold for the original pass I don't care about it.
                if ARGDICT["filter_vcf"]:
                    outvcf.write(line)
            elif line[0:6] == "#CHROM": #need to store all the variable names for the samples and which columns correspond to it.
                header = line.strip().split("\t")
                header_strains = header[9:]
                outfile.write("\t".join(["chrom", "pos"]+header_strains)+"\n")
                if ARGDICT["filter_vcf"]:
                    outvcf.write(line)
            else:
                if line[0] != "#": #need average coverage for parent1, parent2, sel, and unsel
                    SNP_line = line.strip().split("\t")
                    SNP_scaff = SNP_line[0]
                    SNP_pos = SNP_line[1]
                    ref_allele = SNP_line[3]
                    alt_alleles = SNP_line[4].split(",")
                    alt_snp = 1
                    for alt in alt_alleles:
                        if len(alt) > 1:
                            alt_snp = len(alt)
                    if len(ref_allele) == 1 and alt_snp == 1:
                        haplotypes = SNP_line[9:]
                        scores = [SNP_scaff, SNP_pos]
                        index = 0
                        for choice in haplotypes:
                            # is the variant there. 2 = not there
                            if "./" in choice:
                                scores.append("2")
                                index = index + 1
                                continue
                            # is the COVERAGEFILE okay?
                            reads = choice.split(":")[1].split(",")
                            covey = sum([int(i) for i in reads])
                            strain = header_strains[index]
                            avecov = cov[strain]
                            if covey < avecov*(1-ARGDICT["covdiff"]) or covey > avecov*(1+ARGDICT["covdiff"]):
                                scores.append("3")
                                index = index +1
                                continue
                            # this is the read quality filter
                            stuff = SNP_line[7]
                            quals = line_parser(stuff)
                            if ("QD" in quals and "MQ" in quals and "SOR" in quals
                                    and "MQRankSum" in quals and "ReadPosRankSum" in quals):
                                if (quals["QD"] < ARGDICT["qds"]
                                        or quals["MQ"] < ARGDICT["mps"]
                                        or quals["SOR"] >= ARGDICT["sor"]
                                        or quals["MQRankSum"] < -ARGDICT["mqrs"]
                                        or quals["MQRankSum"] > ARGDICT["mqrs"]
                                        or quals["ReadPosRankSum"] < -ARGDICT["rprs"]
                                        or quals["ReadPosRankSum"] > ARGDICT["rprs"]):
                                    scores.append("4")
                                    index = index +1
                                    continue
                            elif ("QD" in quals and "MQ" in quals and "SOR" in quals):
                                if (quals["QD"] < ARGDICT["qds"]
                                        or quals["MQ"] < ARGDICT["mps"]
                                        or quals["SOR"] >= ARGDICT["sor"]):
                                    scores.append("4")
                                    index = index +1
                                    continue 
                            else:
                                scores.append("5")
                                index = index +1
                                continue
                             # now the genotype pass
                            genotype = choice.split(":")[0]
                            numeric_geno = set([int(i) for i in genotype.split("/")])
                            alleles = genotype.split("/")
                            all_uniques = len(list(set(alleles)))
                            if all_uniques == 1:
                                reforalt = int(alleles[0])
                                if reforalt == 0:
                                    scores.append("10")
                                    index = index +1
                                else:
                                    scores.append("-1")
                                    index = index +1
                            else:
                                select_reads = [int(reads[i]) for i in numeric_geno]
                                minor = min(select_reads)/float(sum(select_reads))
                                scores.append("%s"%(minor))
                                index = index + 1
                        outfile.write("\t".join(scores)+"\n")
                        if ARGDICT["filter_vcf"]:
                            passing_scores = [i for i in scores[2:] if
                                              0 <= float(i) <= 1 or int(i) == -1 or int(i) == 10]
                            if len(passing_scores) == len(scores[2:]):
                                outvcf.write(line)
    if ARGDICT["filter_vcf"]:
        outvcf.close()
    outfile.close()

# this is to get the correct overall genomic position
def scale_dict(info_file):
    info_read = open(info_file, "r")
    info_dict = {}
    final_lenny = 0
    lenny = 0
    for line in info_read:
        linetab = (line.rstrip()).split("\t")
        scaffy = linetab[0]
        final_lenny = final_lenny + lenny
        lenny = int(linetab[1])
        info_dict[scaffy] = final_lenny
    info_read.close()
    return(info_dict)

## SLIDING WINDOW FOR HETEROZYGOSITY ACROSS THE GENOMES
def homo():
    subprocess.call("mkdir %s"%(ARGDICT["outdir2"]), shell=True)
    snps = pd.read_table(MATRIX_FILE, sep="\t")
    ### estimates for entire genome first ###
    outhetstat = open(ARGDICT["outdir1"] + "/" + "allhet.txt", "w")
    outhetstat.write("strain"+"\t"+"Genome_wide_heterozygosity"+"\t"+"Genome_wide_minor_freq"+"\t"+"Fraction_windows_with_het>0.1"+"\n")
    outhetdict = {}
    for col in snps.columns[2:]:
        #print("col", col)
        section = snps[col]
        #print(section, "section")
        homos = len(section[section==-1])
        heteros = len(section[(section >= 0) & (section <= 1)])
        if heteros > 0:
            minorfreq = np.mean(section[(section >= 0) & (section <= 1)])
        else:
            minorfreq = 0
        npass = homos+heteros
        #print("Npass", npass)
        passer = npass/float(len(section))
        heterofix = heteros/float(npass)
        strainy = str(col)
        print(strainy, minorfreq)
        outhetdict[strainy] = {}
        outhetdict[strainy]["genome_wide_het"] =  heterofix
        outhetdict[strainy]["genome_wide_minor"] =  minorfreq
        outhetdict[strainy]["n_windows"] = 0
        outhetdict[strainy]["high_het"] = 0
        outfile = open(ARGDICT["outdir2"] + "/" + col + ".txt", "w") # to erase previous files if needed
        outfile.close()
    ### end of estimate for entire genome ###
    scaffolds = snps['chrom'].unique()
    for scaff in scaffolds:
        if (scaff in SHADER) == False:
            continue
        scaff_seg = snps[snps['chrom']==scaff]
        beg = 0
        win = ARGDICT["window"]
        end = ARGDICT["window"]
        slide = ARGDICT["slide"]
        if max(scaff_seg['pos']) < win:
            continue
        while end < max(scaff_seg['pos']):
            chunk = scaff_seg[(scaff_seg['pos'] > beg) & (scaff_seg['pos'] < end)]
            ntotal = len(chunk)
            med = beg + win/2
            medpos = med + SHADER[scaff]
            for col in chunk.columns[2:]:
                #print("col", col)
                section = chunk[col]
                #print(section, "section")
                homos = len(section[section==-1])
                heteros = len(section[(section >= 0) & (section <= 1)])
                if heteros > 0:
                    minorfreq = np.mean(section[(section >= 0) & (section <= 1)])
                else:
                    minorfreq = 0
                npass = homos+heteros
                if len(section) > 0:
                    passer = npass/float(len(section))
                else:
                    passer = 0
                if npass >= ARGDICT["minallele"] and passer >= ARGDICT["minfreqpass"]:
                    heterofix = heteros/float(npass)
                    strainy = str(col)
                    outhetdict[strainy]["n_windows"] = outhetdict[strainy]["n_windows"]+1
                    if heterofix > 0.1:
                        outhetdict[strainy]["high_het"] = outhetdict[strainy]["high_het"] + 1
                    outfile = open(ARGDICT["outdir2"] + "/" + col + ".txt", "a")
                    outfile.write("\t".join([scaff, str(med), str(medpos), str(heterofix), str(minorfreq), str(npass), str(ntotal)])+"\n")
                    outfile.close()
            beg = beg + slide
            end = end + slide
    for strain_dict in sorted(outhetdict):
        high_het_count = str(outhetdict[strain_dict]["high_het"]/float(outhetdict[strain_dict]["n_windows"]))
        outhetstat.write(strain_dict + "\t" + str(outhetdict[strain_dict]["genome_wide_het"]) + "\t" +  str(outhetdict[strain_dict]["genome_wide_minor"]) + "\t" + high_het_count + "\n")
    outhetstat.close()

## PLOTTING
# alternating colors
def plotting():
    subprocess.call("mkdir %s"%(ARGDICT["outdir3"]), shell=True)
    alterna = ['green', 'blue', 'violet', 'purple', 'royalblue', 'orange', 'fuchsia', 'pink']
    plot.rc('axes', prop_cycle=(cycler('color', alterna)))
    fileys = glob.glob(ARGDICT["outdir2"]+"/*.txt")
    # now the plotting
    print("***********************")
    pool = cycle(alterna)
    for filey in fileys:
        fig = plot.figure(figsize=(15, 5), dpi=10000)
        typies = filey.split("/")[-1].split(".")[0]
        tsvet = next(pool)
        print("%s : %s" %(tsvet, filey))
        b = pd.read_table(filey, sep="\t", header=None)
        if ARGDICT["minor_output"]:
             plot.plot(b.iloc[0:,ARGDICT["xx"]],b.iloc[0:,ARGDICT["yy"]+1], lw=1, color="0.75")
        plot.plot(b.iloc[0:,ARGDICT["xx"]],b.iloc[0:,ARGDICT["yy"]], lw=2, color=tsvet)
        # ACCESSORY STUFF TO MAKE IT PRETTY
        # how far will the plot go
        lab.ylim(ARGDICT["ylim_low"], ARGDICT["ylim_high"])
        maxx = max(b.iloc[0:,ARGDICT["xx"]])
        lab.xlim(0, maxx)
        # this are the variables for the shading below
        old = None
        shade = False
        # here comes the shading
        for contig in sorted(SHADER):
            val = SHADER[contig]
            if old != None and shade == True:
                plot.axvspan(old, val, color='0.87', alpha=0.5)
                shade = False
            else:
                if old != None:
                    shade = True
            old = SHADER[contig]
        # the last one
        if SHADER == True:
            plot.axvspan(old, maxx, color='0.87', alpha=0.5)
        # grid and labels
        plot.locator_params(axis='x',nbins=15)
        plot.grid(True, axis="y")
        # for better labeling:
        z = [int(i) for i in plot.xticks()[0]][1:-1] # this is to create labels with numbers
        labs = []
        for i in plot.xticks()[0]:
            j = int(i)
            if len(str(j)) <= 2:
                labs.append(str(round(int(j)/float(1), 3)))
            elif 3 <= len(str(j)) <= 6:
                labs.append(str(round(int(j)/float(1000), 3))+" Kb")
            elif  3 <= len(str(j)) <= 9:
                labs.append(str(int(int(j)/float(1000000))))
            else:
                labs.append(str(round(int(j)/float(1000000000), 1))+" Gb")
        labs = labs[1:-1]
        plot.xticks(z, labs)
        plot.xlabel("Genomic Position (Mb)", fontsize = 22)
        plot.ylabel("Fraction of heterozygous SNPs", fontsize = 22)
        plot.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='on',         # ticks along the top edge are off
            labelbottom='on',
            labelsize = 22) # labels along the bottom edge are off
        plot.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            right='off',      # ticks along the bottom edge are off
            left='on',         # ticks along the top edge are on
            labelbottom='off', # labels along the bottom edge are off
            labelsize = 22)        # changes label size to whatevs
        fig.savefig(ARGDICT["outdir3"]+"/"+typies+".pdf", bbox_inches='tight')
        fig.clf()
    plot.close('all')
    print("***********************")

############################################ THIS PART OF THE SCRIPT ACTUALLY RUNS THE FUNCTIONS ###################

if ARGDICT["plotty"]:
    print("plotting every file")
    SHADER = scale_dict(ARGDICT["sin"])
    plotting()
else:
    if not ARGDICT["sin"]:
        if not ARGDICT["matrix"]:
            # creates the coverage info file
            print("creating coverage info file")
            coverage()
        else:
            if ARGDICT["skip"]:
                COVERAGEFILE = ARGDICT["outdir1"]+"/coverageinfo.txt"
                print("creating matrix")
                get_matrix()
            else:
                # creates the coverage info file
                print("creating coverage info file")
                coverage()
                COVERAGEFILE = ARGDICT["outdir1"]+"/coverageinfo.txt"
                # creates the matrix
                print("creating matrix")
                get_matrix()
    
    else:
        if not ARGDICT["matrix"]: 
            print("performing heterozygosity estimates for each individual")
            MATRIX_FILE = ARGDICT["outdir1"]+"/"+"matrix.txt"
            # process the scaffold information
            SHADER = scale_dict(ARGDICT["sin"])
            homo()
            # plots heterozygosity 
            print("plotting every file")
            plotting()
        else:
            if not ARGDICT["skip"]:
                # creates the coverage info file
                print("creating coverage info file")
                coverage()
            COVERAGEFILE = ARGDICT["outdir1"]+"/coverageinfo.txt"
            print("creating matrix")
            get_matrix()
            MATRIX_FILE = ARGDICT["outdir1"]+"/"+"matrix.txt"
            # processing the scaffold information
            SHADER = scale_dict(ARGDICT["sin"])
            print("performing heterozygosity estimates for each individual")
            homo()
            # plots heterozygosity 
            print("plotting every file")
            plotting()
        
############################################ THE END ##############################################
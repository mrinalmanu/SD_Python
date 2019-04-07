import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
from numpy import random
import os.path
import sys


class Metrics(object):
    """Container of various metrics.
    """
    def __init__(self):
        """Initialize the metrics.
        """
        self.replicates = 0
        self.variant_positions = 0
        self.monomorphic_variant_positions = 0
        self.polymorphic_variant_positions = 0
        self.single_replicate_positions = 0
        self.multiple_replicate_positions = 0


def read_fasta_sequence(fasta_file_path):
    """
    Read a fasta file and return the sequence as a string.
    """
    seq_string = []
    with open(fasta_file_path, "r") as infile:
        for idx, line in enumerate(infile):
            if idx == 0:
                seq_name = line.strip(">").strip("\n")
            elif ">" in line:
                pass
            else:
                seq_string.append(line.strip("\n"))

    seq_str = "".join(seq_string)
    return (seq_name, seq_str)


def write_fasta_sequence(seq_id, file_path, sequence_list, mutations):
    """
    Write a mutated sequence to a fasta file.  The fasta defline is suffixed
    with a mutation description of the form (mutated s=100 i=20 d=20)
    describing the number of substitutions, insertions, and deletions.
    """
    seq_str = "".join(sequence_list)
    seq = Seq(seq_str)
    description = "(mutated s=%i i=%i d=%i)" % mutations
    record = SeqRecord(seq, id=seq_id, description=description)
    SeqIO.write([record], file_path, "fasta")


def get_all_eligible_positions(seq_str):
    """Find all the positions where the original base is ACTG.  Those are the only
    positions eligible for mutation.
    """
    eligible_positions = list()
    pos = 0
    for base in seq_str:
        if base.upper() in "ACGT":
            eligible_positions.append(pos)
        pos += 1
    return eligible_positions


def build_mutated_seq(seq_str, eligible_positions, num_subs, num_insertions, num_deletions):
    """
    Copy a sequence and randomly introduce the specified numbers of
    substitutions, insertions, and deletions.

    """
    substitution_choices = {"A": ["C", "T", "G"],
                            "C": ["A", "T", "G"],
                            "T": ["C", "A", "G"],
                            "G": ["C", "T", "A"],
                            }

    num_mutations = num_subs + num_deletions + num_insertions
    positions = random.choice(eligible_positions, size=num_mutations, replace=False)
    subs_positions = positions[: num_subs]
    deletion_positions = positions[num_subs: (num_subs + num_deletions)]
    insertion_positions = positions[(num_subs + num_deletions): len(positions)]

    # Copy the original sequence in a way that easily allows mutations
    # while preserving position information
    new_indexed_seq = list(seq_str)

    for i in subs_positions:
        original_base = seq_str[i]
        upper_original_base = original_base.upper()
        if upper_original_base in substitution_choices:
            new_base = random.choice(substitution_choices[upper_original_base], size=1)[0]
            new_indexed_seq[i] = new_base
        else:
            print('Warning: unexpected base "%s" at position %i.  There will be no mutation at this position.' % (original_base, i), file=sys.stderr)

    for i in deletion_positions:
        new_indexed_seq[i] = ""

    for i in insertion_positions:
        original_base = seq_str[i]
        insert_base = random.choice(["A", "C", "T", "G"], size=1)[0]
        new_indexed_seq[i] = original_base + insert_base

    return (new_indexed_seq, subs_positions, insertion_positions, deletion_positions, )


def mutate_all(seq_str, eligible_pos, num_subs, num_insertions, num_deletions):
    """
    Generate a set mutation for each position in eligible positions.

    """
    print("Creating monomorphic mutations in %d positions." % len(eligible_pos), file=sys.stderr)

    substitution_choices = {"A": ["C", "T", "G"],
                            "C": ["A", "T", "G"],
                            "T": ["C", "A", "G"],
                            "G": ["C", "T", "A"],
                            }
    sub_d = {}
    ins_d = {}
    del_d = set()

    num_mutations = float(num_subs + num_insertions + num_deletions)
    sub_len = round(num_subs / num_mutations * len(eligible_pos))
    ins_len = round(num_insertions / num_mutations * len(eligible_pos))
    del_len = round(num_deletions / num_mutations * len(eligible_pos))

    for x in eligible_pos:
        original_base = seq_str[x]
        if len(sub_d) < sub_len:
            upper_original_base = original_base.upper()
            sub_d[x] = random.choice(substitution_choices[upper_original_base])
        elif len(del_d) < del_len:
                del_d.add(x)
        else:
            insert_base = random.choice(["A", "C", "T", "G"])
            ins_d[x] = original_base + insert_base

    return (sub_d, ins_d, del_d)


def build_limited_seq(seq_str, eligible_positions, pre_mutated_sub, pre_mutated_ins, pre_mutated_del, num_subs, num_insertions, num_deletions):
    """
    Generate a monomorphic allele sequence.

    """

    # Copy the original sequence in a way that easily allows mutations
    # while preserving position information
    new_indexed_seq = list(seq_str)

    subs_positions = []
    deletion_positions = []
    insertion_positions = []

    if num_subs > 0:
        subs_positions = random.choice(list(pre_mutated_sub.keys()), size=num_subs, replace=False)
        for i in subs_positions:
            new_indexed_seq[i] = pre_mutated_sub[i]

    if num_deletions > 0:
        deletion_positions = random.choice(list(pre_mutated_del), size=num_deletions, replace=False)
        for i in deletion_positions:
            new_indexed_seq[i] = ""

    if num_insertions > 0:
        insertion_positions = random.choice(list(pre_mutated_ins.keys()), size=num_insertions, replace=False)
        for i in insertion_positions:
            new_indexed_seq[i] = pre_mutated_ins[i]

    return (new_indexed_seq, subs_positions, insertion_positions, deletion_positions, )


def run_simulations(seq_str, all_eligible_positions, base_file_name, seq_name, num_sims, num_subs, num_insertions, num_deletions, pool_size, group_size, mono, summary_file_path=None, vcf_file_path=None):
    """Generate multiple random mutations of a reference sequence, repeatedly
    calling build_mutated_seq() to create each of the mutated sequences.
    """
    # Metrics containers
    position_variants = defaultdict(set) # key=pos value=set of variants
    position_replicates = defaultdict(set) # key=pos value=set of replicates

    if summary_file_path:
        snp_list_file = open(summary_file_path, "w")

    if vcf_file_path:
        vcf_writer_obj =  VcfWriter(seq_str, vcf_file_path)

    try:
        if summary_file_path:
            snp_list_file.write("Replicate\tPosition\tOriginalBase\tNewBase\n")

        # When there is no pooling, we create monomorphic mutations only once
        if pool_size == 0:
            eligible_positions = all_eligible_positions
            if mono:
                pre_mutated_sub, pre_mutated_ins, pre_mutated_del = mutate_all(seq_str, eligible_positions, num_subs, num_insertions, num_deletions)

        for replicate in range(1, num_sims + 1):
            # Create a new pool after each group of group_size replicates
            need_new_pool = pool_size > 0 and (replicate - 1) % group_size == 0
            if need_new_pool:
                print("Creating pool of %d positions." % pool_size, file=sys.stderr)
                eligible_positions = random.choice(all_eligible_positions, pool_size, replace=False)
                if mono:
                    pre_mutated_sub, pre_mutated_ins, pre_mutated_del = mutate_all(seq_str, eligible_positions, num_subs, num_insertions, num_deletions)

            print("Creating replicate %i" % replicate, file=sys.stderr)

            replicate_name = base_file_name + "_mutated_" + str(replicate)

            # Mutate
            if not mono:
                new_indexed_seq, subs_positions, insertion_positions, deletion_positions = \
                    build_mutated_seq(seq_str, eligible_positions, num_subs, num_insertions, num_deletions)
            else:
                new_indexed_seq, subs_positions, insertion_positions, deletion_positions = \
                    build_limited_seq(seq_str, eligible_positions, pre_mutated_sub, pre_mutated_ins, pre_mutated_del, num_subs, num_insertions, num_deletions)

            mutations = (num_subs, num_insertions, num_deletions)
            write_fasta_sequence(seq_name, replicate_name + ".fasta", new_indexed_seq, mutations)

            # Collect metrics
            for pos in subs_positions:
                position_variants[pos].add(new_indexed_seq[pos])
                position_replicates[pos].add(replicate_name)
            for pos in insertion_positions:
                position_variants[pos].add(new_indexed_seq[pos])
                position_replicates[pos].add(replicate_name)
            for pos in deletion_positions:
                position_variants[pos].add(new_indexed_seq[pos])
                position_replicates[pos].add(replicate_name)

            # Write summary file
            if summary_file_path:
                summary_list = list()
                summary_list.extend([(pos, new_indexed_seq[pos]) for pos in subs_positions])
                summary_list.extend([(pos, new_indexed_seq[pos] + "_insertion") for pos in insertion_positions])
                summary_list.extend([(pos, "_deletion") for pos in deletion_positions])
                for pos, change in sorted(summary_list):
                    snp_list_file.write("%s\t%i\t%s\t%s\n" % (replicate_name, pos + 1, seq_str[pos], change))

            if vcf_file_path:
                vcf_writer_obj.store_replicate_mutations(replicate_name, new_indexed_seq, subs_positions, insertion_positions, deletion_positions)

        # Prepare metrics
        metrics = Metrics()
        metrics.replicates = num_sims
        metrics.variant_positions = len(position_variants)
        for pos in position_variants:
            if len(position_variants[pos]) == 1:
                metrics.monomorphic_variant_positions += 1
            else:
                metrics.polymorphic_variant_positions += 1
        for pos in position_replicates:
            if len(position_replicates[pos]) == 1:
                metrics.single_replicate_positions += 1
            else:
                metrics.multiple_replicate_positions += 1
        return metrics

    finally:
        if summary_file_path:
            snp_list_file.close()
        if vcf_file_path:
            vcf_writer_obj.write(base_file_name)


def parse_arguments(system_args):
    """Parse command line arguments.
    """
    usage = """I AM MUTATOR KING. I MUTATE. I HAVE GREAT POWER.
    Remember:
    
    "with great power comes great responsibility"-- Uncle Ben"""

    def non_negative_int(value):
        try:
            ivalue = int(value)
        except:
            raise argparse.ArgumentTypeError("Must be a number >= 0")
        if ivalue < 0:
            raise argparse.ArgumentTypeError("Must be >= 0")
        return ivalue

    def positive_int(value):
        try:
            ivalue = int(value)
        except:
            raise argparse.ArgumentTypeError("Must be a number greater than 0")
        if ivalue <= 0:
            raise argparse.ArgumentTypeError("Must be greater than 0")
        return ivalue

    parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(                                             dest="input_fasta_file", type=str,                            help="Input fasta file.")
    parser.add_argument("-o", "--summary",           metavar="FILE", dest="summary_file",     type=str,              default=None, help="Output positional summary file.")
    parser.add_argument("-n", "--num-simulations",   metavar="INT",  dest="num_sims",         type=positive_int,     default=100,  help="Number of mutated sequences to generate.")
    parser.add_argument("-s", "--num-substitutions", metavar="INT",  dest="num_subs",         type=non_negative_int, default=500,  help="Number of substitutions.")
    parser.add_argument("-i", "--num-insertions",    metavar="INT",  dest="num_insertions",   type=non_negative_int, default=20,   help="Number of insertions.")
    parser.add_argument("-d", "--num-deletions",     metavar="INT",  dest="num_deletions",    type=non_negative_int, default=20,   help="Number of deletions.")
    parser.add_argument("-r", "--random-seed",       metavar="INT",  dest="random_seed",      type=int,              default=None, help="Random number seed; if not set, the results are not reproducible.")
    parser.add_argument("-p", "--pool",              metavar="INT",  dest="subset_len",       type=positive_int,     default=0,    help="Choose variants from a pool of eligible positions of the specified size")
    parser.add_argument("-g", "--group",             metavar="INT",  dest="group_size",       type=positive_int,     default=None, help="Group size. When greater than zero, this parameter chooses a new pool of positions for each group of replicates.")
    parser.add_argument('-m', '--mono',         action='store_true', dest="mono",                                                  help="Create monomorphic alleles")
    parser.add_argument("-v", "--vcf",               metavar="FILE", dest="vcf_file",         type=str,              default=None, help="Output VCF file.")
    parser.add_argument("-M", "--metrics",           metavar="FILE", dest="metrics_file",     type=str,              default=None, help="Output metrics file.")

    args = parser.parse_args(system_args)
    return args


def run_from_args(args):
    """Generate multiple random mutations of a reference sequence.
    """
    if args.group_size is None:
        args.group_size = args.num_sims  # just one big group

    # Input file arg
    in_file_path = args.input_fasta_file
    in_file_name = os.path.basename(in_file_path)
    base_file_name, in_file_ext = os.path.splitext(in_file_name)

    # Random seed option
    random.seed(args.random_seed)

    # Read the reference and generate mutations
    seq_name, seq_str = read_fasta_sequence(in_file_path)

    # Find the eligible positions
    all_eligible_positions = get_all_eligible_positions(seq_str)
    if args.subset_len > 0:
        eligible_seq_length = args.subset_len
    else:
        eligible_seq_length = len(all_eligible_positions)

    num_mutations = args.num_subs + args.num_insertions + args.num_deletions
    if num_mutations > eligible_seq_length:
        print("ERROR: You have specified a number of substitutions that is greater than the eligible length of the sequence", file=sys.stderr)
        sys.exit(1)

    metrics = run_simulations(seq_str, all_eligible_positions, base_file_name, seq_name, args.num_sims, args.num_subs, args.num_insertions,
                              args.num_deletions, args.subset_len, args.group_size, args.mono, args.summary_file, args.vcf_file)

    if args.metrics_file:
        with open(args.metrics_file, 'w') as f:
            f.write("%d replicates\n" % metrics.replicates)
            f.write("%d positions having any variants\n" % metrics.variant_positions)
            f.write("%d positions having monomorphic variants\n" % metrics.monomorphic_variant_positions)
            f.write("%d positions having polymorphic variants\n" % metrics.polymorphic_variant_positions)
            f.write("%d variant positions found in exactly one replicate\n" % metrics.single_replicate_positions)
            f.write("%d variant positions found in multiple replicates\n" % metrics.multiple_replicate_positions)


def run_from_line(line):
    """Run a command with a command line.
    """
    argv = line.split()
    args = parse_arguments(argv)
    return run_from_args(args)


def main():
    """This is the main function which is magically turned into an executable
    """
    args = parse_arguments(sys.argv[1:])
    return run_from_args(args)

"""
VCF file writer.
"""

class VcfWriter(object):
    """Class to write VCF files for the mutated replicates.
    """

    def __init__(self, original_seq, file_path):
        """Initialize the VCF writer.
        """
        self.reference = original_seq
        self.file_path = file_path
        self.replicate_names = list()
        self.known_alts = defaultdict(list) # key=pos value=list of alts
        self.alt_dict = dict()              # key=(pos, replicate) value=alt num

    def store_replicate_mutations(self, replicate_name, new_indexed_seq, subs_positions, insertion_positions, deletion_positions):
        """Store all of the mutations for a specified replicate.
        """
        self.replicate_names.append(replicate_name)

        # Store the snps
        for pos in subs_positions:
            alt = new_indexed_seq[pos]
            if alt not in self.known_alts[pos]:
                self.known_alts[pos].append(alt)
            alt_num = 1 + self.known_alts[pos].index(alt) # ALT zero is reference match
            self.alt_dict[(pos, replicate_name)] = alt_num

        # Store the insertions
        for pos in insertion_positions:
            alt = new_indexed_seq[pos]
            if alt not in self.known_alts[pos]:
                self.known_alts[pos].append(alt)
            alt_num = 1 + self.known_alts[pos].index(alt) # ALT zero is reference match
            self.alt_dict[(pos, replicate_name)] = alt_num

        # Store the deletions
        for pos in deletion_positions:
            alt = "*"
            if alt not in self.known_alts[pos]:
                self.known_alts[pos].append(alt)
            alt_num = 1 + self.known_alts[pos].index(alt) # ALT zero is reference match
            self.alt_dict[(pos, replicate_name)] = alt_num

    def write(self, reference_name):
        """Write the stored mutations to the VCF output file.

        Parameters
        ----------
        reference_name : str
            Name of reference to write to the header.
        replicate_names : list of str
            List of sample names.
        """
        print("Writing VCF file.", file=sys.stderr)
        with open(self.file_path, "w") as f:
            # Write the VCF file header sections.
            f.write('##fileformat=VCFv4.2\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##reference=%s\n' % reference_name)
            header = '#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT %s\n' % ' '.join(self.replicate_names)
            header = header.replace(' ', '\t')
            f.write(header)

            # Write the variant calls
            chrom = "1" # TODO: use the real chrome name
            for pos in sorted(self.known_alts):
                ref = self.reference[pos]
                alt_str = ','.join(self.known_alts[pos])
                fields = [chrom, str(pos + 1), '.', ref, alt_str, '.', '.', '.', "GT"]
                # default to ref (0) if alt_dict has no variant for this (pos, replicate_name)
                genotypes = [str(self.alt_dict.get((pos, replicate_name), 0)) for replicate_name in self.replicate_names]
                fields.extend(genotypes)
                f.write('\t'.join(fields))
                f.write("\n")
                
                
if __name__ == '__main__':
    main()    


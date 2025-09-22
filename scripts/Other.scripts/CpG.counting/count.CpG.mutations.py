# -*- coding: utf-8 -*-

from collections import defaultdict

def parse_fasta(filepath):
    sequences = {}
    with open(filepath, 'r') as f:
        current_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(seq_lines).upper()
                current_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_id:
            sequences[current_id] = ''.join(seq_lines).upper()
    return sequences

def extract_cpg_mutations(snps_path, ref_seqs, qry_seqs):
    cpg_mutations = defaultdict(int)
    total_cpg_sites = 0

    with open(snps_path, 'r') as f:
        for line in f:
            if line.startswith('>') or line.startswith('[S1]') or line.startswith('[S2]'):
                continue

            parts = line.strip().split()
            if len(parts) < 12:
                continue

            pos = int(parts[0]) - 1  # Convert to 0-based
            ref_base = parts[1]
            qry_base = parts[2]
            ref_id = parts[10]
            qry_id = parts[11]

            # Skip indels
            if ref_base == '.' or qry_base == '.':
                continue

            ref_seq = ref_seqs.get(ref_id)
            qry_seq = qry_seqs.get(qry_id)
            if not ref_seq or not qry_seq:
                continue

            # Check if the SNP affects a CG site in the reference
            if pos > 0 and ref_seq[pos-1:pos+1] == 'CG':
                total_cpg_sites += 1
                new_dinuc = ref_seq[pos-1] + qry_base
                cpg_mutations[f'CG->{new_dinuc}'] += 1
            elif pos < len(ref_seq) - 1 and ref_seq[pos:pos+2] == 'CG':
                total_cpg_sites += 1
                new_dinuc = qry_base + ref_seq[pos+1]
                cpg_mutations[f'CG->{new_dinuc}'] += 1

    return total_cpg_sites, dict(cpg_mutations)


# === USAGE EXAMPLE ===

# File paths (update these to your actual file paths)
snps_file = "dna.diff.results.from.Mummer/out.snps"
ref_fasta = "B6.allele.specific.DMRs.fa"
qry_fasta = "FVB.allele.specific.DMRs.fa"

# Load sequences
ref_seqs = parse_fasta(ref_fasta)
qry_seqs = parse_fasta(qry_fasta)

# Analyze CpG mutations
total_cpg, cpg_mutation_counts = extract_cpg_mutations(snps_file, ref_seqs, qry_seqs)

# Print results
print(f"Total CpG sites mutated: {total_cpg}")
for mutation, count in sorted(cpg_mutation_counts.items(), key=lambda x: -x[1]):
    print(f"{mutation}: {count}")
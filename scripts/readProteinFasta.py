import csv
import argparse

def read_fasta(file_path):
    """Reads a FASTA file and returns a list of (proteinName, proteinSequence) tuples."""
    sequences = []
    proteinName = None
    sequence_parts = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if proteinName:
                    sequences.append((proteinName, ''.join(sequence_parts)))
                proteinName = line[1:]  # remove '>'
                sequence_parts = []
            else:
                sequence_parts.append(line)

        # add the last record
        if proteinName:
            sequences.append((proteinName, ''.join(sequence_parts)))

    return sequences

def calculate_net_charge(sequence, ph=7.4):
    """Estimates the net charge of a protein sequence at a given pH (default 7.4)."""

    # pKa values for ionizable side chains
    pKa_values = {
        'C_term': 2.2,
        'N_term': 9.6,
        'D': 3.9,  # Aspartic acid
        'E': 4.1,  # Glutamic acid
        'H': 6.0,  # Histidine
        'K': 10.5, # Lysine
        'R': 12.5  # Arginine
    }

    # Count ionizable residues
    counts = {aa: sequence.count(aa) for aa in 'DEHKR'}
    
    # Start with N-terminus and C-terminus
    net_charge = 0
    net_charge += 1 / (1 + 10**(ph - pKa_values['N_term']))   # N-terminus (positive)
    net_charge -= 1 / (1 + 10**(pKa_values['C_term'] - ph))   # C-terminus (negative)

    # Add acidic residues (negative)
    net_charge -= counts['D'] * (1 / (1 + 10**(pKa_values['D'] - ph)))
    net_charge -= counts['E'] * (1 / (1 + 10**(pKa_values['E'] - ph)))

    # Add basic residues (positive)
    net_charge += counts['H'] * (1 / (1 + 10**(ph - pKa_values['H'])))
    net_charge += counts['K'] * (1 / (1 + 10**(ph - pKa_values['K'])))
    net_charge += counts['R'] * (1 / (1 + 10**(ph - pKa_values['R'])))

    return round(net_charge, 2)


def calculate_molecular_weight(sequence):
    """Returns the molecular weight of a protein sequence in Daltons (g/mol)."""
    # Average monoisotopic masses of amino acids (in Daltons)
    aa_weights = {
        'A': 89,  'R': 174, 'N': 132, 'D': 133,
        'C': 121, 'E': 147.04259, 'Q': 146, 'G': 75,
        'H': 155, 'I': 131, 'L': 131, 'K': 146,
        'M': 149, 'F': 165, 'P': 115,  'S': 105,
        'T': 119, 'W': 204, 'Y': 181, 'V': 117.
    }

    weight = 0.0
    for aa in sequence:
        if aa in aa_weights:
            weight += aa_weights[aa]
        else:
            raise ValueError(f"Invalid amino acid: {aa}")

    # Subtract weight of water (18.01528) per peptide bond (n-1) to get actual MW
    # Then add weight of one water for the termini
    if len(sequence) > 0:
        weight -= (len(sequence) - 1) * 18

    return round(weight, 2)

def percent_non_polar(sequence):
    """Calculates the percentage of non-polar amino acids in the sequence."""
    non_polar_residues = set(['A', 'V', 'L', 'I', 'M', 'F', 'W'])
    
    sequence = sequence.upper()
    total = len(sequence)
    if total == 0:
        return 0.0

    non_polar_count = sum(1 for aa in sequence if aa in non_polar_residues)
    percent = (non_polar_count / total) * 100
    return round(percent, 2)


def process_sequence(proteinName, sequence):
    """Placeholder processing function — returns processed info."""
    return {
        'proteinName': proteinName,
        'length': len(sequence),
        'peptideSequence': sequence,
        'molecularWeight': calculate_molecular_weight(sequence),
        'charge': calculate_net_charge(sequence),
        'np%': percent_non_polar(sequence)
    }

def write_to_csv(data, output_path):
    """Writes processed data to a CSV file."""
    with open(output_path, 'w', newline='') as csvfile:
        fieldnames = ['proteinName', 'length', 'peptideSequence', 'molecularWeight', 'charge', 'np%']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in data:
            writer.writerow(row)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a FASTA file and output a CSV.')
    parser.add_argument('input_fasta', help='Path to the input FASTA file')
    parser.add_argument('output_csv', help='Path to the output CSV file')
    args = parser.parse_args()

    sequences = read_fasta(args.input_fasta)
    processed = [process_sequence(h, s) for h, s in sequences]
    write_to_csv(processed, args.output_csv)

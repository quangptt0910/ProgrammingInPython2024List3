import argparse
import json
import pickle
import xml.etree.ElementTree as ET

from Sequence import ProteinSequence
"""
Write a function fasta_to_prot_seq(fasta_file) that opens a UniProt-style FASTA
file and returns a list of ProteinSequence objects.

"""
def fasta_to_prot_seq(fasta_file):
    """
    Parse a FASTA file and convert it to a list of ProteinSequence objects.

    Args:
        fasta_file (str): Path to the FASTA file to be parsed

    Returns:
        list: A list of ProteinSequence objects

    Raises:
        FileNotFoundError: If the input file cannot be found
        ValueError: If the FASTA file is improperly formatted
    """
    identifier = []
    data = []
    data_by_line = ''
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if not line:
                    continue # skip empty line if it has
                if line.startswith('>'):
                    identifier.append(line.strip().replace('>',''))
                    if data_by_line: # Check if there is some line of sequence already be collected
                        data.append(data_by_line)
                        data_by_line = '' # reset the sequence variable
                else:
                    data_by_line += line.strip()

            if data_by_line:
                data.append(data_by_line) # Last sequence data

        list_of_protein_sequence = []
        for i in range(len(identifier)):
            list_of_protein_sequence.append(ProteinSequence(identifier[i], data[i]))

        # Validate input
        if len(identifier) != len(data):
            raise ValueError("Mismatch between identifiers and sequences in FASTA file")

        # Create ProteinSequence objects
        return [ProteinSequence(identifier, sequence)
            for identifier, sequence in zip(identifier, data)]

    except IOError:
        raise FileNotFoundError(f"File not found/ Unable to read file: {fasta_file}")


"""
Write a function fasta_id_to_dict(identifier) that parses (breaks down) the
FASTA file identifier to extract:
• primary access identifier („accession”),
• full name of the protein,
• name of the organism.

https://www.uniprot.org/help/fasta-headers
>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
"""
def fasta_id_to_dict(identifier):
    """
    accession number : UniqueIdentifier
    full name : ProteinName
    name of the organism : OS=OrganismName - scientific name of the organism

    :param identifier:
    :return:
    """
    identifier = identifier.replace('>','')
    parts = identifier.split('|')

    parsed_identifier = {
        'accession': '',
        'protein name': '',
        'organism': ''
    }

    parsed_identifier['accession'] = parts[1]

    details = parts[2] # after the 2nd '|'

    details_section = details.split(' OS=') # split protein name and organism to 2 different str

    entryName_proteinName = details_section[0]
    entry_name = entryName_proteinName.split(' ', 1)[0]
    protein_name = entryName_proteinName.split(' ', 1)[1]
    parsed_identifier['protein name'] = protein_name

    OS_to_end = details_section[1]
    OS_part = OS_to_end.split(' OX=')
    organism_name = OS_part[0]

    parsed_identifier['organism'] = organism_name

    return parsed_identifier


# Argument part
def save_as_pickle(data, filename):
    """Save data to a pickle file."""
    with open(filename, 'wb') as f:
        pickle.dump(data, f)


def save_as_json(data, filename):
    """Save data to a JSON file."""
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)


def save_as_xml(data, filename):
    """Save data to an XML file."""
    root = ET.Element("Proteins")
    for item in data:
        protein = ET.SubElement(root, "Protein")
        for key, val in item.items():
            element = ET.SubElement(protein, key)
            element.text = str(val)
    tree = ET.ElementTree(root)
    tree.write(filename)


def main():
    """Main function to process FASTA file.

    # Basic usage with just the FASTA file
    python fasta.py --fasta uniprot_albumin.fasta

    # Save to all three formats
    python fasta.py --fasta uniprot_albumin.fasta --pickle output.pkl --json output.json --xml output.xml

    # Save to specific formats
    python fasta.py --fasta uniprot_albumin.fasta --pickle output.pkl
    python fasta.py --fasta uniprot_albumin.fasta --json output.json
    python fasta.py --fasta uniprot_albumin.fasta --xml output.xml
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process a UniProt-style FASTA file.")
    parser.add_argument('--fasta', required=True, help="Path to the UniProt-style FASTA file")
    parser.add_argument('--pickle', help="Output pickle file name")
    parser.add_argument('--json', help="Output JSON file name")
    parser.add_argument('--xml', help="Output XML file name")

    # Parse arguments
    args = parser.parse_args()

    # Read protein sequences from FASTA file
    list_seq = fasta_to_prot_seq(args.fasta)

    # Prepare protein data
    protein_data = []
    for seq in list_seq:
        # Extract identifier information
        prot_dict = fasta_id_to_dict(seq.identifier)

        # Add sequence to the dictionary
        prot_dict['sequence'] = seq.data

        protein_data.append(prot_dict)

    # Save data in specified formats
    if args.pickle:
        save_as_pickle(protein_data, args.pickle)
    if args.json:
        save_as_json(protein_data, args.json)
    if args.xml:
        save_as_xml(protein_data, args.xml)

if __name__ == '__main__':
    list_seq = fasta_to_prot_seq('uniprot_albumin.fasta')
    for seq in list_seq:
        print(seq)

    for seq in list_seq:
        print(fasta_id_to_dict(seq.identifier))

    #main()
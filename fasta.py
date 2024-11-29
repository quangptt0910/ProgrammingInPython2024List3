import argparse
import json
import os
import pickle
import sys
import xml.etree.ElementTree as ET
import requests

from Sequence import ProteinSequence

# pip install requests
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

    :param identifier: string of FASTA identifier
    :return: dictionary of parsed identifier
    """
    identifier = identifier.replace('>','')
    parts = identifier.split('|')

    parsed_identifier = {'accession': parts[1],
                         'protein name': '',
                         'organism': ''}

    details = parts[2] # after the 2nd '|'

    details_section = details.split(' OS=') # split protein name and organism to 2 different str

    entryName_proteinName = details_section[0]
    entry_name = entryName_proteinName.split(' ', 1)[0]
    protein_name = entryName_proteinName.split(' ', 1)[1] # Take the proteinName part
    parsed_identifier['protein name'] = protein_name

    ox_to_end = details_section[1]
    ox_part = ox_to_end.split(' OX=') # take the OX organism part
    organism_name = ox_part[0]

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


"""
e) Modify the program written in (d), so it takes as an argument not a file, but the name
of the protein (--protein) and the maximum number of returned sequences --max_seq.
Use the intuitive requests package, known as HTTP for Humans, to access UniProt’s API
(Application Programming Interface) to automatically download data.
"""
# https://www.uniprot.org/help/api_queries
def fetch_uniprot_sequences(protein_name, max_sequences):
    """
    Fetch protein sequences from UniProt API based on protein name.

    Args:
        protein_name (str): Name of the protein to search
        max_sequences (int, optional): Maximum number of sequences to retrieve. Defaults to 10.

    Returns:
        list: A list of ProteinSequence objects

    Raises:
        ValueError: If no sequences are found or API request fails
    """
    # Base URL for UniProt API
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    # Parameters for the API query
    params = {
        'query': protein_name,
        'format': 'fasta',
        'size': max_sequences
    }

    try:
        # Send GET request to UniProt API
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise an exception for bad responses

        # Check if response is empty
        if not response.text.strip():
            raise ValueError(f"No sequences found for protein: {protein_name}")

        # Save the FASTA response to a temporary file
        with open('temp_uniprot_search.fasta', 'w') as f:
            f.write(response.text)

        # Use existing fasta_to_prot_seq function to parse the downloaded FASTA
        return fasta_to_prot_seq('temp_uniprot_search.fasta')

    except requests.RequestException as e:
        raise ValueError(f"API request failed: {e}")


def fetch_main():
    """
    Fetch and process protein sequences from UniProt using either a FASTA file or direct protein search.

    This function provides flexible functionality for retrieving protein sequence data through two primary methods:
    1. Reading sequences from a local UniProt-style FASTA file
    2. Fetching sequences directly from the UniProt API by protein name

    Parameters:
    None (uses command-line arguments)

    Optional Command-Line Arguments:
    --fasta : str, optional
        Path to a UniProt-style FASTA file containing protein sequences.

    --protein : str, optional
        Name of the protein to fetch sequences for from UniProt.

    --max_seq : int, optional
        Maximum number of sequences to retrieve when searching by protein name.
        Default is None, which may return all available sequences.

    --pickle : str, optional
        Output file path for saving protein data in pickle format.

    --json : str, optional
        Output file path for saving protein data in JSON format.

    --xml : str, optional
        Output file path for saving protein data in XML format.

    Usage Examples:
    # Read sequences from a FASTA file
    python fasta.py --fasta uniprot_albumin.fasta --json output.json

    # Fetch up to 5 sequences for albumin
    python fasta.py --protein albumin --max_seq 5 --json albumin.json

    # Fetch and save in multiple formats
    python fasta.py --protein albumin --max_seq 10 --pickle albumin.pkl --json albumin.json --xml albumin.xml

    Raises:
    FileNotFoundError: If the specified FASTA file does not exist.
    Exception: For any errors encountered during UniProt sequence retrieval.

    Note:
    - Requires external libraries: argparse, requests
    - Depends on helper functions: fasta_to_prot_seq(), fetch_uniprot_sequences(),
      fasta_id_to_dict(), save_as_pickle(), save_as_json(), save_as_xml()
    """

    parser = argparse.ArgumentParser(description="Fetch protein sequences from UniProt API.")
    parser.add_argument('--fasta', help="Path to the UniProt-style FASTA file")
    parser.add_argument('--protein', help="Name of the protein")
    parser.add_argument('--max_seq',type=int, help="Maximum number of sequences to return")
    parser.add_argument('--pickle', help="Output pickle file name")
    parser.add_argument('--json', help="Output JSON file name")
    parser.add_argument('--xml', help="Output XML file name")

    # Parse arguments
    args = parser.parse_args()

    # Validate input
    if not args.fasta and not args.protein:
        raise ValueError("Either a FASTA file or a protein name must be provided.")

    # Validate max_seq
    if args.max_seq is not None and (args.max_seq <= 0 or not isinstance(args.max_seq, int)):
        raise ValueError("Maximum sequences must be a positive integer.")

    protein_data = []
    try:
        if args.fasta:
            if not os.path.exists(args.fasta):
                raise FileNotFoundError(f"The specified FASTA file does not exist: {args.fasta}")

            # Read protein sequences from FASTA file
            list_seq = fasta_to_prot_seq(args.fasta)

            # Prepare protein data
            for seq in list_seq:
                # Extract identifier information
                prot_dict = fasta_id_to_dict(seq.identifier)

                # Add sequence to the dictionary
                prot_dict['sequence'] = seq.data

                protein_data.append(prot_dict)

        # Fetch from UniProt
        elif args.protein:
            try:
                # Fetch protein sequences from UniProt
                list_seq = fetch_uniprot_sequences(args.protein, args.max_seq)
                print(f"Retrieved {len(list_seq)} sequences for protein: {args.protein}")
            except Exception as e:
                print(f"Error fetching UniProt sequences: {e}", file=sys.stderr)
                raise RuntimeError(f"Failed to retrieve sequences for protein: {args.protein}")

        # Check if any sequences were retrieved
        if not list_seq:
            raise RuntimeError("No protein sequences were found.")
        # Prepare protein data
        for seq in list_seq:
            try:
                # Extract identifier information
                prot_dict = fasta_id_to_dict(seq.identifier)

                # Add sequence to the dictionary
                prot_dict['sequence'] = seq.data
                protein_data.append(prot_dict)
            except Exception as e:
                print(f"Could not process sequence {seq.identifier}: {e}", file=sys.stderr)

            # Save data in specified formats with error handling
        output_formats = [
            ('pickle', args.pickle, save_as_pickle),
            ('json', args.json, save_as_json),
            ('xml', args.xml, save_as_xml)
        ]

        for format_name, output_file, save_func in output_formats:
            if output_file:
                try:
                    save_func(protein_data, output_file)
                    print(f"Successfully saved {format_name} to {output_file}")
                except Exception as e:
                    print(f"Error saving {format_name} file: {e}", file=sys.stderr)
                    raise IOError(f"Failed to save {format_name} file: {output_file}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        raise

if __name__ == "__main__":
    fetch_main()

from abc import ABC, abstractmethod


class Sequence(ABC):
    """
    An abstract base class representing a biological sequence.

    This class provides the foundation for different types of biological sequences
    (DNA, RNA, protein) with common functionality for sequence manipulation
    and analysis.

    Attributes:
        identifier (str): A unique identifier for the sequence.
        data (str): The actual sequence data.
        valid_chars (set): A set of valid characters allowed in the sequence.

    Methods:
        __len__(): Returns the length of the sequence.
        __str__(): Returns a string representation of the sequence.
        mutate(position, value): Modifies the sequence at a specific position.
        find_motif(motif): Finds all occurrences of a motif in the sequence.
        complement(): Returns the complementary sequence.

    Raises:
        ValueError: If the identifier is empty or if invalid characters are found in the sequence.
        IndexError: If the mutation position is out of range.
    """
    def __init__(self, identifier, data):
        """Initialize a new Sequence instance.

            Args:
                identifier (str): A unique identifier for the sequence.
                data (str): The sequence data.
            Raises:
                ValueError: If identifier is empty or if data contains invalid characters.
        """
        if not identifier:
            raise ValueError("Sequence identifier should not be empty")
        self._identifier = identifier

        # invalid_chars = [c for c in data.upper() if c not in self.valid_chars]
        if not all(c in self.valid_chars for c in data.upper()):
            raise ValueError(f"Invalid characters in sequence data!")
        self._data = data.upper()

    @property
    def identifier(self):
        """
            Get identifier of sequence.
            Returns: (str) identifier of sequence.
        """
        return self._identifier

    @identifier.setter
    def identifier(self, value):
        """
            Set the sequence identifier.
            Args:
                value (str): The new identifier value.
            Raises:
                ValueError: If the identifier is empty.
        """
        if not value:
            raise ValueError('Sequence identifier should not be empty.')
        self._identifier = value

    @property
    def data(self):
        """
            Get data of sequence.
            Returns: (str) data of sequence.
        """
        return self._data

    @data.setter
    def data(self, value):
        """
            Set the sequence data.
            Args:
                value (str): The new sequence data.

            Raises:
                ValueError: If the sequence contains invalid characters.
        """
        invalid_chars = [c for c in value.upper() if c not in self.valid_chars]
        if invalid_chars:
            raise ValueError(f"Invalid characters in sequence: {', '.join(invalid_chars)}")
        self._data = value.upper()

    @property
    @abstractmethod
    def valid_chars(self):
        """
            Get the set of valid characters for this sequence type.
            Returns:
                set: Set of valid characters allowed in the sequence.
        """
        pass

    def __len__(self):
        return len(self.data)

    def __str__(self):
        return '>' + str(self._identifier) + '\n' + str(self._data)

    def mutate(self, position, value):
        """Modify sequence data at given position to given value.
        Args:
            position (int): Zero-based index where modification should occur.
            value (str): New character to insert at the specified position.

        Raises:
            IndexError: If position is out of range.
            ValueError: If the new value contains invalid characters.
        """
        if not (0 <= position < len(self._data)):
            raise IndexError('Position out of range.')
        if value.upper() not in self.valid_chars:
            raise ValueError('Invalid characters in sequence.')

        data_list = list(self._data)
        data_list[position] = value.upper()
        self._data = "".join(data_list)
        # self.data[position] = value WONT work since Python string are immutable

    def find_motif(self, motif):
        """
        Find all occurrences of a motif in the sequence.
            Args:
                motif (str): The sequence motif to search for.

            Returns:
                list: Zero-based positions where the motif occurs in the sequence.
        """
        positions = []
        motif = motif.upper()
        for i in range(len(self.data) - len(motif) + 1):
            if self.data[i:i+len(motif)] == motif:
                positions.append(i)
        return positions


class DNASequence(Sequence):
    """A class representing DNA sequences.

    This class implements the Sequence interface specifically for DNA molecules,
    supporting standard DNA operations including complementation and transcription.

    Attributes:
        identifier (str): A unique identifier for the DNA sequence.
        data (str): The DNA sequence using A, C, G, T nucleotides.
        valid_chars (set): The set {'A', 'C', 'G', 'T'}.

    Methods:
        complement(): Returns the complementary DNA sequence following base pairing rules
            (A↔T, C↔G).
        transcribe(): Converts the DNA sequence to its corresponding RNA sequence
            by replacing T with U.
        find_motif(motif): Finds all occurrences of a specific DNA motif.

    Example:
        >>> dna = DNASequence("gene1", "ATCG")
        >>> str(dna.complement())
        'TAGC'
        >>> rna = dna.transcribe()
        >>> str(rna)
        'AUCG'
    """
    _valid_chars = {'A', 'C', 'G', 'T'}

    @property
    def valid_chars(self):
        return self._valid_chars

    def complement(self):
        """
        Generate the complementary DNA sequence.
            Returns:
                str: The complementary sequence following base pairing rules (A↔T, C↔G).
            Raises:
                ValueError: If the sequence contains invalid characters.
            Example:
                >>> dna = DNASequence("example", "ATCG")
                >>> dna.complement()
                'TAGC'
        """
        complement_pair_up = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        try:
            return ''.join(complement_pair_up[base] for base in self._data)
        except KeyError as e:
            raise ValueError(f"Invalid character in DNA strand: {e}")

    def transcribe(self) -> 'RNASequence':
        """
        Transcribe the DNA sequence into RNA.
            Returns:
                RNASequence: A new RNASequence instance where all 'T' are replaced with 'U'.
            Example:
                >>> dna = DNASequence("example", "ATCG")
                >>> rna = dna.transcribe()
                >>> str(rna)
                'AUCG'
        """
        rna_data = self._data.replace('T', 'U')
        return RNASequence(f"{self.identifier}_RNA", rna_data)


class RNASequence(Sequence):
    """A class representing RNA sequences.

        This class implements the Sequence interface specifically for RNA molecules,
        supporting RNA-specific operations including complementation and translation
        to protein sequences.

        Attributes:
            identifier (str): A unique identifier for the RNA sequence.
            data (str): The RNA sequence using A, C, G, U nucleotides.
            valid_chars (set): The set {'A', 'C', 'G', 'U'}.

        Methods:
            complement(): Returns the complementary RNA sequence following base pairing rules
                (A↔U, C↔G).
            translate(): Converts the RNA sequence to its corresponding protein sequence
                using the standard genetic code.

        Example:
            >>> rna = RNASequence("transcript1", "AUGCUAUGA")
            >>> protein = rna.translate()
            >>> str(protein)
            'ML'
    """

    _valid_chars = {'A', 'C', 'G', 'U'}

    @property
    def valid_chars(self):
        return self._valid_chars

    def complement(self):
        complement_pair_up = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        try:
            return ''.join(complement_pair_up[base] for base in self._data)
        except KeyError as e:
            raise ValueError(f"Invalid character in RNA strand: {e}")

    def translate(self) -> 'ProteinSequence':
        """Translate the RNA sequence into a protein sequence.

        The translation starts from the beginning of the sequence and continues
        until either a stop codon is encountered or the sequence ends. Codons
        are read in groups of three nucleotides.

        Returns:
            ProteinSequence: A new ProteinSequence instance containing the translated amino acid sequence.

        Notes:
            - Translation stops at the first stop codon ('UAA', 'UAG', or 'UGA')
            - Incomplete codons at the end of the sequence are ignored
            - Unknown codons are represented as 'X' in the protein sequence

        Example:
            >>> rna = RNASequence("example", "AUGGCCUAA")
            >>> protein = rna.translate()
            >>> str(protein)
            'MA'
        """
        genetic_code = {
            'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

        protein = []
        for i in range(0, len(self._data) - 2, 3):
            codon = self._data[i:i + 3]
            if len(codon) == 3:  # ensure we have a complete codon
                amino_acid = genetic_code.get(codon, 'X')  # 'X' for unknown codons
                if amino_acid == '*':  # Stop codon
                    break
                protein.append(amino_acid)

        return ProteinSequence(f"{self.identifier}_protein", "".join(protein))


class ProteinSequence(Sequence):
    """A class representing protein sequences.

       This class implements the Sequence interface for protein sequences,
       supporting standard amino acid sequence operations.

       Attributes:
           identifier (str): A unique identifier for the protein sequence.
           data (str): The protein sequence using single-letter amino acid codes.
           valid_chars (set): The set of valid amino acid codes (A-Z) plus stop codon (*).

       Methods:
           find_motif(motif): Finds all occurrences of a specific amino acid motif.

       Notes:
           - Uses standard one-letter amino acid codes (A, C, D, E, F, G, H, I, K, L, M, N,
             P, Q, R, S, T, V, W, Y) plus '*' for stop codons.
           - Sequence is case-insensitive but stored in upper case.

       Example:
           >>> protein = ProteinSequence("protein1", "MGKL")
           >>> protein.find_motif("GK")
           [1]
       """

    # Class attribute: valid amino acid single letter codes plus stop codon
    _valid_chars = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                    'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'}

    @property
    def valid_chars(self):
        """Set of valid characters for protein sequence"""
        return self._valid_chars


if __name__ == "__main__":
    dnaSeq = DNASequence("DNA1", "ATGGCCTAA")
    print(dnaSeq)
    rnaSeq = RNASequence("RNA1","AUGGCCUAA") # "GAUGGAACUUGACUACGUAAAUU")
    print(rnaSeq)
    print(rnaSeq.find_motif("GAA"))
    print(rnaSeq.translate())
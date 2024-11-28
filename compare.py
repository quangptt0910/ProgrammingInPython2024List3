import argparse
import sys
import string
from collections import Counter

from shingle import shingles

# Write program compare.py that compares two text files, provided as input arguments
# --query and --target, by calculating and printing the Jaccard similarity between the
# texts represented by multisets (a.k.a. bag) of their k-shingles. The program should also
# have argument -k to define the size of shingle and switch --remove_punctuation to
# enable or disable the removal of punctuation marks.


#https://stackoverflow.com/questions/34293875/how-to-remove-punctuation-marks-from-a-string-in-python-3-x-using-translate
def remove_punctuation(text):
    """
    Remove punctuation marks from the input text.

    Args:
        text (str): Input text to remove punctuation from.

    Returns:
        str: Text with all punctuation marks removed.
    """
    # Remove all punctuation using string.punctuation
    return text.translate(str.maketrans('', '', string.punctuation))

# https://intro2ml.pages.doc.ic.ac.uk/autumn2021/modules/lab-java/files
def get_shingle_from_file(file_name, k, remove_punct):
    """
    Read text from file and generate k-shingles.

    Args:
        file_name (str): Path to the input file
        k (int): Size of shingles
        remove_punct (bool): Whether to remove punctuation

    Returns:
        list: List of k-shingles
    """
    try:
        with open(file_name, "r") as infile:
            content = infile.read()
    except IOError as e:
        print(f"Error reading file {file_name}: {e}", file=sys.stderr)
        sys.exit(1)

    # Preprocess text
    # https://stackoverflow.com/questions/34293875/how-to-remove-punctuation-marks-from-a-string-in-python-3-x-using-translate
    if remove_punct:
        content = remove_punctuation(content)

    # Generate and return k-shingles
    tokens = content.split()
    return [tuple(tokens[i:i + k]) for i in range(len(tokens) - k + 1)]



# https://en.wikipedia.org/wiki/Jaccard_index
# Jaccard similarity
# https://www.geeksforgeeks.org/how-to-calculate-jaccard-similarity-in-python/
def calc_jaccard(query, target):
    """
    function calculate the Jaccard similarity between two lists
    :param query: first shingle list from shingle()
    :param target: second shingle list from shingle()
    :return: the Jaccard similarity between the two lists
    """
    query_counter = Counter(query) # Counter() for query shingle -> Counter({key: value})
    target_counter = Counter(target) # Counter() for target shingle -> Counter({key: value})
    intersect_counter = {} # dict to store key and value of intersection between query_counter and target_counter

    # find the intersection between two shingles
    for key in query_counter.keys():
        if key in target_counter:
            if key not in intersect_counter:
                intersect_counter[key] = min(query_counter[key], target_counter[key])

    print(intersect_counter)

    intersect = sum(intersect_counter.values())
    union = sum(query_counter.values()) + sum(target_counter.values()) - intersect
    # formula : intersect / union
    return f'{intersect}/{union}'



def get_jaccard_similarity():
    """
        Calculate Jaccard similarity between two text files using k-shingles.

        This function compares two text files by breaking them down into k-shingles
        (subsequences of k characters) and calculating their Jaccard similarity.

        Usage in terminal:
        python script_name.py --query path/to/query_file.txt --target path/to/target_file.txt -k 3 [optional_flags]

        Arguments:
        --query (str, required): Path to the first text file (query file)
        --target (str, required): Path to the second text file (target file)
        -k (int, required): Size of the shingles (character subsequence length)
        --remove_punctuation (optional): If present, removes punctuation before creating shingles

        Example usages:
        1. Basic usage:
           python compare.py --query file1.txt --target file2.txt -k 4

        2. With punctuation removal:
           python compare.py --query file1.txt --target file2.txt -k 4 --remove_punctuation

        Returns:
        String (of fraction): Jaccard similarity between the two text files' shingle sets

        Raises:
        argparse.ArgumentTypeError: If required arguments are missing or invalid
        """
    parser = argparse.ArgumentParser(description='Function that compares two text files, '
                                                 'by calculating and printing the Jaccard similarity between the '
                                                 'texts represented by multisets (a.k.a. bag) of their k-shingles')

    parser.add_argument('--query', required=True, help='First text file (as query)')
    parser.add_argument('--target', required=True, help='Second text file (as target)')
    parser.add_argument('-k', required=True, type=int, help='size of the shingle')
    parser.add_argument('--remove_punctuation', action='store_true', help='a switch to enable or disable the removal of punctuation marks')

    args = parser.parse_args()
    # get the shingles from the text files
    query = get_shingle_from_file(args.query, args.k, args.remove_punctuation)
    target = get_shingle_from_file(args.target, args.k, args.remove_punctuation)

    # calculate and print the Jaccard similarity
    print(calc_jaccard(query, target))

if __name__ == "__main__":
    get_jaccard_similarity()



import argparse
import sys
from collections import Counter

"""
Function to calculate n most common k-shingles
Example terminal/console input:
$ python shingle.py -n 3 -k 2
Let it be, let it be, let it be
whisper words of wisdom, let it be
>>> Result:
let it: 3
it be,: 2
be, let: 2
"""
# function shingles(t, k)
def shingles(t, k):
    """
        Generate k-shingles from a list of tokens.

        Args:
            t (list): List of tokens to generate shingles from
            k (int): Size of each shingle

        Returns:
            list: List of k-shingles as tuples
        """
    # Handle edge cases
    if k <= 0 or k > len(t):
        return []

    list_of_kShingle = []
    if k == 0: return []
    for i in range(0, len(t) - k + 1):
        list_of_kShingle.append(tuple(t[i:i + k])) # tuple() will work
    return list_of_kShingle

# https://realpython.com/python-counter/
# >>> word = "mississippi"
# >>> counter = {}
#
# >>> for letter in word:
# ...     counter[letter] = counter.get(letter, 0) + 1
#
# >>> counter
# {'m': 1, 'i': 4, 's': 4, 'p': 2}
def most_common_k_shingles(t, k, n):
    """
        Find the n most common k-shingles in the given tokens.

        Args:
            t (list): List of tokens
            k (int): Size of shingles
            n (int): Number of most common shingles to return

        Returns:
            list: List of (shingle, count) tuples
        """
    # Generate k-shingles and count their occurrences
    shingle_list = shingles(t, k)
    counter = Counter(shingle_list)

    # Return n most common shingles
    return counter.most_common(n)

# function to input and print the most common k shingles with arguments
def get_shingle():
    parser = argparse.ArgumentParser(description='Function to calculate n most common k-shingles')
    parser.add_argument('-n', required=True, type=int, help='n')
    parser.add_argument('-k', required=True, type=int, help='k')

    args = parser.parse_args()

    # Read all input
    input_tokens = []
    for line in sys.stdin:
        input_tokens.extend(line.split())

    # Compute and print most common k-shingles
    try:
        result = most_common_k_shingles(input_tokens, args.k, args.n)

        # print the result
        print("Most common k-shingles:")
        for shingle, count in result:
            print(f"{shingle}: {count}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    get_shingle()
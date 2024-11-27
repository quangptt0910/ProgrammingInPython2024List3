import argparse
import sys
def get_student_info():
    parser = argparse.ArgumentParser(description='Collect student information.')
    parser.add_argument('--name', required=True, help='Student name')
    parser.add_argument('--age', type=int, required=True, help='Student age')
    parser.add_argument('--major', required=True, help='Student major')

    args = parser.parse_args()

    print(f"\nStudent Information:")
    print(f"Name: {args.name}")
    print(f"Age: {args.age}")
    print(f"Major: {args.major}")

    # Read additional comments or notes interactively using sys.stdin.read()
    print("\nEnter additional comments or notes (Ctrl-D/Z to end):")
    try:
        notes = sys.stdin.read()
        if notes.strip():  # Check if there's any non-whitespace content
            print(f"\nAdditional Comments/Notes:")
            print(notes)
    except EOFError:
        pass

if __name__ == "__main__":
    get_student_info()

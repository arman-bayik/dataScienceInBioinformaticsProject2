import numpy
import sys
import random


def main():
    filename = "SeaCucumberGlobin.txt"
    create_hp_model(filename)


# method to generate the initial population
def create_initial_population(size, aa_sequence_source):
    # # Create the population of directions; F,L,and R.
    # letters = "FLR"
    # # stores the population of directions
    # directions = []
    # # initialize the population to 5 children.
    # for j in range(100):
    #     direction = ''.join(random.choice(letters) for i in range(size))
    #     directions.append(direction)

    initial_population = []  # variable to store the initial population
    while len(initial_population) < size:
        member = []  # variable to represent a member of the initial population in the form [traverse_sequence, d2_lattice, score]
        d2_lattice = numpy.ndarray([2*size+2, 2*size+2], dtype=object)
        for r in range(0, len(d2_lattice), 1):
            for c in range(0, len(d2_lattice[r]), 1):
                # format is: [index, amino acid, previous index, next index]
                d2_lattice[r, c] = ["None", "None", "None", "None"]
        r_index = int(size)
        c_index = int(size)  # indices in the d2_lattice matrix
        orientation = 'E'  # current orientation in units of cardinal direction
        traverse_sequence = []  # sequence of traversal moves
        bump = False  # flag for if there's a bump
        for i in range(0, len(aa_sequence_source), 1):
            bump = True
            # insert amino acids into the positions they should go
            if i == 0:
                d2_lattice[r_index, c_index] = [i, aa_sequence_source[i], "None", i+1]
            elif i == len(aa_sequence_source)-1:
                d2_lattice[r_index, c_index] = [i, aa_sequence_source[i], i-1, "None"]
            else:
                d2_lattice[r_index, c_index] = [i, aa_sequence_source[i], i-1, i+1]

            # loop to catch a bump and try to avoid it
            tries_remaining = 27  # has 27 tries to choose a traversal that doesn't result in a bump
            while bump and tries_remaining != 0:
                # randomly determine which direction to traverse
                current_traversal = 'X'
                traverse_type = random.randint(0, 2)
                if traverse_type == 0:
                    current_traversal = 'L'
                elif traverse_type == 1:
                    current_traversal = 'R'
                elif traverse_type == 2:
                    current_traversal = 'F'

                # determine index of new insertion
                current_traversal_info = get_traversal_info(current_traversal, orientation)
                new_r_index = r_index + current_traversal_info[0]
                new_c_index = c_index + current_traversal_info[1]

                # determines if the current traversal will result in a bump based on contents of new index
                # if not, add current traversal to traverse_sequence, update orientation, and update indices
                if d2_lattice[new_r_index, new_c_index][0] == "None":
                    traverse_sequence.append(current_traversal)
                    orientation = current_traversal_info[2]
                    r_index = new_r_index
                    c_index = new_c_index
                    bump = False

                # decrements the number of tries remaining
                tries_remaining = tries_remaining - 1

            # if the sequence contains a bump that could not be corrected, disregard the rest of the sequence
            if bump:
                break

        # if the sequence doesn't have a bump, add it to the initial population
        if not bump:
            member.append(traverse_sequence)
            member.append(d2_lattice)
            member.append(score_lattice(d2_lattice))
            initial_population.append(member)
            print(len(initial_population))

    return initial_population


# method that takes in a traversal and an orientation and returns [delta_r, delta_c, new_orientation]
def get_traversal_info(traversal, orientation):
    delta_r = 9
    delta_c = 9
    new_orientation = 'X'
    if traversal == 'L':
        if orientation == 'N':
            new_orientation = 'W'
            delta_r = 0
            delta_c = -1
        elif orientation == 'S':
            new_orientation = 'E'
            delta_r = 0
            delta_c = 1
        elif orientation == 'E':
            new_orientation = 'N'
            delta_r = -1
            delta_c = 0
        elif orientation == 'W':
            new_orientation = 'S'
            delta_r = 1
            delta_c = 0
        else:
            new_orientation = 'X'

    elif traversal == 'R':
        if orientation == 'N':
            new_orientation = 'E'
            delta_r = 0
            delta_c = 1
        elif orientation == 'S':
            new_orientation = 'W'
            delta_r = 0
            delta_c = -1
        elif orientation == 'E':
            new_orientation = 'S'
            delta_r = 1
            delta_c = 0
        elif orientation == 'W':
            new_orientation = 'N'
            delta_r = -1
            delta_c = 0
        else:
            new_orientation = 'X'

    elif traversal == 'F':
        new_orientation = orientation
        if orientation == 'N':
            delta_r = -1
            delta_c = 0
        elif orientation == 'S':
            delta_r = 1
            delta_c = 0
        elif orientation == 'E':
            delta_r = 0
            delta_c = 1
        elif orientation == 'W':
            delta_r = 0
            delta_c = -1
        else:
            new_orientation = 'X'

    return [delta_r, delta_c, new_orientation]


# def get_cardinal_direction(path):
#     count = 0
#
#     for letter in path:
#         if letter == "L":
#             count = count - 1
#         if letter == "R":
#             count = count + 1
#         if letter
#


def score_lattice(d2_lattice):
    # IMPLEMENT LATER
    return -sys.maxsize


def create_hp_model(filename):
    # H-Hydrophobic, F-Hydrophilic
    amino_dict = {'X': 'H', 'A': 'H', 'R': 'F', 'N': 'F', 'D': 'F',
                  'C': 'H', 'Q': 'F', 'E': 'F', 'G': 'F', 'H': 'F',
                  'I': 'H', 'L': 'H', 'K': 'F', 'M': 'H', 'F': 'H',
                  'P': 'H', 'S': 'F', 'T': 'F', 'W': 'H', 'Y': 'H', 'V': 'H'}

    # open and read amino acid sequence from the fasta file
    try:
        input_file = open(filename, "r")
    except FileNotFoundError:
        print("Error: specified input file, '" + filename + "', does not exist.")
        sys.exit(1)

    aa_sequence_raw = []
    for line in input_file:
        if '>' not in line:
            aa_sequence_raw.append(line.rstrip('\r\n'))
    aa_sequence = []
    for sequence_segment in aa_sequence_raw:
        for aa in sequence_segment:
            aa_sequence.append(aa)

    init_pop = create_initial_population(100, aa_sequence)

    # # H-Hydrophobic, F-Hydrophilic
    # amino_dict = {"X": "H", "A": "H", "R": "F", "N": "F", "D": "F",
    #               "C": "H", "Q": "F", "E": "F", "G": "F", "H": "F",
    #               "I": "H", "L": "H", "K": "F", "M": "H", "F": "H",
    #               "P": "H", "S": "F", "T": "F", "W": "H", "Y": "H", "V": "H"}
    # # This shows the proteins that are hydrophobic and hydrophilic in the input sequence
    # h_f_list = []

    # open and read from the fasta file
    # try:
    #     input_file = open(filename, "r")
    # except FileNotFoundError:
    #     print("Error: specified input file, '" + filename + "', does not exist.")
    #     sys.exit(1)
    #
    # raw_sequences = []
    # amino_acid_sequences = []
    # single_aminos = []
    # # Trim the header from the input file
    # for raw_sequence in input_file:
    #     if raw_sequence[0] != ">":
    #         raw_sequences.append(raw_sequence)
    # # Trim the new line character from the
    # for sequence in raw_sequences:
    #     amino_acid_sequence = sequence[0:-1]
    #     amino_acid_sequences.append(amino_acid_sequence)
    # for polypeptide in amino_acid_sequences:
    #     for i in range(len(polypeptide)):
    #         single_aminos.append(polypeptide[i])
    #
    # sequence_length = len(single_aminos)
    # print(single_aminos)
    # print(sequence_length)
    #
    # # fill up the h_f_list for the input file.
    # for x in single_aminos:
    #     amino_property = amino_dict.get(x)
    #     h_f_list.append(amino_property)
    # print(h_f_list)
    #
    # # |--------------------CREATING THE POPULATION--------------------|
    # directions = create_population(sequence_length)
    #
    # # |--------------------SELECTION--------------------|
    # # 1. Evaluating fitness
    # # If the direction contains LLL or FFF then we can't consider it fit
    # has_kids = False
    # while not has_kids:
    #     # This contains the fit children
    #     fit_directions = []
    #     for direction in directions:
    #         if not direction.__contains__("LLLL") | direction.__contains__("RRRR"):
    #             fit_directions.append(direction)
    #     # This happens if all children have bumps. In this case try making a new population until you have kids that
    #     # don't have bumps
    #     if len(fit_directions) == 0:
    #         has_kids = False
    #         directions = create_population(sequence_length)
    #     else:
    #         has_kids = True
    #
    # # At this point fit_directions contains the directions that don't cause bumps.
    # # Use the directions to plot the Hydrophobic or Hydrophillic properties in a 2-D array.
    # # First get the cardinal direction based off of each direction.
    # count = 0
    # d = "W"
    # facing = []
    # fd = []
    # for direction in fit_directions:
    #     for letter in direction:
    #         # remove fd later
    #         fd.append(letter)
    #
    #         if letter == "L":
    #             count = count - 1
    #         elif letter == "R":
    #             count = count + 1
    #         else:
    #             count = count + 0
    #
    #         if count % 4 == 0:
    #             d = "E"
    #         elif count % 4 == 1:
    #             d = "S"
    #         elif count % 4 == 2:
    #             d = "W"
    #         elif count % 4 == 3:
    #             d = "N"
    #         facing.append(d)
    #
    #     # At this point we have the correct cardinal directions based off of each direction we had to go in.
    #     # Fill out a 2D matrix that is sequence_length by sequence_length; starting from the middle;
    #     # in the lattice represent H(Hydrophobic) as 1 and F(Hydrophyllic) as 0
    #     x = int(sequence_length / 2)
    #     y = int(sequence_length / 2)
    #     lattice = [[2] * sequence_length] * sequence_length
    #     matrix = numpy.array(lattice)
    #
    #     # If the sequence starts with a Hydrophobic protein insert a "1" at the center of the matrix insert "0"
    #     # otherwise
    #     if h_f_list[0] == "H":
    #         matrix[x][y] = 1
    #     else:
    #         matrix[x][y] = 0
    #     for i in range(sequence_length):
    #         if facing[i - 1] == "N":
    #             if fd[i] == "L":
    #                 y = y - 1
    #             elif fd[i] == "R":
    #                 y = y + 1
    #             elif fd[i] == "F":
    #                 x = x + 1
    #
    #         elif facing[i - 1] == "E":
    #             if fd[i] == "L":
    #                 x = x + 1
    #             elif fd[i] == "R":
    #                 x = x - 1
    #             elif fd[i] == "F":
    #                 y = y + 1
    #
    #         elif facing[i - 1] == "S":
    #             if fd[i] == "L":
    #                 y = y + 1
    #             elif fd[i] == "R":
    #                 y = y - 1
    #             elif fd[i] == "F":
    #                 x = x - 1
    #
    #         elif facing[i - 1] == "W":
    #             if fd[i] == "L":
    #                 x = x - 1
    #             elif fd[i] == "R":
    #                 x = x + 1
    #             elif fd[i] == "F":
    #                 y = y - 1
    #
    #         if h_f_list[i] == "H":
    #             matrix[x][y] = 1
    #         else:
    #             matrix[x][y] = 0
    #
    #
    #
    #
    #
    #
    #
    #     print(fd)
    #     print(matrix)
    # print(facing)
    return 0


if __name__ == '__main__':
    main()

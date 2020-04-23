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


def score_lattice(d2_lattice):
    # dictionary of all amino acids
    # H-Hydrophobic, F-Hydrophilic
    amino_dict = {'X': 'H', 'A': 'H', 'R': 'F', 'N': 'F', 'D': 'F',
                  'C': 'H', 'Q': 'F', 'E': 'F', 'G': 'F', 'H': 'F',
                  'I': 'H', 'L': 'H', 'K': 'F', 'M': 'H', 'F': 'H',
                  'P': 'H', 'S': 'F', 'T': 'F', 'W': 'H', 'Y': 'H', 'V': 'H'}
    score = 0
    considered_aas = []
    # traverse d2_lattice matrix
    for r in range(0, len(d2_lattice), 1):
        for c in range(0, len(d2_lattice), 1):
            is_current_aa_considered = False
            current_aa = d2_lattice[r, c]
            # if the current index contains an amino acid, check for hydrophobia
            if not current_aa[0] == "None":
                if amino_dict[current_aa[1]] == 'H':
                    # check if the current amino acid has already been considered for an H-H bond
                    for i in range(0, len(considered_aas), 1):
                        # if it has been, set flag
                        if considered_aas[i] == current_aa[0]:
                            is_current_aa_considered = True
                            break
                    # if the current amino acid has not already been considered for an H-H bond, check for nearby H's
                    if not is_current_aa_considered:
                        is_above_aa_considered = False
                        is_left_aa_considered = False
                        is_right_aa_considered = False
                        is_below_aa_considered = False
                        above_aa = d2_lattice[r-1, c]
                        left_aa = d2_lattice[r, c-1]
                        right_aa = d2_lattice[r, c+1]
                        below_aa = d2_lattice[r+1, c]
                        # check if nearby amino acids have already been considered for an H-H bond
                        for i in range(0, len(considered_aas), 1):
                            # if above amino acid has been considered, set flag
                            if considered_aas[i] == above_aa[0]:
                                is_above_aa_considered = True
                            # if left amino acid has been considered, set flag
                            if considered_aas[i] == left_aa[0]:
                                is_left_aa_considered = True
                            # if right amino acid has been considered, set flag
                            if considered_aas[i] == right_aa[0]:
                                is_right_aa_considered = True
                            # if below amino acid has been considered, set flag
                            if considered_aas[i] == below_aa[0]:
                                is_below_aa_considered = True
                        # if right amino acid hasn't been considered, check if it is an H
                        if not is_right_aa_considered and not right_aa[1] == "None":
                            if amino_dict[right_aa[1]] == 'H':
                                # if right amino acid is an H, check if it's linked to the current amino acid
                                # if it is an unconsidered H that isn't linked to the current amino acid, it is an H-H connection
                                if not (current_aa[2] == right_aa[0] or current_aa[3] == right_aa[0]):
                                    score = score + 1  # increment score
                                    # add both amino acids to the list of considered amino acids
                                    considered_aas.append(current_aa[0])
                                    considered_aas.append(right_aa[0])
                        # if left amino acid hasn't been considered, check if it is an H
                        if not is_left_aa_considered and not left_aa[1] == "None":
                            if amino_dict[left_aa[1]] == 'H':
                                # if left amino acid is an H, check if it's linked to the current amino acid
                                # if it is an unconsidered H that isn't linked to the current amino acid, it is an H-H connection
                                if not (current_aa[2] == left_aa[0] or current_aa[3] == left_aa[0]):
                                    score = score + 1  # increment score
                                    # add both amino acids to the list of considered amino acids
                                    considered_aas.append(current_aa[0])
                                    considered_aas.append(left_aa[0])
                        # if above amino acid hasn't been considered, check if it is an H
                        if not is_above_aa_considered and not above_aa[1] == "None":
                            if amino_dict[above_aa[1]] == 'H':
                                # if above amino acid is an H, check if it's linked to the current amino acid
                                # if it is an unconsidered H that isn't linked to the current amino acid, it is an H-H connection
                                if not (current_aa[2] == above_aa[0] or current_aa[3] == above_aa[0]):
                                    score = score + 1  # increment score
                                    # add both amino acids to the list of considered amino acids
                                    considered_aas.append(current_aa[0])
                                    considered_aas.append(above_aa[0])
                        # if below amino acid hasn't been considered, check if it is an H
                        if not is_below_aa_considered and not below_aa[1] == "None":
                            if amino_dict[below_aa[1]] == 'H':
                                # if below amino acid is an H, check if it's linked to the current amino acid
                                # if it is an unconsidered H that isn't linked to the current amino acid, it is an H-H connection
                                if not (current_aa[2] == below_aa[0] or current_aa[3] == below_aa[0]):
                                    score = score + 1  # increment score
                                    # add both amino acids to the list of considered amino acids
                                    considered_aas.append(current_aa[0])
                                    considered_aas.append(below_aa[0])

    return score


def create_hp_model(filename):
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

    return 0


if __name__ == '__main__':
    main()

# Simulates DNA sequence partitioning
# Inputs:
#      string filename: input file name

# Sample call: python genetic_algorithm.py "SeaCucumberGlobin.txt"

import numpy
import sys
import random
from matplotlib import pyplot


def main():
    if len(sys.argv) < 2:  # Verifying there aren't too few inputs
        print("Error: Too few inputs.\n")
        sys.exit(1)
    if len(sys.argv) > 2:  # Verifying there aren't too many inputs
        print("Error: Too many inputs.\n")
        sys.exit(1)

    try:  # Verifying input of filename is correct and storing it in variable filename
        filename = sys.argv[1]
    except ValueError:
        print("Error: input for filename must be of type string.")
        sys.exit(1)

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

    size = 100
    population = create_initial_population(aa_sequence, size)
    epoch = 0
    average_score = 0
    avg_scores = []
    max_score = 0
    when_to_break = 0

    # This while loop evolves the fold sequences
    while epoch < 27:

        scores = []
        parents = []

        for i in range(len(population)):
            scores.append(population[i][2])
            parents.append(population[i][0])

        breeders = select_parents(scores, parents)
        children = []
        for i in range(0, len(population), 1):
            child = create_child(breeders, aa_sequence)
            children.append(child)

        # print(numpy.average(scores))

        if max_score < average_score:
            max_score = average_score
            when_to_break = 0
        else:
            when_to_break = when_to_break + 1

        # If we have gotten to the maximum score for 10 iterations we stop running the algorithm.
        if when_to_break == 10:
            print("The maximum score is: " + str(max_score))
            break

        epoch += 1
        population.clear()
        population = children
        average_score = sum(scores) / len(scores)
        print("The score of the fold for epoch " + str(epoch) + " is: ")
        print(average_score)
        avg_scores.append(average_score)
    plot(avg_scores, max_score)

    return 0


# method to generate the initial population
def create_initial_population(aa_sequence, size):
    # Create the population of directions; F,L,and R.
    initial_population = []  # variable to store the initial population
    while len(initial_population) < size:
        member = []  # variable to represent a member of the initial population in the form [traverse_sequence, d2_lattice, score]
        d2_lattice = numpy.ndarray([2 * size + 2, 2 * size + 2], dtype=object)
        for r in range(0, len(d2_lattice), 1):
            for c in range(0, len(d2_lattice[r]), 1):
                # format is: [index, amino acid, previous index, next index]
                d2_lattice[r, c] = ["None", "None", "None", "None"]
        r_index = int(size)
        c_index = int(size)  # indices in the d2_lattice matrix
        orientation = 'E'  # current orientation in units of cardinal direction
        traverse_sequence = []  # sequence of traversal moves
        bump = False  # flag for if there's a bump
        for i in range(0, len(aa_sequence), 1):
            bump = True
            # insert amino acids into the positions they should go
            if i == 0:
                d2_lattice[r_index, c_index] = [i, aa_sequence[i], "None", i + 1]
            elif i == len(aa_sequence) - 1:
                d2_lattice[r_index, c_index] = [i, aa_sequence[i], i - 1, "None"]
            else:
                d2_lattice[r_index, c_index] = [i, aa_sequence[i], i - 1, i + 1]

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
    # traverse d2_lattice matrix
    for r in range(0, len(d2_lattice), 1):
        for c in range(0, len(d2_lattice), 1):
            current_aa = d2_lattice[r, c]
            # if the current index contains an amino acid, check for hydrophobia
            if not current_aa[0] == "None":
                if amino_dict[current_aa[1]] == 'H':
                    # check for nearby H's
                    above_aa = d2_lattice[r - 1, c]
                    left_aa = d2_lattice[r, c - 1]
                    right_aa = d2_lattice[r, c + 1]
                    below_aa = d2_lattice[r + 1, c]
                    # check if right amino acid is an H
                    if not right_aa[1] == "None":
                        if amino_dict[right_aa[1]] == 'H':
                            # if right amino acid is an H, check if it's linked to the current amino acid
                            # if it is an H that isn't linked to the current amino acid, it is an H-H connection
                            if not (current_aa[2] == right_aa[0] or current_aa[3] == right_aa[0]):
                                score = score + 1  # increment score
                    # check if left amino acid is an H
                    if not left_aa[1] == "None":
                        if amino_dict[left_aa[1]] == 'H':
                            # if left amino acid is an H, check if it's linked to the current amino acid
                            # if it is an H that isn't linked to the current amino acid, it is an H-H connection
                            if not (current_aa[2] == left_aa[0] or current_aa[3] == left_aa[0]):
                                score = score + 1  # increment score
                    # check if above amino acid is an H
                    if not above_aa[1] == "None":
                        if amino_dict[above_aa[1]] == 'H':
                            # if above amino acid is an H, check if it's linked to the current amino acid
                            # if it is an H that isn't linked to the current amino acid, it is an H-H connection
                            if not (current_aa[2] == above_aa[0] or current_aa[3] == above_aa[0]):
                                score = score + 1  # increment score
                    # check if below amino acid is an H
                    if not below_aa[1] == "None":
                        if amino_dict[below_aa[1]] == 'H':
                            # if below amino acid is an H, check if it's linked to the current amino acid
                            # if it is an H that isn't linked to the current amino acid, it is an H-H connection
                            if not (current_aa[2] == below_aa[0] or current_aa[3] == below_aa[0]):
                                score = score + 1  # increment score
    return score


def select_parents(scores, parents):
    # parents is a list that stores lists of parent path sequences
    # scores is a list that stores the score of each parent.

    # Normalize the scores and place them in a list; normalized scores
    normalized_scores = []
    sum_scores = sum(scores)
    for i in range(len(scores)):
        normalized_scores.append(scores[i] / sum_scores)

    # We need to create a population pool based off of the normalized scores. A parents score is proportional to the
    # number of copies of itself in the mating pool In this example parent1 has a 15% probability of being a parent.
    # If we initialize the mating pool to 100. 15 of the parents will be parent1
    pool_size = 100
    pool = []

    for i in range(len(normalized_scores)):
        number_in_pop = int(normalized_scores[i] * pool_size)
        for n in range(number_in_pop):
            pool.append(parents[i])

    # Now the pool contains the parents in their corresponding proportions relative to the pool size
    # Randomize the pool of parents and then randomly select two parents.
    # random.shuffle(pool)
    parent_one_index = 1
    parent_two_index = 1

    # This code allows us to have two unique parents
    while pool[parent_one_index] == pool[parent_two_index]:
        parent_one_index = random.randint(0, len(pool) - 1)
        parent_two_index = random.randint(0, len(pool) - 1)

    selected_parents = [pool[parent_one_index], pool[parent_two_index]]
    return selected_parents


# After selecting the parents we create a child using these two parents
def create_child(selected_parents, aa_sequence):
    # spawn is the child of the selected_parents.
    spawn = []

    # Use the 50-50 method to create a new child
    midpoint = int(len(selected_parents[0]) / 2)
    first_half = selected_parents[0][0:midpoint]
    second_half = selected_parents[1][midpoint:len(selected_parents[1])]

    child_traversal_sequence = first_half + second_half
    child = mutate(child_traversal_sequence, aa_sequence)
    return child


# We have to have a function that mutates the children to introduce variation into the population.
def mutate(child_traversal_sequence, aa_sequence):
    # We have to have a mutation rate as well.
    # 1% mutation rate
    mutation_rate = 1
    size = len(child_traversal_sequence)
    child_d2_lattice = numpy.ndarray([2 * size + 2, 2 * size + 2], dtype=object)
    score = 0
    for r in range(0, len(child_d2_lattice), 1):
        for c in range(0, len(child_d2_lattice[r]), 1):
            # format is: [index, amino acid, previous index, next index]
            child_d2_lattice[r, c] = ["None", "None", "None", "None"]
    for i in range(len(child_traversal_sequence)):
        random_number = random.randint(0, 100)
        letter = ['F', 'L', 'R']
        random_letter_index = random.randint(0, 2)
        if random_number <= mutation_rate:
            child_traversal_sequence[i] = letter[random_letter_index]

    r_index = int(size)
    c_index = int(size)  # indices in the d2_lattice matrix
    orientation = 'E'  # current orientation in units of cardinal direction
    for i in range(0, len(child_traversal_sequence), 1):
        if i == 0:
            child_d2_lattice[r_index, c_index] = [i, aa_sequence[i], "None", i + 1]
        elif i == len(aa_sequence) - 1:
            child_d2_lattice[r_index, c_index] = [i, aa_sequence[i], i - 1, "None"]
        else:
            child_d2_lattice[r_index, c_index] = [i, aa_sequence[i], i - 1, i + 1]

        current_traversal = child_traversal_sequence[i]
        current_traversal_info = get_traversal_info(current_traversal, orientation)
        r_index = r_index + current_traversal_info[0]
        c_index = c_index + current_traversal_info[1]
        if not child_d2_lattice[r_index, c_index][0] == "None":
            score = score - 2000
        orientation = current_traversal_info[2]
        # if child_d2_lattice[new_r_index, new_c_index][0] == "None":
        #     orientation = current_traversal_info[2]
        #     r_index = new_r_index
        #     c_index = new_c_index
    # print("Mutation complete")
    # print(score)
    return [child_traversal_sequence, child_d2_lattice, score_lattice(child_d2_lattice)]


def plot(avg_scores, max_score):
    figure = pyplot.figure()
    font = "Comic Sans MS"
    graph_area = figure.add_subplot(111)
    x_values = []
    for i in range(0, len(avg_scores), 1):
        x_values.append(i)
    graph_area.get_xaxis().set_ticks(x_values)
    pyplot.ylim(0, 100)
    pyplot.xlim(1, len(avg_scores))
    pyplot.plot(x_values, avg_scores)
    pyplot.axhline(y=max_score, color="red", label="Maximum Score", font=font, linestyle="-")
    figure.suptitle('Average Lattice Scores (Pop. = 100)', font=font)
    pyplot.xlabel('Epoch', font=font)
    pyplot.ylabel('Avg. Score', font=font)

    pyplot.show()
    return 0


if __name__ == '__main__':
    main()

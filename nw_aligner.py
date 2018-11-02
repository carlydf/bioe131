#!/usr/bin/python

"""

Needleman-Wunsch Aligner
Bioengineering 131/231, Fall 2018

Command-line script to read in the contents of a multi-FASTA containing two
sequences and a score matrix, perform a global alignment, and print the
resulting alignment matrix and optimal alignment to STDOUT.

"""

import os
import sys

class NWAligner:
    def __init__(self, score_matrix_fname):
        self.score_matrix, self.gap_penalty = self.load_score_matrix(score_matrix_fname)

    @staticmethod
    def load_score_matrix(fname):
        """
        Input: (String) A path to a scoring matrix file.
        Output: (Dictionary) A nested dictionary where keys are strings
                and elements are scores as integers.
    
        Example:
    
        >>> matrix, gap_penalty = NWAligner.load_score_matrix('/home/bioe131/BLOSUM62')
        >>> matrix['A']['A']
        4
        >>> matrix['W']['W']
        11
        >>> gap_penalty
        -4

        """

        score_matrix = {}
        gap_penalty = None
        
        
        """with open(fname) as fp:
            for line_num, line in enumerate(fp):
                # ignore comments in matrix file
                if line.startswith("#"):
                    continue"""

                # Parse matrix file line-by-line and load into nested dictionaries.
                #
                # Last line of matrix contains the gap penalty which must be pulled
                # out and returned.
        file = open(fname).read()
                
        lst1 = file.split()
        gap_penalty = [lst1[len(lst1)-1]]
    
        count = 0
        cols = int(sqrt(len(lst1)))
        rows = cols
        for r in range(rows):
            rkey = lst1[r]
            #gives each row in the scoring matrix an empty dictionary to be filled with scores keyed by ckey
            score_matrix[rkey] = {}
            for c in range(cols):
                ckey = lst1[c]
                #adds a key to the nested dictionary for each column, that points to that columns's score in row r
                score_matrix[rkey][ckey] = lst1[cols*(1+r) + c]
        #adds every row to the matrix
            #lst = [lst1[x:x+23] for x in range(0, len(lst1), 23)]
        #letterlist = list(lst[0])
        #lst.remove(lst[0])
        
        #score_matrix = dict(zip(letterlist, lst))
        
        return score_matrix, gap_penalty

    @staticmethod
    def load_FASTA(fname):
        """
        Input: (String) A path to a FASTA file containing exactly two sequences.
        Output: (List) A list containing two strings: one for each sequence.

        Example:

        >>> seqs = NWAligner.load_FASTA('example.fa')
        >>> seqs[0]
        'YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        >>> seqs[1]
        'WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        >>> len(seqs)
        2

        """
        Fasta_file = open(fname).read()
        seqs = Fasta_file.split()
        seqs.remove(seqs[0])
        seqs.remove(seqs[1])
        
        if len(seqs) > 2:
            raise ValueError("list must have no more then 2 sequences")
        ### TODO ###
        # Load FASTA file and return list of sequences.
        # Throw an error if there are more than two sequences in the file.

        return seqs

    def align(self, seq_x, seq_y, print_matrix = False):
        """
        Input: (Strings) Two sequences to be aligned (seq_x and seq_y).
               (Boolean) If print_matrix is True, print the dynamic programming
                         matrix before traceback.
        Output: (Tuple of strings) Two sequences, aligned.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = aligner.load_FASTA('example.fa')
        >>> aligner.align(seqs[0], seqs[1])
        ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---',
         'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')

        """

        ###
        ### INITIALIZATION
        ###

        # create two empty matrices with sizes based on the input sequences.
        # one contains the dynamic programming matrix, the other contains
        # pointers we'll use during traceback
        matrix = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        pointers = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]

        ### TODO ###
        #scoring scheme:
        #match = +5, mismatch = -5, indel = -10
        rows = len(seq_x) + 1
        cols = len(seq_y) + 1
        up = 'U'
        left = 'L'
        diag = 'D'
        start = 'S'
        # Fill the top row of the matrix with scores for gaps. indexing matrix[r][c]
        for i in range(rows):
            matrix[i][0] = i*(-10) #set first column to 0, -10, -20,...
            pointers[i][0] = up
        # Fill the first column of the matrix with scores for gaps.
        for i in range(cols):
            matrix[0][i] = i*(-10) #set first row to 0, -10, -20,...
            pointers[0][i] = left
            
        pointers[0][0] = start

        ###
        ### RECURSION
        ###

        # fill the dynamic programming and pointer matrices
        for x in range(1, rows):
            for y in range(1, cols):
                match_score = int(self.score_matrix[seq_x[x - 1]][seq_y[y - 1]])
                gap_penalty = int(self.gap_penalty[0])
                qd = matrix[x-1][y-1] + match_score
                qu = matrix[x-1][y] + gap_penalty
                ql = matrix[x][y-1] + gap_penalty
                
                #tracking pointers:
                qs = [qd, qu, ql]
                pointer_options = [diag, up, left]
                #puts the max value in our matrix
                matrix[x][y] = max(qs)
                #finds the index of the max value
                max_index = qs.index(max(qs))
                #gets the desired pointer value from pointer_options and puts it in the pointers matrix
                pointers[x][y] = pointer_options[max_index]

                ### TODO ###
                # Take the maximum score of three possibilities:
                #   1) The element in the matrix diagonal from this one
                #      plus the score of an exact match
                #   2) The element to the left plus a gap penalty
                #   3) The element above plus a gap penalty
                # ... and set the current element (matrix[x][y]) equal to that
                #
                # Keep track of which of these choices you made by setting
                #   the same element (i.e., pointers[x][y]) to some value that
                #   has meaning to you.

        # print the dynamic programming matrix
        if print_matrix:
            for x in range(len(seq_x) + 1):
                print(" ".join(map(lambda i: str(int(i)), matrix[x])))

        ###
        ### TRACEBACK
        ###

        # starting from the bottom right corner, follow the pointers back
        x, y = len(seq_x), len(seq_y)

        # fill these lists with the aligned sequences
        align_x = []
        align_y = []

        while x > 0 or y > 0:
            move = pointers[x][y]
            if move == diag:
                align_x += [seq_x[x-1]]
                align_y += [seq_y[y-1]]
                x -= 1
                y -= 1
            elif move == up:
                align_x += [seq_x[x-1]]
                align_y += ['-']
                x -= 1
            elif move == left:
                align_x += ['-']
                align_y += [seq_y[y-1]]
                y -= 1

            ### TODO ###
            # Follow pointers back through the matrix to the origin.
            # Depending on which "move" you made at each element in the
            #   matrix, you'll either align seq_x to seq_y, seq_x to a gap, or
            #   seq_y to a gap.

        # flip the alignments, as they're reversed
        return ("".join(align_x[::-1]), "".join(align_y[::-1]))

###                                      ###
### NO NEED TO EDIT CODE BELOW THIS LINE ###
###                                      ###

"""if __name__ == '__main__':
    def usage():
        print('usage: %s matrixfilename stringfilename')
        sys.exit(1)

    if len(sys.argv) != 3:
        usage()

    for fname in sys.argv[1:]:
        if not os.path.isfile(fname):
            print('Can not open %s' % (fname,))
            usage()

    aligner = NWAligner(sys.argv[1])
    seqs = aligner.load_FASTA(sys.argv[2])
    result = aligner.align(seqs[0], seqs[1])

    print('>seq1\n%s\n>seq2\n%s' % (result[0], result[1]))
"""
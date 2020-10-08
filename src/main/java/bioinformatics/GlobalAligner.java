package bioinformatics;

/**
 * Global Aligner class
 */
public class GlobalAligner {

    private final String sequence1, sequence2;
    private final SubstitutionMatrix substitutionMatrix;
    private final int[][] scoreMatrix;
    private final int gapPenalty;
    private final int alignmentScore;

    GlobalAligner(String sequence1, String sequence2, SubstitutionMatrix substitutionMatrix, int gapPenalty) {
        this.sequence1 = sequence1;
        this.sequence2 = sequence2;
        this.substitutionMatrix = substitutionMatrix;
        this.scoreMatrix = new int[sequence1.length()+1][sequence2.length()+1];
        this.gapPenalty = gapPenalty;

        initializeScoreMatrix();
        this.alignmentScore = scoreMatrix[scoreMatrix.length-1][scoreMatrix[0].length-1];

    }

    public String getSequence1() { return this.sequence1; }

    public String getSequence2() { return this.sequence2; }

    public SubstitutionMatrix getSubstitutionMatrix() { return this.substitutionMatrix; }

    public int[][] getScoreMatrix() { return this.scoreMatrix; }

    public int getGapPenalty() { return gapPenalty; }

    public int getAlignmentScore() { return this.alignmentScore; }

    /**
     * Calculates the optimal score of a single position in the score matrix based
     * on Global Alignment recurrence relation
     * @param i,j scoreMatrix indices, integer
     * @param a,b amino acid base pair, String
     * @return optimal score, integer
     */
    private int calculateScore(int i, int j, String a, String b) {
        // diagonal score C(i-1, j-1) + w(S_i, T_j)
        int score = scoreMatrix[i-1][j-1] + substitutionMatrix.getScore(a, b);
        // left score C(i, j-1) - gapPenalty
        int leftScore = scoreMatrix[i][j-1] + gapPenalty;
        // top score C(i-1, j) - gapPenalty
        int topScore = scoreMatrix[i-1][j] + gapPenalty;
        // find max of scores
        if (leftScore > score)
            score = leftScore;
        if (topScore > score)
            score = topScore;
        return score;
    }

    private void initializeScoreMatrix() {
        for (int i = 0; i < scoreMatrix.length; i++) {
            for (int j = 0; j < scoreMatrix[0].length; j++) {
                if (i < 1) {
                    scoreMatrix[i][j] = j * gapPenalty;
                }
                else if (j < 1) {
                    scoreMatrix[i][j] = i * gapPenalty;
                }
                else {
                    String base1 = String.valueOf(sequence1.charAt(i-1));
                    String base2 = String.valueOf(sequence2.charAt(j-1));
                    scoreMatrix[i][j] = calculateScore(i, j, base1, base2);
                }
            }
        }
    }

    /**
     * Calculates the optimal aligned sequences using highroad traceback to construct
     * 3 alignment strings: sequence 1 with gaps, matches denoted by '|', sequence 2 with gaps
     * @return array of aligned sequences [sequence1, matches, sequence2], String array
     */
    public String[] performHighroadTraceback() {
        String[] optimalAlignmentStrings = new String[] {"", "", ""};
        // begin at bottom right of scoreMatrix
        int i = scoreMatrix.length - 1;
        int j = scoreMatrix[0].length - 1;
        // trace until no elements of the sequences remain
        while (i > 0 || j > 0) {
            // sequence 1 base
            String base1 = String.valueOf(sequence1.charAt(i-1));
            // sequence 2 base
            String base2 = String.valueOf(sequence2.charAt(j-1));
            // retrieve scores from adjacent indices
            int score = scoreMatrix[i][j];
            int diag = scoreMatrix[i-1][j-1] + substitutionMatrix.getScore(base1, base2);
            int left = scoreMatrix[i][j-1] + gapPenalty;
            int top = scoreMatrix[i-1][j] + gapPenalty;
            // if top score matches, add base 1 to aligned sequence 1 and add gap to aligned sequence 2
            if (score == top) {
                optimalAlignmentStrings[0] = base1 + optimalAlignmentStrings[0];
                optimalAlignmentStrings[2] = "-" + optimalAlignmentStrings[2];
                optimalAlignmentStrings[1] = " " + optimalAlignmentStrings[1];

                i--;
            }
            // if diagonal score matches, add both bases to aligned strings and check for match
            else if (score == diag) {
                optimalAlignmentStrings[0] = base1 + optimalAlignmentStrings[0];
                optimalAlignmentStrings[2] = base2 + optimalAlignmentStrings[2];
                if (base1.equals(base2))
                    optimalAlignmentStrings[1] = "|" + optimalAlignmentStrings[1];
                else
                    optimalAlignmentStrings[1] = " " + optimalAlignmentStrings[1];

                i--;
                j--;
            }
            // otherwise add gap to aligned sequence 1 and add base 2 to aligned sequence 2
            else {
                optimalAlignmentStrings[0] = "-" + optimalAlignmentStrings[0];
                optimalAlignmentStrings[2] = base2 + optimalAlignmentStrings[2];
                optimalAlignmentStrings[1] = " " + optimalAlignmentStrings[1];

                j--;
            }
        }

        return optimalAlignmentStrings;
    }

}

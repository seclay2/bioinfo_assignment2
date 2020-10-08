package bioinformatics;

import java.io.File;
import java.util.Scanner;

/**
 * Assignment 2: Global alignment with linear gap penalty
 *
 * Requires 3 command line arguments:
 *  - Sequence filename
 *  - Substitution matrix filename
 *  - Gap penalty value
 *
 * Performs global alignment on 2 DNA or Protein sequences using Highroad Traceback.
 * Prints formatted results.
 *
 */
public class App
{
    public static void main( String[] args )
    {
        // Check for appropriate arguments
        if (args.length < 3) {
            System.out.println("Missing required arguments: [sequence file] [substitution matrix file] [gap penalty]");
        }
        else {

            File sequenceFile = new File(args[0]);
            File matrixFile = new File(args[1]);
            int gapPenalty = Integer.parseInt(args[2]);

            try {
                // parse sequence file to get 2 separate sequences
                Scanner scanner = new Scanner(sequenceFile);
                scanner.useDelimiter("\n\n");
                String sequence1 = scanner.next();
                String sequence2 = scanner.next();
                scanner.close();

                // create substitution matrix object
                SubstitutionMatrix substitutionMatrix = new SubstitutionMatrix(matrixFile);

                // create Global Aligner object
                GlobalAligner globalAligner = new GlobalAligner(sequence1, sequence2, substitutionMatrix, gapPenalty);

                // get aligned sequences based on highroad traceback
                String[] alignmentStrings = globalAligner.performHighroadTraceback();

                // Print results
                System.out.println("Global Alignment");
                if (sequence1.toLowerCase().matches("[actg]*")) {
                    System.out.println("DNA Sequence");
                } else {
                    System.out.println("Protein Sequence");
                }
                System.out.println(globalAligner.getGapPenalty());
                System.out.println(alignmentStrings[0]);
                System.out.println(alignmentStrings[1]);
                System.out.println(alignmentStrings[2]);
                System.out.println(globalAligner.getAlignmentScore());
            }
            catch (Exception e) {
                System.out.println("An error occurred.");
                e.printStackTrace();
            }
        }
    }

}


import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.Random;


public class Main {

    //number of Monte Carlo trials in sample
    public static final long MC = 2000;
    //learning rate for Boltzmann learning
    public static final double ETA = .1;

    /**
     * main function
     *
     * @param args
     */
    public static void main(String args[]) {
        Scanner scan;
        try {
            scan = new Scanner(new File("cytochrome.fasta"));
        } catch (FileNotFoundException e) {
            System.out.println(e.toString());
            return;
        }
        ArrayList<String> alignments = processMSA(scan);
        String consensus = "--------GDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSY";
        ArrayList<Byte[]> bitAlignments = new ArrayList<>();
        bitAlignments.add(proteinToBitString(consensus, consensus));
        alignments.stream().forEach((s) -> {
            bitAlignments.add(proteinToBitString(consensus, s));
        });
        double h[] = new double[consensus.length()];
        double J[][] = new double[consensus.length()][consensus.length()];
        double singleProbs[] = calculateSingleMutationalProbabilities(bitAlignments, consensus.length());
        System.out.println("Done with single probs");
        double doubleProbs[][] = calculateDoubleMutationalProbabilities(bitAlignments, consensus.length());
        System.out.println("Done with double probs");
        HamiltonianConstants hc = boltzmannLearn(singleProbs, doubleProbs, consensus.length(), h, J);
        
        //this is just a toy example that has nothing to do with HIV
        System.out.println("TUNA: " + evaluateHamiltonian(bitAlignments.get(1), hc.getH(), hc.getJ()));
        System.out.println("GRAY WHALE: " + evaluateHamiltonian(bitAlignments.get(2), hc.getH(), hc.getJ()));
        System.out.println("SNAPPING TURTLE: " + evaluateHamiltonian(bitAlignments.get(3), hc.getH(), hc.getJ()));
        System.out.println("RHESUS MONKEY: " + evaluateHamiltonian(bitAlignments.get(4), hc.getH(), hc.getJ()));
        //this is just to test that I am correctly calculating energy
    }

    /**
     * Process MSA (Multiple Sequence Alignment)
     *
     * This processes an MSA from an input stream in FASTA format. First, it
     * essentially just converts the stream MSA into a list of sequences
     *
     * @param scan Scanner object defining the input stream to use (usually from
     * a file)
     * @return ArrayList<String> containing each sequence from the MSA in order
     */
    public static ArrayList<String> processMSA(Scanner scan) {
        ArrayList<String> out = new ArrayList<>();
        String curString = "";
        while (scan.hasNextLine()) {
            String a = scan.nextLine();
            if (a.substring(0, 1).equals(">")) {
                out.add(curString);
                curString = "";
            } else {
                curString += removeWhitespace(a);
            }
        }
        out.remove(0);
        return out;
    }

    /**
     * Remove Whitespace function
     *
     * Helper function for processMSA to remove whitespace
     *
     * @param a String from which whitespace must be removed
     * @return String without whitespace
     */
    private static String removeWhitespace(String a) {
        String out = "";
        for (int i = 0; i < a.length(); i++) {
            if (!a.substring(i, i + 1).equals(" ") && !a.substring(i, i + 1).equals("\n") && !a.substring(i, i + 1).equals("\t")) {
                out += a.substring(i, i + 1);
            }
        }
        return out;
    }

    /**
     * Protein to Bit String Function
     *
     * Converts an amino acid sequence into a binary number using a consensus
     * sequence. Converts each amino acid to {0,1}. 0 denotes the "wild-type"
     * sequence detailed in the consensus. 1 denotes a mutation of any amino
     * acid or a gap (-) in the sequence.
     *
     * @param consensus Consensus sequence denoting the "wild-type"
     * @param protein Mutated protein
     * @return Byte[] of 1's and 0's representing the protein (there is no bit
     * class and booleans can be unhelpful for math)
     */
    public static Byte[] proteinToBitString(String consensus, String protein) {
        Byte[] a = new Byte[consensus.length()];
        for (int i = 0; i < consensus.length(); i++) {
            a[i] = booleanToByte(!consensus.substring(i, i + 1).equals(protein.substring(i, i + 1)));
        }
        return a;
    }

    /**
     * Boolean to Byte Function
     *
     * Helper method for Protein to BitString
     *
     * @param a boolean
     * @return byte
     */
    private static Byte booleanToByte(boolean a) {
        return a ? new Byte((byte) 1) : new Byte((byte) 0);
    }

    /**
     * Evaluate Hamiltonian Function
     *
     * Evaluates the Ising model Hamiltonian given its input. The Ising model
     * Hamiltonian is used in the fitness model for HIV. Given constant vector h
     * and matrix J, the Hamiltonian evaluates a function for a bit string x
     * that, if constants are accurate, should accurately be correlated with the
     * mutational probabilities of each amino acid.
     *
     * @param s BitString protein
     * @param h Lagrangian multiplier
     * @param J Lagrangian multiplier
     * @return
     */
    public static double evaluateHamiltonian(Byte[] s, double[] h, double[][] J) {
        double h_sum = 0;
        for (int i = 0; i < s.length; i++) {
            h_sum += s[i] * h[i];
        }
        double J_sum = 0;
        for (int i = 0; i < s.length; i++) {
            for (int j = i + 1; j < s.length; j++) {
                J_sum += s[i] * s[j] * J[i][j];
            }
        }
        return h_sum + J_sum;
    }

    /**
     * This class encapsulates the Hamiltonian constant vectors h and J.
     */
    public static class HamiltonianConstants {

        private final double[] h;
        private final double[][] J;

        public HamiltonianConstants(double[] h, double[][] J) {
            this.h = h;
            this.J = J;
        }

        public double[] getH() {
            return h;
        }

        public double[][] getJ() {
            return J;
        }
    }

    /**
     * Boltzmann Learning function
     *
     * Generates constant vectors h_i and J_ij for an MSA. The process is to run
     * through the Hamiltonian and in each step, edit h and J such that they
     * predict closer to the observed mutational probabilities
     *
     * @param singleProbs Observed mutational probabilities at each locus
     * @param doubleProbs Observed double mutational probabilities at each pair
     * of loci
     * @param length
     * @param h Initial value for h, the constant vector
     * @param J Initial value for J, the constant vector
     * @return HamiltonianConstants object containing corrected constant values
     */
    public static HamiltonianConstants boltzmannLearn(double[] singleProbs, double[][] doubleProbs, int length, double[] h, double[][] J) {
        int count = 0;
        double[] delta_h = new double[h.length];
        double[][] delta_J = new double[J.length][J[0].length];
        for (int i = 0; i < length; i++) {
            delta_h[i] = ETA * (singleProbs[i]);
            for (int j = i + 1; j < length; j++) {
                delta_J[i][j] = ETA * (doubleProbs[i][j]);
            }
        }
        for (int i = 0; i < h.length; i++) {
            h[i] += delta_h[i];
        }
        for (int i = 0; i < length; i++) {
            for (int j = i + 1; j < length; j++) {
                J[i][j] = J[i][j] + delta_J[i][j];
            }
        }
        while (count++ < 5) { //TODO: make it so this only terminates when delta_h and delta_J are small enough
            for (int i = 0; i < length; i++) {
                delta_h[i] = ETA * (singleProbs[i] - thermalAverage(i, length, h, J));
                for (int j = i + 1; j < length; j++) {
                    delta_J[i][j] = ETA * (doubleProbs[i][j] - thermalAverage(i, j, length, h, J));
                }
            }
            for (int i = 0; i < h.length; i++) {
                h[i] += delta_h[i];
            }
            for (int i = 0; i < length; i++) {
                for (int j = i + 1; j < length; j++) {
                    J[i][j] = J[i][j] + delta_J[i][j];
                }
            }
            System.out.println("Learning step: "+count);
        }
        return new HamiltonianConstants(h, J);
    }

    /**
     * Thermal Averaging function (one position)
     *
     * Finds the average over the Ising distribution of a position. This finds,
     * essentially, a weighted sum of all possible sequences in the Hamiltonian
     * in which the the position specified is a mutation. This is extremely
     * computationally intensive, so we instead take a Monte Carlo sample of the
     * possible sequences (20,000 trials) and take that average instead. This
     * sacrifices accuracy but gains efficiency.
     *
     * @param position position of amino acid that is mutated
     * @param length length of bit string protein (number of amino acids)
     * @param h Hamiltonian constant h
     * @param J Hamiltonian constant J
     * @return The thermal average over the Ising distribution
     */
    public static double thermalAverage(int position, int length, double[] h, double[][] J) {
        double sum = 0;
        Random r = new Random();
        for (int i = 0; i < MC; i++) {
            Byte[] rand = new Byte[length];
            for (int j = 0; j < length; j++) {
                rand[j] = r.nextBoolean() ? new Byte((byte) 1) : new Byte((byte) 0);
            }
            rand[position] = 1;
            sum += Math.exp(-evaluateHamiltonian(rand, h, J));
        }

        return sum / MC;
    }

    /**
     * Thermal Averaging function (two positions)
     *
     * Finds the average over the Ising distribution of a position. This finds,
     * essentially, a weighted sum of all possible sequences in the Hamiltonian
     * in which the the position specified is a mutation. This is extremely
     * computationally intensive, so we instead take a Monte Carlo sample of the
     * possible sequences (20,000 trials) and take that average instead. This
     * sacrifices accuracy but gains efficiency.
     * @param position1 position of one mutated amino acid
     * @param position2 position of other mutated amino acid
     * @param size length of bit string protein (number of amino acids)
     * @param h Hamiltonian constant h
     * @param J Hamiltonian constant J
     * @return The thermal average over the Ising distribution
     */
    public static double thermalAverage(int position1, int position2, int size, double[] h, double[][] J) {
        double sum = 0;
        Random r = new Random();
        for (int i = 0; i < MC; i++) {
            Byte[] rand = new Byte[size];
            for (int j = 0; j < size; j++) {
                rand[j] = r.nextBoolean() ? new Byte((byte) 1) : new Byte((byte) 0);
            }
            rand[position1] = 1;
            rand[position2] = 1;
            sum += Math.exp(-evaluateHamiltonian(rand, h, J));
        }
        return sum / MC;
    }

    /**
     * Calculates mutational probabilities by averaging them over the MSA
     * 
     * @param msa       Bytestring MSA
     * @param length    length of each
     * @return          Probability of each mutation
     */
    public static double[] calculateSingleMutationalProbabilities(ArrayList<Byte[]> msa, int length) {
        double[] avgs = new double[length];
        msa.stream().forEach((b) -> {
            for (int i = 0; i < length; i++) {
                avgs[i] += b[i].byteValue();
            }
        });
        for (int i = 0; i < length; i++) {
            avgs[i] /= msa.size();
        }

        return avgs;
    }

    /**
     * Calculates double mutational probabilities by averaging them over the MSA
     * 
     * @param msa       Bytestring MSA
     * @param length    length of each
     * @return          Probability of each double mutation
     */
    public static double[][] calculateDoubleMutationalProbabilities(ArrayList<Byte[]> msa, int length) {
        double[][] avgs = new double[length][length];
        msa.stream().forEach((b) -> {
            for (int i = 0; i < length; i++) {
                for (int j = i + 1; j < length; j++) {
                    avgs[i][j] += b[i].byteValue() * b[j].byteValue();
                }
            }
        });
        for (int i = 0; i < length; i++) {
            for (int j = i + 1; j < length; j++) {
                avgs[i][j] /= msa.size();
            }
        }
        return avgs;
    }
}

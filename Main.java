
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Main {

    //number of Boltzmann learning steps
    public static final long BOLTZ = 10;
    //learning rate for Boltzmann learning
    public static final double ETA = 2;
    //number of steps in walking for MMC
    public static final long MC_WALK = 100000;
    //number of rejections per sample for MMC
    public static final long MC_RATE = 500;
    //number of Monte Carlo trials in thermalAverage; not currently in use
    public static final long MC = 10000;
    //number of Monte Carlo trials in partition
    public static final long MC_PARTITION = 100000; //134217728;

    //value by which I have to multiply the score in order to convert to energy
    public static double score_multiplier;

    /**
     * main function
     *
     * @param args
     */
    public static void main(String args[]) {
        // // the below code scans a file for resistance mutation scores
        // // <editor-fold defaultstate="collapsed" desc="Scan File">
        // Scanner fileScanner;
        // try {
        //     fileScanner = new Scanner(new File("nrti.txt"));
        // } catch (FileNotFoundException e) {
        //     System.out.println(e.toString());
        //     return;
        // }
        // // </editor-fold>
        //
        // // the below code processes the scanner made above to create a list of resistance mutations
        // // <editor-fold defaultstate="collapsed" desc="Make Resistance Mutation Database">
        // ArrayList<ResistanceMutations> scores = new ArrayList<>();
        // while (fileScanner.hasNext()) {
        //     int position = Integer.parseInt(fileScanner.nextLine());
        //     char cons = fileScanner.nextLine().charAt(0);
        //     char aa = fileScanner.nextLine().charAt(0);
        //     byte threetc = Byte.parseByte(fileScanner.nextLine());
        //     byte ftc = Byte.parseByte(fileScanner.nextLine());
        //     byte abc = Byte.parseByte(fileScanner.nextLine());
        //     byte azt = Byte.parseByte(fileScanner.nextLine());
        //     byte d4t = Byte.parseByte(fileScanner.nextLine());
        //     byte ddi = Byte.parseByte(fileScanner.nextLine());
        //     byte tdf = Byte.parseByte(fileScanner.nextLine());
        //     if (aa=='i' || aa=='d') continue;
        //     scores.add(new ResistanceMutations(position, cons, aa, threetc, ftc, abc, azt, d4t, ddi, tdf));
        // }
        // // </editor-fold>
        //
        // System.out.println("pos\tcons\taa\t3tc\tftc\tabc\tazt\td4t\tddi\ttdf");
        // for (ResistanceMutations rm : scores) {
        //     System.out.println(rm);
        // }
        // constantGeneration(args);

        // <editor-fold defaultstate="collapsed" desc="Scan File">
        // Scanner fileScanner;
        // try {
        //     fileScanner = new Scanner(new File("hiv-ddi.fasta"));
        // } catch (FileNotFoundException e) {
        //     System.out.println(e.toString());
        //     return;
        // }
        // // </editor-fold>
        //
        // // the below code processes the scanner made above to create an MSA (multiple sequence alignment)
        // // <editor-fold defaultstate="collapsed" desc="Create MSA">
        // ArrayList<String> names = new ArrayList<>();
        // ArrayList<String> alignments = processMSA(fileScanner, names);
        // String consensus = "-----------------------------KIKAL-EICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLD--FRKYTAFTIPS-NNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQY-DDLYVGSDLEIGQHR-KIEELR-HLL-WGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPI-----------------------------------------------------------------------------------------";
        // System.out.println("MSA SIZE: " + alignments.size());
        // System.out.println("CONSENSUS LENGTH: " + consensus.length());
        // ArrayList<Byte[]> bitAlignments = new ArrayList<>();
        // alignments.stream().forEach((s) -> {
        //     System.out.println(s);
        //     bitAlignments.add(proteinToBitString(consensus, s));
        // });
        // </editor-fold>
        // constantGeneration(args);
        // HamiltonianConstants hivUntreated = readConstants("hiv-untreated.hiv");
        // HamiltonianConstants hivTreated = readConstants("hiv-treated.hiv");
        //
        // for (int i = 0; i < bitAlignments.size(); i++) {
        //     System.out.println(names.get(i) + "\t" + evaluateHamiltonian(bitAlignments.get(i), hivUntreated.getH(), hivUntreated.getJ()) + "\t" + evaluateHamiltonian(bitAlignments.get(i), hivTreated.getH(), hivTreated.getJ()));
        // }
        // constantGeneration(new String[]{"hiv-untreated-FINAL.hiv"}, "hiv-untreated_subalign.fasta");
        // constantGeneration(new String[]{"hiv-treated-FINAL.hiv"}, "hiv-rt-treated_subalign.fasta");

        HamiltonianConstants hc_treated = readConstants("hiv-treated-FINAL.hiv");
        HamiltonianConstants hc_untreated = readConstants("hiv-untreated-FINAL.hiv");

        double untreated_partition = partition(410, hc_untreated.getH(), hc_untreated.getJ(), false);
        double treated_partition = partition(410, hc_treated.getH(), hc_treated.getJ(), false);

        System.out.println("I\tTreated\tUntreated");

        for (int i = 0; i < 410; i++) {
            double treated = MMC(i, 410, hc_treated.getH(), hc_treated.getJ(), treated_partition);
            double untreated = MMC(i, 410, hc_untreated.getH(), hc_untreated.getJ(), untreated_partition);
            System.out.println(i+"\t"+treated+"\t"+untreated);
        }
    }

    /**
     * Constant Generation.
     *
     * This function does all legwork for constant generation, running necessary
     * functions and scanning necessary files.
     * @param args Array from main() because args[0] should be file to save constants.
     */
    public static void constantGeneration(String args[], String inputFile) {
        // the below code scans a file for an MSA in FASTA format
        // <editor-fold defaultstate="collapsed" desc="Scan File">
        Scanner fileScanner;
        try {
            fileScanner = new Scanner(new File(inputFile));
        } catch (FileNotFoundException e) {
            System.out.println(e.toString());
            return;
        }
        // </editor-fold>

        // the below code processes the scanner made above to create an MSA (multiple sequence alignment)
        // <editor-fold defaultstate="collapsed" desc="Create MSA">
        ArrayList<String> names = new ArrayList<>();
        ArrayList<String> alignments = processMSA(fileScanner, names);
        String consensus = "-----------------------------KIKAL-EICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLD--FRKYTAFTIPS-NNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQY-DDLYVGSDLEIGQHR-KIEELR-HLL-WGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPI----------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        System.out.println("MSA SIZE: " + alignments.size());
        System.out.println("ALIGNMENT LENGTH: " + alignments.get(0).length());
        System.out.println("CONSENSUS LENGTH: " + consensus.length());
        ArrayList<Byte[]> bitAlignments = new ArrayList<>();
        alignments.stream().forEach((s) -> {
            System.out.println(s);
            bitAlignments.add(proteinToBitString(consensus, s));
        });
        // </editor-fold>

        // the below code creates Hamiltonian Constants h_i and J_ij
        // <editor-fold defaultstate="collapsed" desc="Generate Hamiltonian Constants">
        System.out.println(alignments.get(0).substring(0,1));
        double h[] = new double[consensus.length()];
        double J[][] = new double[consensus.length()][consensus.length()];
        for (int i = 0; i < J.length; i++) {
          for (int j = 0; j < J.length; j++) {
            J[i][j] = 0;
          }
        }
        double singleProbs[] = calculateSingleMutationalProbabilities(bitAlignments, consensus.length());
        System.out.println("Done with single probs");
        double doubleProbs[][] = calculateDoubleMutationalProbabilities(bitAlignments, consensus.length());
        System.out.println("Done with double probs");
        long timer = System.currentTimeMillis();
        // J = estimateJ(singleProbs, doubleProbs, consensus.length());
        HamiltonianConstants hc = boltzmannLearn(singleProbs, doubleProbs, consensus.length(), h, J);
        // </editor-fold>

        // the below code evaluates the Hamiltonian for each input sequence
        //<editor-fold defaultstate="collapsed" desc="Print Hamiltonians">
        double partition = partition(consensus.length(), h, J, false);
        for (int i = 0; i < bitAlignments.size(); i++) {
            System.out.println(names.get(i) + ": " + Math.exp(-1*evaluateHamiltonian(bitAlignments.get(i), hc.getH(), hc.getJ()))/partition);
        }
        //</editor-fold>

        System.out.println("TIME: " + (long) (System.currentTimeMillis() - timer));

        // the below code saves the Hamiltonian constants to a ".hiv" file
        //<editor-fold defaultstate="collapsed" desc="Save Hamiltonian Constants">
        FileOutputStream f_out = null;
        try {
            f_out = new FileOutputStream(args[0] + ".hiv");
        } catch (ArrayIndexOutOfBoundsException ex) {
            try {
                f_out = new FileOutputStream(System.currentTimeMillis() + ".hiv");
            } catch (FileNotFoundException ex1) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex1);
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
        try {
            ObjectOutputStream obj_out = new ObjectOutputStream(f_out);
            obj_out.writeObject(hc);
        } catch (IOException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
        //</editor-fold>
    }

    /**
     * Reading Hamiltonian Constants from a file.
     *
     * @param path  Path of file to read (should be in code-saved .hiv format)
     * @return      HamiltonianConstants object if read or null if failed
     */
    public static HamiltonianConstants readConstants (String path) {
        FileInputStream f_in = null;
        HamiltonianConstants toReturn = null;
        try {
            f_in = new FileInputStream(path);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
        try {
            ObjectInputStream obj_in = new ObjectInputStream(f_in);
            toReturn = (HamiltonianConstants) obj_in.readObject();
        } catch (IOException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
        return toReturn;
    }

    /**
     * Process MSA (Multiple Sequence Alignment).
     *
     * This processes an MSA from an input stream in FASTA format. First, it
     * essentially just converts the stream MSA into a list of sequences
     *
     * @param scan Scanner object defining the input stream to use (usually from
     * a file)
     * @param names Empty ArrayList to which the function will add the names of
     * each sequence (as directed by FASTA)
     * @return ArrayList<String> containing each sequence from the MSA in order
     */
    public static ArrayList<String> processMSA(Scanner scan, ArrayList<String> names) {
        ArrayList<String> out = new ArrayList<>();
        String curString = "";
        while (scan.hasNextLine()) {
            String a = scan.nextLine();
            if (a.equals("")) {
                //skip this line
            } else if (a.substring(0, 1).equals(">")) {
                out.add(curString);
                curString = "";
                System.out.println(a);
                names.add(a.substring(1, a.length()));
            } else {
                curString += removeWhitespace(a);
            }
        }
        out.add(curString);
        out.remove(0);
        return out;
    }

    /**
     * Remove Whitespace function.
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
     * Protein to Bit String Function.
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
     * Boolean to Byte Function.
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
     * Evaluate Hamiltonian Function.
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
    public static class HamiltonianConstants implements Serializable {

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
     * Boltzmann Learning function.
     *
     * Generates constant vectors h_i and J_ij for an MSA. The process is to run
     * through the Hamiltonian and in each step, edit h and J such that they
     * predict closer to the observed mutational probabilities
     *
     * @param singleProbs Observed mutational probabilities at each locus
     * @param doubleProbs Observed double mutational probabilities at each pair
     * of loci
     * @param length length of protein
     * @param h Initial value for h, the constant vector
     * @param J Initial value for J, the constant vector
     * @return HamiltonianConstants object containing corrected constant values
     */
    public static HamiltonianConstants boltzmannLearn(double[] singleProbs, double[][] doubleProbs, int length, double[] h, double[][] J) {
        int count = 0;
        double[] delta_h = new double[h.length];
        double[][] delta_J = new double[J.length][J[0].length];
        System.out.println("Entering the Boltzmann loop");
        while (count++ < BOLTZ) {
            double p = partition(length, h, J, false);
            System.out.println("Partition: "+p);
            for (int i = 0; i < length; i++) {
                double adjustedSingleProbs = singleProbs[i];
                double adjustedMMC = MMC(i, length, h, J, p);
                // double adjustedMMC = thermalAverage(i, length, h, J, p);
                delta_h[i] = ETA * (adjustedSingleProbs - adjustedMMC);
                // System.out.println(i+"\tSingle Prob: "+adjustedSingleProbs+"\tMMC: "+adjustedMMC+"\tdelta h: "+delta_h[i]+"\th: "+(double)(delta_h[i]+h[i]));
            }
            for (int i = 0; i < h.length; i++) {
                h[i] += delta_h[i];
            }
            System.out.println("Learning step: " + count);
        }
        return new HamiltonianConstants(h, J);
    }

    /**
     * Thermal Averaging function (one position). Deprecated.
     *
     * Finds the average over the Ising distribution of a position. This finds,
     * essentially, a weighted sum of all possible sequences in the Hamiltonian
     * in which the the position specified is a mutation. This is extremely
     * computationally intensive, so we instead take a Monte Carlo sample of the
     * possible sequences (20,000 trials) and take that average instead. This
     * sacrifices accuracy but gains efficiency.
     *
     * Thermal averaging of one position is not currently in use, replaced by
     * MMC (Metropolis Monte Carlo). That function does essentially the same
     * thing, except that it should be computationally much faster than the
     * thermal average function with equal or greater accuracy.
     *
     * @param position position of amino acid that is mutated
     * @param length length of bit string protein (number of amino acids)
     * @param h Hamiltonian constant h
     * @param J Hamiltonian constant J
     * @return The thermal average over the Ising distribution
     */
    public static double thermalAverage(int position, int length, double[] h, double[][] J, double partition) {
        double sum = 0;
        Random r = new Random();
        for (int i = 0; i < MC; i++) {
            Byte[] rand = new Byte[length];
            for (int j = 0; j < length; j++) {
                rand[j] = r.nextBoolean() ? new Byte((byte) 1) : new Byte((byte) 0);
            }
            rand[position] = 1;
            sum += Math.exp(-evaluateHamiltonian(rand, h, J)) / partition;
        }

        return sum / MC;
    }

    /**
     * Thermal Averaging function (two positions). Deprecated.
     *
     * Finds the average over the Ising distribution of a position. This finds,
     * essentially, a weighted sum of all possible sequences in the Hamiltonian
     * in which the the position specified is a mutation. This is extremely
     * computationally intensive, so we instead take a Monte Carlo sample of the
     * possible sequences (20,000 trials) and take that average instead. This
     * sacrifices accuracy but gains efficiency.
     *
     * Thermal averaging of two positions is not currently in use because we can
     * estimate J with reasonable accuracy in a different way.
     *
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
     * Calculates mutational probabilities by averaging them over the MSA.
     *
     * @param msa Bytestring MSA
     * @param length length of each
     * @return Probability of each mutation
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
     * Calculates double mutational probabilities by averaging them over the
     * MSA.
     *
     * @param msa Bytestring MSA
     * @param length length of each
     * @return Probability of each double mutation
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

    /**
     * Estimates J using the Independent-Pair Approximation.
     *
     * This is an estimation of the J constants (not great accuracy but decent
     * for my purposes) from Roudi et. al. 2009. Added this because using
     * Boltzmann learning to calculate J is extremely computationally intensive
     * and this has (apparently) a decent level of accuracy.
     *
     * @param singleMutations Single mutational probabilities vector
     * @param doubleMutations Double mutational probabilities matrix
     * @param length Length of protein
     * @return Matrix of all J_ij where j > i
     */
    public static double[][] estimateJ(double[] singleMutations, double[][] doubleMutations, int length) {
        double[][] J = new double[length][length];
        for (int i = 0; i < length; i++) {
            for (int j = i + 1; j < length; j++) {
                double sub = doubleMutations[i][j] / ((1 + singleMutations[i]) * (1 + singleMutations[j]));
                J[i][j] = .25 * Math.log(1 + sub);
            }
        }
        System.out.println("J estimated");
        return J;
    }

    /**
     * MMC (Metropolis Monte Carlo) Function.
     *
     * Follows basic sampling procedures found in the Supplementary Methods of
     * Ferguson et. al. This is essentially the same as the thermal average
     * function, but instead of generating random sequences each time, it
     * "steps" through all possible states (as in a Markov Chain) and records
     * the CHANGE in energy each time. That way, evaluateHamiltonian() is only
     * called once. Then we just record every MC_RATE state and average.
     *
     * This used to be done differently, rejecting based on whether or not it
     * increased the energy. I used this algorithm from Beichl and Sullivan,
     * 2000. However, this was a mistake to use because it is useful for finding
     * a local minimum of maximum but not for taking a straight average.
     *
     * @param position position of amino acid that is mutated
     * @param length length of protein
     * @param h current h vector
     * @param J current J matrix
     * @return double representing the thermal average over the Ising
     * distribution
     */
    public static double MMC(int position, int length, double[] h, double[][] J, double partition) {
        Random rand = new Random();
        Byte[] sigma = new Byte[length];
        for (int i = 0; i < length; i++) {
            sigma[i] = 0;
        }
        sigma[position] = 1; //position is the place that must be mutated for a proper average
        double sum = 0.0;
        int count = 1;
        double E = evaluateHamiltonian(sigma, h, J);
        sum += E;
        while (count < MC_WALK) {
            int toChange = rand.nextInt(length);

            while (toChange == position) {
                toChange = rand.nextInt(length); //ensures that we don't change position
            }
            if (sigma[toChange] == 1) {
                sigma[toChange] = 0;
            } else {
                sigma[toChange] = 1;
            }

            E += calculateDeltaE(toChange, length, h, J, sigma);
            if (count % MC_RATE == 0) {
                sum += Math.exp(-E)/partition;
            }
            count++;
        }
        // System.out.println(sum / ((int) (count / MC_RATE)));
        return sum / ((int) (count / MC_RATE));
    }

    /**
     * Helper method for MMC to calculate the deltaE (change in energy).
     *
     * This essentially allows a much faster sampling because we don't have to
     * go through the entire evaluateHamiltonian() step in order to find the new
     * energy of the system after a single flip in the Monte Carlo MMC().
     *
     * @param change amino acid that was changed
     * @param length length of protein
     * @param h h vector of Ising
     * @param J J matrix of Ising
     * @param sigma Byte Array of protein containing elements from the set {0,1}
     * @return delta E of the new state
     */
    private static double calculateDeltaE(int change, int length, double[] h, double[][] J, Byte[] sigma) {
        double delta_E = 0.0;
        if (change == 0) { //this means that it has been changed to a 0 from a 1
            delta_E -= h[change]; //sigma[change]*h[change] used to be in E, but it is no longer
            for (int j = change + 1; j < length; j++) { //sigma[j]*sigma[change]*J[change][j] is no longer in E
                delta_E -= J[change][j] * sigma[j];
            }
            for (int i = 0; i < change; i++) { //sigma[i]*sigma[change]*J[i][change] is no longer in E
                delta_E -= J[i][change] * sigma[i];
            }
        } else { //this means that is has been changed to a 1 from a 0
            delta_E += h[change]; //we now include h[change]*sigma[change]
            for (int j = change + 1; j < length; j++) { //sigma[j]*sigma[change]*J[change][j] is added
                delta_E += J[change][j] * sigma[j];
            }
            for (int i = 0; i < change; i++) { //sigma[i]*sigma[change]*J[i][change] is added
                delta_E += J[i][change] * sigma[i];
            }
        }
        // System.out.println(delta_E);
        return delta_E;
    }

    /**
     * Partition Function.
     *
     * Finds the partition function as an average over the Ising distribution times
     * the number of amino acids in the distribution. Takes MC_PARTITION totally
     * random samples (no MC-walk for now).
     */
    public static Double partition(int length, double[] h, double[][] J, boolean verbose) {
        double sum = 0;
        Random r = new Random();
        double exponent = -1;
        for (int i = 0; i < MC_PARTITION; i++) {
            Byte[] rand = new Byte[length];
            for (int j = 0; j < length; j++) {
                rand[j] = r.nextBoolean() ? new Byte((byte) 1) : new Byte((byte) 0);
            }
            double x = evaluateHamiltonian(rand,h,J);
            // if (exponent == -1) exponent = Math.getExponent(x);
            // x /= Math.exp(exponent);
            if (verbose) System.out.println("Partition Hamiltonian: "+x);
            sum += Math.exp(-1*x);
            if(verbose) System.out.println(i+" partition, sum is "+sum);
            if (i%10000==0) {
                // for (int j = 0; j < length; j++) {
                    // System.out.println(j+" H VALUE "+h[j]);
                // }
                System.out.println(i+" out of "+MC_PARTITION+"; sum is "+sum);
            }
        }
        // sum *= Math.exp(Math.exp(exponent));
        System.out.println("Average value is "+(double)(sum / MC_PARTITION));
        return new Double(sum / MC_PARTITION * Math.pow(2,length));
    }
}

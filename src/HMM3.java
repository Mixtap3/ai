import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Locale;
import java.util.Scanner;

/**
 * Created by Leonard on 2017-09-11.
 */
public class HMM3 {

    Scanner sc;
    double[][] A;
    double[][] B;
    double[][] pi;
    int[] obs;
    double[][] alpha;
    double[][] beta;
    double[][] digamma;
    double[][] gamma;
    int timeSteps;
    int stateLength;
    double[] c;
    double logProb;
    double oldLogProb;
    int iterLimit = 500;
    int iters = 0;

    HMM3(Scanner sc) {
        this.sc = sc;
        this.createMatrix();
    }

    void createMatrix() {

        int arow = sc.nextInt();
        int acol = sc.nextInt();
        this.A = new double[arow][acol];
        for (int i = 0; i < arow; i++) {
            for (int j = 0; j < acol; j++) {
                this.A[i][j] = sc.nextDouble();
            }
        }

        int brow = sc.nextInt();
        int bcol = sc.nextInt();
        this.B = new double[brow][bcol];
        for (int i = 0; i < brow; i++) {
            for (int j = 0; j < bcol; j++) {
                this.B[i][j] = sc.nextDouble();
            }
        }

        int pirow = sc.nextInt();
        int picol = sc.nextInt();
        this.pi = new double[pirow][picol];
        for (int i = 0; i < pirow; i++) {
            for (int j = 0; j < picol; j++) {
                this.pi[i][j] = sc.nextDouble();
            }
        }

        stateLength = pi[0].length;

        timeSteps = sc.nextInt();
        obs = new int[timeSteps];

        for (int i = 0; i < obs.length; i++) {
            obs[i] = sc.nextInt();
        }

        c = new double[timeSteps];

        this.baum();
    }

    void forward(){
        alpha = new double[stateLength][timeSteps];
        this.initAlpha();

        for (int t = 1; t < timeSteps; t++) {
            c[t] = 0.0;
            double[][] alphaVec = new double[1][stateLength];
            for (int i = 0; i < stateLength; i++) {
                alphaVec[0][i] = alpha[i][t-1];
            }
            double[][] newAlpha = multiplicar(alphaVec, A);
            for (int i = 0; i < stateLength; i++) {
                alpha[i][t] = newAlpha[0][i] * B[i][obs[t]];
                c[t] += alpha[i][t];
            }

            //scale
            c[t] = 1/c[t];
            for (int i = 0; i < stateLength; i++) {
                alpha[i][t] = alpha[i][t] * c[t];
            }
        }
    }

    void initAlpha(){
        c[0] = 0.0;
        for (int i = 0; i < stateLength; i++) {
            alpha[i][0] = pi[0][i] * B[i][obs[0]];
            c[0] = c[0] + alpha[i][0];
        }

        //scale
        c[0] = 1/c[0];
        for (int i = 0; i < stateLength; i++) {
            alpha[i][0] = c[0]*alpha[i][0];
        }
    }

    void backward(){
        beta = new double[stateLength][timeSteps];
        for (int i = 0; i < stateLength; i++) {
            beta[i][timeSteps-1] = c[timeSteps-1];
        }

        for (int t = timeSteps - 2; t >= 0; t--) {
            for (int i = 0; i < stateLength; i++) {
                beta[i][t] = 0.0;
                for (int j = 0; j < stateLength; j++) {
                    beta[i][t] += A[i][j] * B[j][obs[t+1]] * beta[j][t+1];
                }
                beta[i][t] = beta[i][t] * c[t];
            }
        }
    }

    void createDigamma(){
        digamma = new double[stateLength*stateLength][timeSteps - 1];
        gamma = new double[stateLength][timeSteps];


        for (int t = 0; t < timeSteps - 1; t++) {
            double denom = 0.0;
            for (int i = 0; i < stateLength; i++) {
                for (int j = 0; j < stateLength; j++) {
                    denom += alpha[i][t]*A[i][j]*B[j][obs[t+1]]*beta[j][t+1];
                }
            }

            for (int i = 0; i < stateLength; i++) {
                gamma[i][t] = 0.0;
                for (int j = 0; j < stateLength; j++) {
                    digamma[stateLength*i + j][t] = alpha[i][t] * A[i][j] * B[j][obs[t+1]] * beta[j][t+1] / denom;
                    gamma[i][t] += digamma[stateLength*i + j][t];
                }
            }
        }

        double denom = 0.0;
        for (int i = 0; i < stateLength; i++) {
            denom += alpha[i][timeSteps-1];
        }
        for (int i = 0; i < stateLength; i++) {
            gamma[i][timeSteps - 1] = alpha[i][timeSteps-1] / denom;
        }


    }

    void createGamma(){
        gamma = new double[stateLength][timeSteps - 1];
        for (int t = 0; t < timeSteps - 1; t++) {
            for (int i = 0; i < stateLength; i++) {
                gamma[i][t] = 0.0;
                for (int j = 0; j < stateLength; j++) {
                    gamma[i][t] += digamma[4*i + j][t];
                }
            }
        }
    }

    void transitionEstimates(){
        for (int i = 0; i < stateLength; i++) {
            for (int j = 0; j < stateLength; j++) {
                double digamSum = 0.0;
                double gamSum = 0.0;
                for (int t = 0; t < timeSteps - 1; t++) {
                    digamSum += digamma[stateLength*i + j][t];
                    gamSum += gamma[i][t];
                }
                A[i][j] = digamSum / gamSum;
            }
        }
    }

    void emissionEstimates(){

        for (int i = 0; i < stateLength; i++) {
            for (int k = 0; k < B[0].length; k++) {
                double nom = 0.0;
                double denom = 0.0;
                for (int t = 0; t < timeSteps; t++) {
                    denom += gamma[i][t];
                    if (obs[t] == k) {
                        nom += gamma[i][t];
                    }
                }
                B[i][k] = nom / denom;
            }
        }
    }

    void newInitialState(){
        for (int i = 0; i < stateLength; i++) {
            pi[0][i] = gamma[i][0];
        }
    }

    void baum(){
        oldLogProb = Double.NEGATIVE_INFINITY;
        while (true) {

            iters++;

            if(B[0].length == 1){
                throw new IllegalArgumentException("hej");
            }

            this.forward();
            this.backward();
            this.createDigamma();
            //this.createGamma();
            this.transitionEstimates();
            this.emissionEstimates();
            this.newInitialState();

            this.computeLogProb();

            if(!(logProb > oldLogProb)){
                this.writeResult();
                break;
            }
            oldLogProb = logProb;
        }
    }

    void computeLogProb(){
        logProb = 0.0;
        for (int i = 0; i < timeSteps; i++) {
            logProb += Math.log(c[i]);
        }
        logProb = -logProb;
    }

    void writeResult(){
        StringBuilder sb = new StringBuilder();
        sb.append(A.length + " " + A[0].length);
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                double r = (double)Math.round(A[i][j] * 1000000d) / 1000000d;
                sb.append(" " + r);
            }
        }
        System.out.println(sb.toString());

        StringBuilder sb2 = new StringBuilder();
        sb2.append(B.length + " " + B[0].length);
        for (int i = 0; i < B.length; i++) {
            for (int j = 0; j < B[0].length; j++) {
                double r = (double)Math.round(B[i][j] * 1000000d) / 1000000d;
                sb2.append(" " + r);
            }
        }
        System.out.println(sb2.toString());
    }

    public static double[][] multiplicar(double[][] A, double[][] B) {
        int aRows = A.length;
        int aColumns = A[0].length;
        int bRows = B.length;
        int bColumns = B[0].length;
        if (aColumns != bRows) {
            throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
        }

        double[][] C = new double[aRows][bColumns];
        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bColumns; j++) {
                C[i][j] = 0.00000;
            }
        }

        for (int i = 0; i < aRows; i++) { // aRow
            for (int j = 0; j < bColumns; j++) { // bColumn
                for (int k = 0; k < aColumns; k++) { // aColumn
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }


    public static void main(String[] args) throws IOException {
//        File file = new File("testCase3.txt");
//        System.setIn(new FileInputStream(file));

        Scanner sc = new Scanner(System.in).useLocale(Locale.US);
        HMM3 m = new HMM3(sc);

    }
}

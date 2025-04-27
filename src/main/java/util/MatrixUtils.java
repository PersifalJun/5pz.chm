package util;

public class MatrixUtils {
    public static double[][] generateMatrix(int N) {
        double[][] A = new double[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = 1.0 / (1 + 0.8 * (i + 1) + 5 * (j + 1));
            }
        }
        return A;
    }

    public static double[] multiply(double[][] A, double[] x) {
        int N = x.length;
        double[] result = new double[N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                result[i] += A[i][j] * x[j];
            }
        }
        return result;
    }

    public static double relativeError(double[] x, double[] xStar) {
        double num = 0, den = 0;
        for (int i = 0; i < x.length; i++) {
            num += Math.pow(x[i] - xStar[i], 2);
            den += Math.pow(xStar[i], 2);
        }
        return Math.sqrt(num / den);
    }

    public static boolean checkApproxEquals(double[][] A, double[][] B, double eps) {
        int N = A.length;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (Math.abs(A[i][j] - B[i][j]) > eps) {
                    return false;
                }
            }
        }
        return true;
    }

    public static double[][] multiply(double[][] A, double[][] B) {
        int N = A.length;
        double[][] result = new double[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return result;
    }
}

package org.entity;

import util.MatrixUtils;
import java.util.Arrays;

public class Main {

    public static void main(String[] args) {
        int[] sizes = {5, 10, 20}; // Размерности N
        int trials = 100;          // Количество прогонов для усреднения времени

        for (int N : sizes) {
            System.out.println("Размерность N = " + N);

            double[] xStar = new double[N];
            Arrays.fill(xStar, 1.0); // Вектор точного решения (x*)

            double[][] A = MatrixUtils.generateMatrix(N);
            double[] f = MatrixUtils.multiply(A, xStar);

            // LU-разложение
            double[][] L = new double[N][N];
            double[][] U = new double[N][N];
            double[] xLU = LUDecomposition.solve(A, f, L, U);
            double errLU = MatrixUtils.relativeError(xLU, xStar);
            boolean isLUCorrect = MatrixUtils.checkApproxEquals(MatrixUtils.multiply(L, U), A, 1e-10);
            long avgTimeLU = measureTimeLU(N, trials);

            //  QR-разложение (Гивенс)
            double[][] Q = new double[N][N];
            double[][] R = new double[N][N];
            double[] xQR = GivensQR.solve(A, f, Q, R);
            double errQR = MatrixUtils.relativeError(xQR, xStar);
            boolean isQRCorrect = MatrixUtils.checkApproxEquals(MatrixUtils.multiply(Q, R), A, 1e-10);
            long avgTimeQR = measureTimeQR(N, trials);

            //  SVD-разложение
            double[] xSVD = null;
            double[][] Umat = null, Sigmamat = null, Vmat = null;
            long totalTimeSVD = 0;
            double condSVD = 0;

            for (int t = 0; t < trials; t++) {
                double[][] Atrial = MatrixUtils.generateMatrix(N);
                double[] ftrial = MatrixUtils.multiply(Atrial, xStar);

                long start = System.nanoTime();
                SVD.SVDResult result = SVD.decompose(Atrial);
                Umat = result.U;
                Sigmamat = result.Sigma;
                Vmat = result.V;
                xSVD = solveSVD(Umat, Sigmamat, Vmat, ftrial);
                long end = System.nanoTime();

                totalTimeSVD += (end - start);
            }

            condSVD = conditionNumber(Sigmamat);
            long avgTimeSVD = totalTimeSVD / trials / 1000;
            double errSVD = MatrixUtils.relativeError(xSVD, xStar);

            //  Вывод результатов
            System.out.println("LU-разложение:");
            System.out.println("Среднее время (мкс): " + avgTimeLU);
            System.out.println("Погрешность δ = " + errLU);
            System.out.println("LU ≈ A: " + isLUCorrect);
            System.out.println();

            System.out.println("QR-разложение (Гивенс):");
            System.out.println("Среднее время (мкс): " + avgTimeQR);
            System.out.println("Погрешность δ = " + errQR);
            System.out.println("QR ≈ A: " + isQRCorrect);
            System.out.println();

            System.out.println("SVD-разложение:");
            System.out.println("Среднее время (мкс): " + avgTimeSVD);
            System.out.println("Погрешность δ = " + errSVD);
            System.out.println("Число обусловленности Cond(A) = " + condSVD);
            System.out.println();
        }
    }

    // Функции для замера времени

    public static long measureTimeLU(int N, int trials) {
        double[] x = new double[N];
        Arrays.fill(x, 1.0);
        long totalTime = 0;

        for (int t = 0; t < trials; t++) {
            double[][] A = MatrixUtils.generateMatrix(N);
            double[] f = MatrixUtils.multiply(A, x);
            double[][] L = new double[N][N];
            double[][] U = new double[N][N];

            long start = System.nanoTime();
            LUDecomposition.solve(A, f, L, U);
            long end = System.nanoTime();

            totalTime += (end - start);
        }
        return totalTime / trials / 1000;
    }

    public static long measureTimeQR(int N, int trials) {
        double[] x = new double[N];
        Arrays.fill(x, 1.0);
        long totalTime = 0;

        for (int t = 0; t < trials; t++) {
            double[][] A = MatrixUtils.generateMatrix(N);
            double[] f = MatrixUtils.multiply(A, x);
            double[][] Q = new double[N][N];
            double[][] R = new double[N][N];

            long start = System.nanoTime();
            GivensQR.solve(A, f, Q, R);
            long end = System.nanoTime();

            totalTime += (end - start);
        }
        return totalTime / trials / 1000;
    }

    // Функции для SVD-псевдорешения

    public static double[] solveSVD(double[][] U, double[][] Sigma, double[][] V, double[] f) {
        int N = f.length;
        double[] y = new double[N];
        double[] temp = new double[N];
        double[] x = new double[N];

        // Находим максимальное сингулярное число
        double maxSigma = Sigma[0][0];
        double threshold = maxSigma * 1e-12; // Относительный порог

        // y = U^T * f
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                y[i] += U[j][i] * f[j];
            }
        }

        // temp = Sigma^-1 * y (но только если sigma[i] не слишком маленькая)
        for (int i = 0; i < N; i++) {
            if (Math.abs(Sigma[i][i]) > threshold) {
                temp[i] = y[i] / Sigma[i][i];
            } else {
                temp[i] = 0;
            }
        }

        // x = V * temp
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                x[i] += V[i][j] * temp[j];
            }
        }

        return x;
    }

    //  Функция для расчета числа обусловленности

    public static double conditionNumber(double[][] Sigma) {
        double sigmaMax = Sigma[0][0];
        double sigmaMin = Double.MAX_VALUE;
        int N = Sigma.length;

        for (int i = 0; i < N; i++) {
            if (Math.abs(Sigma[i][i]) > 1e-15) {
                sigmaMin = Math.min(sigmaMin, Sigma[i][i]);
            }
        }

        return sigmaMax / sigmaMin;
    }
}

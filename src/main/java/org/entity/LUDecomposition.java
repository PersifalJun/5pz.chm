package org.entity;

public class LUDecomposition {
    static final double EPS = 1e-20; // малое число для защиты от деления на 0 (если потребуется)

    public static double[] solve(double[][] A, double[] b, double[][] L_out, double[][] U_out) {
        int N = A.length;

        // Инициализируем матрицы L (нижняя треугольная) и U (верхняя треугольная)
        double[][] L = new double[N][N];
        double[][] U = new double[N][N];

        //  LU-разложение A = L * U
        for (int i = 0; i < N; i++) {
            // Заполнение U (верхняя треугольная матрица)
            for (int j = i; j < N; j++) {
                U[i][j] = A[i][j];
                for (int k = 0; k < i; k++) {
                    // Вычитаем скалярное произведение соответствующих элементов L и U
                    U[i][j] -= L[i][k] * U[k][j];
                }
            }

            // Заполнение L (нижняя треугольная матрица)
            for (int j = i; j < N; j++) {
                if (i == j) {
                    L[i][i] = 1; // диагональные элементы L = 1 по определению
                } else {
                    L[j][i] = A[j][i];
                    for (int k = 0; k < i; k++) {
                        // Вычитаем накопленную сумму из предыдущих столбцов
                        L[j][i] -= L[j][k] * U[k][i];
                    }
                    L[j][i] /= U[i][i]; // делим на диагональный элемент U
                }
            }
        }

        //  копируем L и U во внешние массивы (если переданы)
        if (L_out != null && U_out != null) {
            for (int i = 0; i < N; i++) {
                System.arraycopy(L[i], 0, L_out[i], 0, N);
                System.arraycopy(U[i], 0, U_out[i], 0, N);
            }
        }

        // прямой ход — решаем L * y = b
        double[] y = new double[N];
        for (int i = 0; i < N; i++) {
            y[i] = b[i];
            for (int j = 0; j < i; j++) {
                y[i] -= L[i][j] * y[j]; // подставляем уже найденные y[j]
            }
        }

        // обратный ход — решаем U * x = y
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            x[i] = y[i];
            for (int j = i + 1; j < N; j++) {
                x[i] -= U[i][j] * x[j]; // подставляем уже найденные x[j]
            }
            x[i] /= U[i][i]; // делим на диагональный элемент
        }

        // Возвращаем решение x
        return x;
    }
}
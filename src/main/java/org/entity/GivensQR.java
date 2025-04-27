package org.entity;

public class GivensQR {
    static final double EPS = 1e-20;

    public static double[] solve(double[][] A, double[] b, double[][] Q_out, double[][] R_out) {
        int N = A.length;

        // Инициализация матрицы R (копия A) и Q (единичная матрица)
        double[][] R = new double[N][N];
        double[][] Q = new double[N][N];
        for (int i = 0; i < N; i++) {
            Q[i][i] = 1; // Единичная матрица Q
            System.arraycopy(A[i], 0, R[i], 0, N); // Копируем A в R
        }

        // Проход по всем столбцам (QR-разложение с помощью поворотов Гивенса)
        for (int j = 0; j < N - 1; j++) {
            for (int i = j + 1; i < N; i++) {
                // Если элемент под диагональю не ноль, надо его занулить
                if (Math.abs(R[i][j]) > EPS) {
                    // Вычисляем гипотенузу (длину вектора из двух элементов) для нормализации
                    double r = Math.hypot(R[j][j], R[i][j]);

                    // Вычисляем cos(θ) и sin(θ) для матрицы поворота Гивенса
                    double c = R[j][j] / r;
                    double s = R[i][j] / r;

                    // Применяем поворот Гивенса к R (Gt * R)
                    for (int k = j; k < N; k++) {
                        double Rjk = R[j][k];
                        double Rik = R[i][k];
                        R[j][k] = c * Rjk + s * Rik;
                        R[i][k] = -s * Rjk + c * Rik;
                    }

                    // Применяем тот же поворот к Q (Q = Q * G)
                    for (int k = 0; k < N; k++) {
                        double Qkj = Q[k][j];
                        double Qki = Q[k][i];
                        Q[k][j] = c * Qkj + s * Qki;
                        Q[k][i] = -s * Qkj + c * Qki;
                    }
                }
            }
        }

        // Если переданы внешние ссылки на Q и R — скопировать результат
        if (Q_out != null && R_out != null) {
            for (int i = 0; i < N; i++) {
                System.arraycopy(Q[i], 0, Q_out[i], 0, N);
                System.arraycopy(R[i], 0, R_out[i], 0, N);
            }
        }

        // Умножаем транспонированную Q на f (b), получаем Q^T * b
        double[] Qtb = new double[N];
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < N; k++) {
                Qtb[i] += Q[k][i] * b[k]; // Q^T * b
            }
        }

        // Обратный ход: решаем R * x = Q^T * b
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            x[i] = Qtb[i];
            for (int j = i + 1; j < N; j++) {
                x[i] -= R[i][j] * x[j];
            }
            x[i] /= R[i][i]; // делим на диагональный элемент
        }

        return x; // возвращаем найденное решение x
    }
}

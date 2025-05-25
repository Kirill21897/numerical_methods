using System;
using ConsoleChMethod;

class EquationSolver
{
    // Делегат для функции, возвращающей производную (например, dx/dt)
    public delegate Vector DerivativeFunction(double t, Vector x);

    // Пример: гармонические колебания (x'' = -x)
    public static Vector Oscillator(double t, Vector x)
    {
        Vector result = new Vector(2);
        result[0] = x[1];     // dx/dt = v (скорость)
        result[1] = -x[0];    // dv/dt = -x (ускорение)
        return result;
    }

    // Метод Эйлера — простейший метод численного интегрирования
    public static Matrix SolveEuler(double t0, double tEnd, Vector initialState, int steps, DerivativeFunction f)
    {
        int dim = initialState.Size;  // Размерность системы уравнений
        Matrix solution = new Matrix(dim + 1, steps + 1);  // Решение: время + координаты/скорости
        double h = (tEnd - t0) / steps;  // Шаг по времени

        Vector currentRow = new Vector(dim + 1);
        currentRow[0] = t0;
        for (int i = 0; i < dim; i++) currentRow[i + 1] = initialState[i];
        solution.SetColumn(0, currentRow);  // Сохраняем начальное состояние

        Vector state = initialState.Copy();  // Текущее состояние системы
        double time = t0;

        // Основной цикл по шагам
        for (int step = 1; step <= steps; step++)
        {
            Vector derivative = f(time, state);  // Вычисляем производные
            state = state + derivative * h;      // Обновляем состояние
            time += h;

            currentRow[0] = time;
            for (int i = 0; i < dim; i++) currentRow[i + 1] = state[i];

            solution.SetColumn(step, currentRow);  // Сохраняем результат
        }

        return solution;
    }

    // Метод Рунге–Кутты второго порядка
    public static Matrix SolveRungeKutta2(double t0, double tEnd, Vector initialState, int steps, DerivativeFunction f)
    {
        int dim = initialState.Size;
        Matrix solution = new Matrix(dim + 1, steps + 1);
        double h = (tEnd - t0) / steps;

        Vector currentRow = new Vector(dim + 1);
        currentRow[0] = t0;
        for (int i = 0; i < dim; i++) currentRow[i + 1] = initialState[i];
        solution.SetColumn(0, currentRow);

        Vector state = initialState.Copy();
        double time = t0;

        for (int step = 1; step <= steps; step++)
        {
            Vector k1 = f(time, state);                     // Первый шаг
            Vector k2 = f(time + h, state + k1 * h);        // Второй шаг

            state = state + (k1 + k2) * (h / 2.0);         // Обновление состояния
            time += h;

            currentRow[0] = time;
            for (int i = 0; i < dim; i++) currentRow[i + 1] = state[i];

            solution.SetColumn(step, currentRow);
        }

        return solution;
    }

    // Метод Рунге–Кутты четвёртого порядка — точный и популярный
    public static Matrix SolveRungeKutta4(double t0, double tEnd, Vector initialState, int steps, DerivativeFunction f)
    {
        int dim = initialState.Size;
        Matrix solution = new Matrix(dim + 1, steps + 1);
        double h = (tEnd - t0) / steps;

        Vector currentRow = new Vector(dim + 1);
        currentRow[0] = t0;
        for (int i = 0; i < dim; i++) currentRow[i + 1] = initialState[i];
        solution.SetColumn(0, currentRow);

        Vector state = initialState.Copy();
        double time = t0;

        for (int step = 1; step <= steps; step++)
        {
            Vector k1 = f(time, state);
            Vector k2 = f(time + h / 2, state + k1 * h / 2);  // 2-й шаг
            Vector k3 = f(time + h / 2, state + k2 * h / 2);  // 3-й шаг
            Vector k4 = f(time + h, state + k3 * h);           // 4-й шаг

            state = state + (k1 + 2 * k2 + 2 * k3 + k4) * (h / 6.0);  // Обновление
            time += h;

            currentRow[0] = time;
            for (int i = 0; i < dim; i++) currentRow[i + 1] = state[i];

            solution.SetColumn(step, currentRow);
        }

        return solution;
    }

    // Метод Адамса (четвертого порядка) — многошаговый метод
    public static Matrix SolveAdams(double t0, double tEnd, Vector initialState, int steps, DerivativeFunction f)
    {
        int dim = initialState.Size;
        Matrix solution = new Matrix(dim + 1, steps + 1);
        double h = (tEnd - t0) / steps;

        Vector currentRow = new Vector(dim + 1);
        currentRow[0] = t0;
        for (int i = 0; i < dim; i++) currentRow[i + 1] = initialState[i];
        solution.SetColumn(0, currentRow);

        Vector[] states = new Vector[steps + 1];   // Хранение состояний
        double[] times = new double[steps + 1];    // Хранение времён
        Vector[] derivatives = new Vector[steps + 1];  // Производные

        times[0] = t0;
        states[0] = initialState.Copy();
        derivatives[0] = f(times[0], states[0]);

        // Используем метод Рунге–Кутты для первых 4 шагов
        for (int i = 1; i <= Math.Min(3, steps); i++)
        {
            times[i] = times[i - 1] + h;
            derivatives[i] = f(times[i], states[i - 1]);
            states[i] = states[i - 1] + derivatives[i] * h;

            currentRow[0] = times[i];
            for (int j = 0; j < dim; j++) currentRow[j + 1] = states[i][j];
            solution.SetColumn(i, currentRow);
        }

        // Основной цикл метода Адамса
        for (int i = 3; i < steps; i++)
        {
            times[i + 1] = times[i] + h;
            // Формула Адамса четвёртого порядка
            states[i + 1] = states[i] +
                (55 * derivatives[i] - 59 * derivatives[i - 1] +
                 37 * derivatives[i - 2] - 9 * derivatives[i - 3]) * h / 24.0;
            derivatives[i + 1] = f(times[i + 1], states[i + 1]);

            currentRow[0] = times[i + 1];
            for (int j = 0; j < dim; j++) currentRow[j + 1] = states[i + 1][j];
            solution.SetColumn(i + 1, currentRow);
        }

        return solution;
    }

    // Модифицированный метод Рунге–Кутты второго порядка
    // Использует предиктор-корректор
    public static Matrix SolveRungeKutta2_Modified(double t0, double tEnd, Vector initialState, int steps, DerivativeFunction f)
    {
        int dim = initialState.Size;
        Matrix solution = new Matrix(dim + 1, steps + 1);
        double h = (tEnd - t0) / steps;

        Vector currentRow = new Vector(dim + 1);
        currentRow[0] = t0;
        for (int i = 0; i < dim; i++) currentRow[i + 1] = initialState[i];
        solution.SetColumn(0, currentRow);

        Vector state = initialState.Copy();
        double time = t0;

        for (int step = 1; step <= steps; step++)
        {
            // Предиктор
            Vector k1 = f(time, state);
            Vector yPredicted = state + k1 * h;

            // Корректор
            Vector k2 = f(time + h, yPredicted);
            state = state + (k1 + k2) * (h / 2.0);
            time += h;

            currentRow[0] = time;
            for (int i = 0; i < dim; i++) currentRow[i + 1] = state[i];

            solution.SetColumn(step, currentRow);
        }

        return solution;
    }
}
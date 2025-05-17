using System;
namespace ConsoleChMethod;

class EquationSolver
{
    // Делегат для функции, возвращающей производную
    public delegate Vector DerivativeFunction(double t, Vector x);

    // Пример: гармонические колебания (x'' = -x)
    public static Vector Oscillator(double t, Vector x)
    {
        Vector result = new Vector(2);
        result[0] = x[1];     // dx/dt = v
        result[1] = -x[0];    // dv/dt = -x
        return result;
    }

    // Метод Эйлера
    public static Matrix SolveEuler(double t0, double tEnd, Vector initialState, int steps, DerivativeFunction f)
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
            Vector derivative = f(time, state);
            state = state + derivative * h;
            time += h;

            currentRow[0] = time;
            for (int i = 0; i < dim; i++) currentRow[i + 1] = state[i];

            solution.SetColumn(step, currentRow);
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
            Vector k1 = f(time, state);
            Vector k2 = f(time + h, state + k1 * h);

            state = state + (k1 + k2) * (h / 2.0);
            time += h;

            currentRow[0] = time;
            for (int i = 0; i < dim; i++) currentRow[i + 1] = state[i];

            solution.SetColumn(step, currentRow);
        }

        return solution;
    }

    // Метод Рунге–Кутты четвёртого порядка
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
            Vector k2 = f(time + h / 2, state + k1 * h / 2);
            Vector k3 = f(time + h / 2, state + k2 * h / 2);
            Vector k4 = f(time + h, state + k3 * h);

            state = state + (k1 + 2 * k2 + 2 * k3 + k4) * (h / 6.0);
            time += h;

            currentRow[0] = time;
            for (int i = 0; i < dim; i++) currentRow[i + 1] = state[i];

            solution.SetColumn(step, currentRow);
        }

        return solution;
    }

    // Метод Адамса (четырёхшаговый)
    public static Matrix SolveAdams(double t0, double tEnd, Vector initialState, int steps, DerivativeFunction f)
    {
        int dim = initialState.Size;
        Matrix solution = new Matrix(dim + 1, steps + 1);
        double h = (tEnd - t0) / steps;

        Vector currentRow = new Vector(dim + 1);
        currentRow[0] = t0;
        for (int i = 0; i < dim; i++) currentRow[i + 1] = initialState[i];
        solution.SetColumn(0, currentRow);

        Vector[] states = new Vector[steps + 1];
        double[] times = new double[steps + 1];
        Vector[] derivatives = new Vector[steps + 1];

        times[0] = t0;
        states[0] = initialState.Copy();
        derivatives[0] = f(times[0], states[0]);

        for (int i = 1; i <= Math.Min(3, steps); i++)
        {
            times[i] = times[i - 1] + h;
            derivatives[i] = f(times[i], states[i - 1]);
            states[i] = states[i - 1] + derivatives[i] * h;

            currentRow[0] = times[i];
            for (int j = 0; j < dim; j++) currentRow[j + 1] = states[i][j];
            solution.SetColumn(i, currentRow);
        }

        for (int i = 3; i < steps; i++)
        {
            times[i + 1] = times[i] + h;
            states[i + 1] = states[i] + (55 * derivatives[i] - 59 * derivatives[i - 1] + 37 * derivatives[i - 2] - 9 * derivatives[i - 3]) * h / 24.0;
            derivatives[i + 1] = f(times[i + 1], states[i + 1]);

            currentRow[0] = times[i + 1];
            for (int j = 0; j < dim; j++) currentRow[j + 1] = states[i + 1][j];
            solution.SetColumn(i + 1, currentRow);
        }

        return solution;
    }
}
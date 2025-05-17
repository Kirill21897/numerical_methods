public class Spline
{
    private readonly double[] x;
    private readonly double[] y;

    private int intervalCount;
    private double[] coeffA, coeffB, coeffC, coeffD;

    public Spline(double[] x, double[] y)
    {
        this.x = (double[])x.Clone(); 
        this.y = (double[])y.Clone();

        if (x.Length != y.Length)
            throw new ArgumentException("Количество аргументов не совпадает с количеством значений функции");

        intervalCount = x.Length - 1;
        coeffC = new double[intervalCount];

        CreateSpline();
    }

    private void CreateSpline()
    {
        double[] intervalLengths = new double[intervalCount];
        for (int i = 0; i < intervalCount; i++)
            intervalLengths[i] = x[i + 1] - x[i];

        coeffA = new double[intervalCount];
        for (int i = 0; i < intervalCount; i++)
            coeffA[i] = y[i + 1];

        double[] diagA = new double[intervalCount - 1];
        diagA[0] = 0;
        for (int i = 1; i < intervalCount - 1; i++)
            diagA[i] = intervalLengths[i];

        double[] diagC = new double[intervalCount - 1];
        for (int i = 0; i < intervalCount - 1; i++)
            diagC[i] = 2 * (intervalLengths[i] + intervalLengths[i + 1]);

        double[] diagB = new double[intervalCount - 1];
        for (int i = 0; i < intervalCount - 2; i++)
            diagB[i] = intervalLengths[i + 1];
        diagB[intervalCount - 2] = 0;

        double[] rhsVector = new double[intervalCount - 1];
        for (int i = 0; i < intervalCount - 1; i++)
        {
            double slopeNext = (y[i + 2] - y[i + 1]) / intervalLengths[i + 1];
            double slopePrev = (y[i + 1] - y[i]) / intervalLengths[i];
            rhsVector[i] = 6 * (slopeNext - slopePrev);
        }

        double[] resultVector = TridiagonalMatrixSolver(diagA, diagC, diagB, rhsVector);

        coeffC = new double[intervalCount + 1]; 

        for (int i = 0; i < intervalCount - 1; i++)
            coeffC[i + 1] = resultVector[i];

        coeffC[0] = 0;
        coeffC[intervalCount] = 0; 

        coeffD = new double[intervalCount];
        for (int i = 0; i < intervalCount; i++)
        {
            if (i == 0)
                coeffD[i] = coeffC[i] / intervalLengths[i];
            else
                coeffD[i] = (coeffC[i] - coeffC[i - 1]) / intervalLengths[i];
        }

        coeffB = new double[intervalCount];
        for (int i = 0; i < intervalCount; i++)
        {
            coeffB[i] = (intervalLengths[i] / 2) * coeffC[i + 1]
                        - (Math.Pow(intervalLengths[i], 2) / 6) * coeffD[i]
                        + (y[i + 1] - y[i]) / intervalLengths[i];
        }
    }

    public double Interpolate(double xt)
    {
        if (xt < x[0] || xt > x[intervalCount])
            return double.NaN;

        int interval = 0;
        for (int i = 0; i < intervalCount; i++)
        {
            if (xt >= x[i] && xt <= x[i + 1])
            {
                interval = i;
                break;
            }
        }

        double dx = xt - x[interval];
        double value = y[interval] +
                       coeffB[interval] * dx +
                       (coeffC[interval] / 2) * Math.Pow(dx, 2) +
                       (coeffD[interval] / 6) * Math.Pow(dx, 3);

        return value;
    }

    private static double[] TridiagonalMatrixSolver(double[] a, double[] b, double[] c, double[] d)
    {
        int n = d.Length;
        double[] p = new double[n];
        double[] q = new double[n];
        double[] solution = new double[n];

        p[0] = -c[0] / b[0];
        q[0] = d[0] / b[0];

        for (int i = 1; i < n; i++)
        {
            double denom = b[i] + a[i - 1] * p[i - 1];
            p[i] = -c[i] / denom;
            q[i] = (d[i] - a[i - 1] * q[i - 1]) / denom;
        }

        solution[n - 1] = q[n - 1];

        for (int i = n - 2; i >= 0; i--)
        {
            solution[i] = p[i] * solution[i + 1] + q[i];
        }

        return solution;
    }
}
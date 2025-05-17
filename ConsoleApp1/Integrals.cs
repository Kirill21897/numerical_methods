public static class IntegralMethods
    {
        // Метод средних прямоугольников
        public static double MidpointIntegral(double a, double b, double eps, Func<double, double> f)
        {
            int n = 1;
            double prev = 0;

            do
            {
                n *= 2;
                double h = (b - a) / n;
                double sum = 0;

                for (int i = 0; i < n; i++)
                {
                    double mid = a + h * (i + 0.5);
                    sum += f(mid);
                }

                double current = sum * h;

                if (Math.Abs(current - prev) < eps || n > 1_000_000)
                    return current;

                prev = current;

            } while (true);
        }

        // Метод трапеций
        public static double TrapezoidIntegral(double a, double b, double eps, Func<double, double> f)
        {
            int n = 1;
            double prev = 0;

            do
            {
                n *= 2;
                double h = (b - a) / n;
                double sum = 0.5 * (f(a) + f(b));

                for (int i = 1; i < n; i++)
                {
                    sum += f(a + i * h);
                }

                double current = sum * h;

                if (Math.Abs(current - prev) < eps || n > 1_000_000)
                    return current;

                prev = current;

            } while (true);
        }

        // Метод Симпсона
        public static double SimpsonIntegral(double a, double b, double eps, Func<double, double> f)
        {
            int n = 2;
            double prev = 0;

            do
            {
                n *= 2;
                double h = (b - a) / n;
                double sumOdd = 0, sumEven = 0;

                for (int i = 1; i < n; i++)
                {
                    double x = a + i * h;
                    if (i % 2 == 1)
                        sumOdd += f(x);
                    else
                        sumEven += f(x);
                }

                double current = (h / 3) * (f(a) + f(b) + 4 * sumOdd + 2 * sumEven);

                if (Math.Abs(current - prev) < eps || n > 1_000_000)
                    return current;

                prev = current;

            } while (true);
        }

        // Двойной интеграл (метод Симпсона по обоим направлениям)
        public static double DoubleIntegral(
            double x1, double x2,
            double y1, double y2,
            double eps,
            Func<double, double, double> f)
        {
            int nx = 2, ny = 2;
            double prev = 0;

            while (true)
            {
                double hx = (x2 - x1) / nx;
                double hy = (y2 - y1) / ny;
                double total = 0;

                for (int i = 0; i <= nx; i++)
                {
                    double x = x1 + i * hx;
                    for (int j = 0; j <= ny; j++)
                    {
                        double y = y1 + j * hy;

                        int wx = (i == 0 || i == nx) ? 1 : (i % 2 == 0) ? 2 : 4;
                        int wy = (j == 0 || j == ny) ? 1 : (j % 2 == 0) ? 2 : 4;

                        total += wx * wy * f(x, y);
                    }
                }

                double current = total * hx * hy / 9;

                if (Math.Abs(current - prev) < eps)
                    return current;

                prev = current;
                nx *= 2;
                ny *= 2;
            }
        }
    }
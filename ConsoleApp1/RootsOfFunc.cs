using System;

public class EquationSolution 
{
    public static double BisectionMethod(double leftBorder, double rightBorder, double eps, Func<double, double> func)
    {
        double f_left = func(leftBorder);
        double f_right = func(rightBorder);
        if (f_left * f_right > 0)
            throw new ArgumentException("Нет корня на отрезке");

        double midle;
        do
        {
            midle = (leftBorder + rightBorder) / 2;
            double funcFromMidle = func(midle);
          //  Console.WriteLine("a={0}  b={1} c={2} fl={3} fr={4} fc={5}",leftBorder,rightBorder,midle,f_left,f_right,funcFromMidle);
            if (f_left * funcFromMidle < 0)
            {
                rightBorder = midle;
                //f_right = funcFromMidle; 
            }
            else
            {
                leftBorder = midle;
                f_left = funcFromMidle; 
            }
        } while (Math.Abs(rightBorder - leftBorder) > eps);
        return midle;
    }

    public static double NewtonMethod(double x0, double eps, Func<double, double> func)
    {
        double x = x0;
        int count_of_iter = 1;
        double delta = eps / 2;
        double dxold = 100000.0;
        while (true)
        {
            double f_x = func(x);
            //double f_x_d = derivative(x); 
            double f_x_d = (func(x + delta) - f_x) / delta;
            if (Math.Abs(f_x_d) < eps) 
            {
                throw new InvalidOperationException("Производная близка к нулю, метод может расходиться.");
            }

            double x_new = x - f_x / f_x_d;
            double dx = Math.Abs(x_new - x);

            if(Math.Abs(dx) > Math.Abs(dxold)){
                Console.WriteLine("Функция расходится");
                throw new InvalidOperationException("Функция расходится");
            }
            //Console.WriteLine($"{dxold} {dx}");
            dxold = dx;

            if (dx < eps)
            {
                Console.WriteLine($"\nКол-во итераций методом Ньютона: {count_of_iter}");
                return x_new;
            }

            x = x_new;
            count_of_iter += 1;
        }
        //return -100000.0;
    }

    public static double SimpleIterations(double x0, double eps, Func<double, double> func)
    {
        double x = x0;
        int count_of_iter = 1;
        double dxold=100000.0;
       // double delta = eps / 2;
        while (true)
        {
            double f_x = func(x);
            //double f_x_d = derivative(x); 
           // double f_x_d = (func(x + delta)-f_x) / delta;
           
            double x_new = f_x;
            double dx=x_new - x;
            //Console.WriteLine($"{dxold} {dx}");

            if(Math.Abs(dx) > Math.Abs(dxold)){
                //Console.WriteLine("Функция расходиться");
                throw new InvalidOperationException("Функция расходиться");
            }
            
            dxold = dx;

            if (Math.Abs(dx) < eps)
            {
                Console.WriteLine($"\nКол-во итераций методом SimpleIterations: {count_of_iter}");
                return x_new;
            }

            x = x_new;
            count_of_iter += 1;
        }
       // return -100000.0;
    }
    
    public static double SteffMethod(double x0, double eps, Func<double, double> func)
    {
        double x = x0;
        int count_of_iter = 1;
        double delta = eps / 2;
        double dxold = 100000.0;
        while (true)
        {
            double f_x = func(x);
            //double f_x_d = derivative(x); 
            double f_x_d = (func(x + f_x) - f_x) / f_x - 1.0; 
            if (Math.Abs(f_x_d) < eps) 
            {
                throw new InvalidOperationException("Производная близка к нулю, метод может расходиться.");
            }

            double x_new = x - f_x / f_x_d;
            double dx = Math.Abs(x_new - x);

            if(Math.Abs(dx) > Math.Abs(dxold)){
                Console.WriteLine("Функция расходится");
                throw new InvalidOperationException("Функция расходится");
            }
            //Console.WriteLine($"{dxold} {dx}");
            dxold = dx;

            if (dx < eps)
            {
                Console.WriteLine($"\nКол-во итераций методом Стеффенсена: {count_of_iter}");
                return x_new;
            }

            x = x_new;
            count_of_iter += 1;
        }
        //return -100000.0;
    }

}


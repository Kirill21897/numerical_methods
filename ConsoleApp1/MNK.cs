using ConsoleChMethod;
namespace ConsoleApp1;

class NaimQuadrat
{
    private readonly Vector argumentX;
    private readonly Vector argumentY;
    public Vector parameters;
    public delegate double PsiFunction(double value);
    public PsiFunction[] psiFunctions;
    public int argumentCount;
    public int psiFunctionCount;

    public NaimQuadrat(Vector x, Vector y, PsiFunction[] functions)
    {
        if (x.Size != y.Size)
            throw new Exception("Количество аргументов не совпадает с количеством значений функции");

        this.argumentX = x;
        this.argumentY = y;
        this.argumentCount = x.Size;
        this.psiFunctions = functions;
        this.psiFunctionCount = functions.Length;

        CalculateParameters();
    }

    private void CalculateParameters()
    {
        // Формируем матрицу H (Hankel?)
        Matrix H = new Matrix(argumentCount, psiFunctionCount);
        for (int i = 0; i < argumentCount; i++)
        {
            Vector row = GetFunctionValues(argumentX[i]);
            H.SetRow(i, row);
        }

        // H^T
        Matrix Ht = H.Transpose();

        // D = Ht * H
        Matrix D = Ht * H;

        // D_inverse = D^(-1)
        Matrix D_Inverse = Matrix.obrSolveGaussian(D); // или ObrSolveUpperTriangular(), если треугольная

        // Q = D_inverse * Ht
        Matrix Q = D_Inverse * Ht;

        // parameters = Q * argumentY
        parameters = Q * argumentY;
    }

    public Vector GetFunctionValues(double x)
    {
        Vector result = new Vector(psiFunctionCount);
        for (int i = 0; i < psiFunctionCount; i++)
        {
            result[i] = psiFunctions[i](x);
        }
        return result;
    }

    public double GetCriterion()
    {
        Vector residuals = new Vector(argumentCount);
        for (int i = 0; i < argumentCount; i++)
        {
            Vector phi = GetFunctionValues(argumentX[i]);
            residuals[i] = argumentY[i] - (parameters * phi); // скалярное произведение
        }
        return residuals.Norma1(); // L2 норма невязки
    }
}
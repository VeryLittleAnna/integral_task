#include <stdio.h>
#include <stdlib.h>
#include <math.h>

long double f1(long double x)
{
    return x;
}

long double f1_derivative(long double x)
{
    return 1;
}

long double f2(long double x)
{
    return 3 - 0.5 * x;
}

long double f2_derivative(long double x)
{
    return -0.5;
}

long double f3(long double x)
{
    return x * x / 16;
}

long double f3_derivative(long double x)
{
    return x / 8;
}

long double f4 (long double x)
{
    return x*x;
}

long double f4_derivative(long double x)
{
    return 2 * x;
}
long double f5(long double x)
{
    return 1 / x;
}

long double f5_derivative(long double x)
{
    return -1 / (x * x);
}



long double root(long double (* f) (long double), long double (*g) (long double), long double (* f_derivative) (long double), long double (* g_derivative) (long double), long double a, long double b, long double eps1, int * cnt)
{ ///function F(x) = f(x) - g(x)
    *cnt = 0;
    long double c = a, d = b, func_value_c = f(a) - g(a), func_value_d = f(b) - g(b); /// [c, d] - current section, func_value_c = F(c), func_value_d = F(d)
    int type = (f_derivative(a) - g_derivative(a)) * (f_derivative(b) - g_derivative(b) - f_derivative(a) + g_derivative(a)) >= 0; /// F'(a) * (F'(b) - F'(a)) > 0 <=> F' * F" > 0f_derivative(a) - g_derivative(a)
    while (d - c > eps1)
    {
        if (*cnt % 2 == 0) ///iteration of secant method
        {
            if (type) /// x[n] < x0 , x0 - exact root, x[n] - x after nth iteration
            {
                c -= (d - c) / (func_value_d - func_value_c) * func_value_c;
                func_value_c = f(c) - g(c);
            }
            else /// x[n] > x0
            {
                d -= (d - c) / (func_value_d - func_value_c) * func_value_d;
                func_value_d = f(d) - g(d);
            }
        }
        else ///iteration of Newton's method
        {
            if (type)/// x[n] > x0
            {
                d -= func_value_d / (f_derivative(d) - g_derivative(d));
                func_value_d = f(d) - g(d);
            }
            else /// x[n] < x0
            {
                c -= func_value_c / (f_derivative(c) - g_derivative(c));
                func_value_c = f(c) - g(c);
            }
        }
        (*cnt)++;
    }
    return (c + d) / 2;
}



long double integral( long double (*f) (long double), long double a, long double b, long double eps2)
{
    int n = 20; /// n0 - number of sections at the beginning
    long double cur_delta = 2 * eps2, I = 0, next_I = 0, tmp; /// I - current integral
    long double h = (b - a) / n, sum_of_cur_F = 0, sum_of_next_F = 0, f_a = f(a), f_b = f(b); /// h - length of sections, sum_of_cur_F = F(2) + F(4) + ... F(2 * n - 2), sum_of_next_F = F(1] + F(3) + ... + F(2 * N - 1)
    for (int i = 1; i < n; ++i)
    {
        tmp = f(a + i * h);
        sum_of_cur_F += tmp;
        I += tmp * (i % 2 ? 2 : 4);
    }
    I = (I + f_a + f_b) * h / 3;
    cur_delta = I;
    while (15 * cur_delta > eps2) ///Runge's rule
    {
        next_I = 0;
        h /= 2;
        n *= 2;
        sum_of_next_F = 0;
        for (int i = 1; i < n; i += 2)
        {
            sum_of_next_F += f(a + i * h);
        }
        next_I = (f_a + f_b + sum_of_cur_F * 2 + sum_of_next_F * 4) * h / 3;
        sum_of_cur_F += sum_of_next_F;
        cur_delta = fabs(next_I - I);
        I = next_I;
    }
    return fabs(I);
}

#define A -1
#define B 6
#define eps 0.000001

int main(int argc, char * argv[])
{
    long double ans = 0;
    int cnt = 0;
    ans = root(f1, f2, f1_derivative, f2_derivative, A, B, eps, &cnt); ///2
    printf("%Lf %d\n", ans, cnt);
    ans = root(f2, f3, f2_derivative, f3_derivative, A, B, eps, &cnt); ///4
    printf("%Lf %d\n", ans, cnt);
    ans = root(f1, f3, f1_derivative, f3_derivative, 10, 19, eps, &cnt); ///16
    printf("%Lf %d\n", ans, cnt);
    ans = root(f4, f1, f4_derivative, f1_derivative, 0.5, B, eps, &cnt); ///1
    printf("%Lf %d\n", ans, cnt);
    ans = root(f4, f1, f4_derivative, f1_derivative,-4, 0.6, eps, &cnt); /// 0
    printf("%Lf %d\n", ans, cnt);
    ans = root(f5, f4, f5_derivative, f4_derivative, 0.2, 2, eps, &cnt); ///1
    printf("%Lf %d\n", ans, cnt);

    printf("\n\n");
    printf("%Lf\n", integral(f1, A, B, eps)); ///17.5
    printf("%Lf\n", integral(f2, A, B, eps)); ///12.25
    printf("%Lf\n", integral(f3, -4, 3, eps)); ///1.8958
    printf("%Lf\n", integral(f4, 3, 5, eps)); ///32.667
    printf("%Lf\n", integral(f5, 1, 2.71828, eps)); ///1
    return 0;
}

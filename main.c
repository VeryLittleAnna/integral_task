#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef test
long double f1(long double x)
{
    return 0.6 * x + 3;
}

long double f1_derivative(long double x)
{
    return 0.6;
}

long double f2(long double x)
{
    return (x - 2) * (x - 2) * (x - 2) - 1;
}

long double f2_derivative(long double x)
{
    return 3 * (x - 2) * (x - 2);
}

long double f3(long double x)
{
    return 3.0 / x;
}

long double f3_derivative(long double x)
{
    return - 3.0 / x / x;
}
#else
extern long double f1(long double x);
extern long double f1_derivative(long double x);
extern long double f2(long double x);
extern long double f2_derivative(long double x);
extern long double f3(long double x);
extern long double f3_derivative(long double x);
#endif



long double root(long double (* f) (long double), long double (*g) (long double), long double (* f_derivative) (long double), long double (* g_derivative) (long double), long double a, long double b, long double eps1, int * cnt)
{ ///function F(x) = f(x) - g(x)
    *cnt = 0;
    long double c = a, d = b, func_value_c = f(a) - g(a), func_value_d = f(b) - g(b); /// [c, d] - current section, func_value_c = F(c), func_value_d = F(d)
    int type = (f_derivative(a) - g_derivative(a)) * (f_derivative(b) - g_derivative(b) - f_derivative(a) + g_derivative(a)) >= 0; /// F'(a) * (F'(b) - F'(a)) > 0 <=> F' * F" > 0f_derivative(a) - g_derivative(a)
    while (fabs(d - c) > eps1)
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

#define eps 0.000001

int main(int argc, char * argv[])
{
	//printf("%d\n", argc);
    long double ans = 0;
    int cnt = 0;

    printf("%Lf %Lf %Lf %Lf %Lf\n", f1(1), f1(2), f1(0), f1(4), f1_derivative(1)); ///3.6  4.2  3  5.4  0.6
    printf("%Lf %Lf %Lf %Lf %Lf %Lf %Lf\n", f2(1), f2(2), f2(0), f2(4), f2_derivative(0), f2_derivative(1), f2_derivative(-1)); ///
    ans = root(f1, f2, f1_derivative, f2_derivative, 2.5, 4, eps, &cnt);
    printf("%Lf %d\n", ans, cnt);
    return 0;
}

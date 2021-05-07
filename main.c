#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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

#ifdef COMBINE

long double root(long double (* f) (long double), long double (*g) (long double), long double (* f_derivative) (long double), long double (* g_derivative) (long double), long double a, long double b, long double eps1, int * cnt)
{ ///function F(x) = f(x) - g(x)
    *cnt = 0;
    long double c = a, d = b, func_value_c = f(a) - g(a), func_value_d = f(b) - g(b); /// [c, d] - current section, func_value_c = F(c), func_value_d = F(d)
    int type = (f_derivative(a) - g_derivative(a)) * (f_derivative(b) - g_derivative(b) - f_derivative(a) + g_derivative(a)) >= 0; /// F'(a) * (F'(b) - F'(a)) > 0 <=> F' * F" > 0f_derivative(a) - g_derivative(a)
    while (fabsl(d - c) > eps1)
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
 #else

long double root(long double (* f) (long double), long double (*g) (long double), long double (* f_derivative) (long double), long double (* g_derivative) (long double), long double a, long double b, long double eps1, int * cnt)
{ ///function F(x) = f(x) - g(x)
    long double c = a, d = b, func_c = f(a) - g(a), func_d = f(b) - g(b), x, func_x; ///current section [c, d], F(c) * F(d) < 0
    *cnt = 0;
    while (d - c > eps1)
    {
        x = (c + d) / 2;
        func_x = f(x) - g(x);
        if (func_x * func_c >= 0) c = x;
        else d = x;
        (*cnt)++;
    }
    return (c + d) / 2;
}
#endif


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
    cur_delta = fabsl(I);
    int cnt = 0;
    while (cur_delta > eps2 / 15) ///Runge's rule
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
        cur_delta = fabsl(next_I - I);
        I = next_I;
        ++cnt;
    }
    return I;
}

#define eps_2 0.0003
#define eps_1 0.00001
#define KEYS_NUM 5

char * keys[KEYS_NUM] = {"-help", "-test_i", "-test_r", "-cross", "-steps"};
char * meanings[KEYS_NUM] = {"\n    print all keys", "INDEX_OF_FUNCTION L R\n    find integral under function with entered index (1/2/3) on [L, R]", \
"INDEX_OF_FUNCTION INDEX_OF_FUNCTION L R\n    find intersections of two functions with entered indices (1/2/3) on [L, R]", \
"\n    print points of function intersection", "\n    print number of iterations during root calculations"};

void print_keys(void)
{
    for (int i = 0; i < KEYS_NUM; ++i)
    {
        printf("%s %s\n", keys[i], meanings[i]);
    }
}

void test_integral(int index, long double a, long double b)
{
    long double result = 0;
    switch (index)
    {
        case 1: result = integral(f1, a, b, eps_2);
            break;
        case 2: result = integral(f2, a, b, eps_2);
            break;
        case 3: result = integral(f3, a, b, eps_2);
            break;
    }
    printf("integral of function %d on [%Lf, %Lf] = %Lf\n", index, a, b, result);
}

void test_root(int index1, int index2, long double a, long double b)
{
    int indices = 10 * (index1 < index2 ? index1 : index2) + (index1 < index2 ? index2 : index1); ///min(index1, index2) * 10 + max(index1, index2)
    long double result = 0;
    int cnt = 0;
    switch (indices)
    {
        case 12: result = root(f1, f2, f1_derivative, f2_derivative, a, b, eps_1, &cnt); /// 1 and 2
            break;
        case 13: result = root(f1, f3, f1_derivative, f3_derivative, a, b, eps_1, &cnt); /// 1 and 3
            break;
        case 23: result = root(f2, f3, f2_derivative, f3_derivative, a, b, eps_1, &cnt); /// 2 and 3
            break;
    }
    printf("root of F = f%d-f%d on [%Lf, %Lf] = %Lf\n with %d iterations\n", index1, index2, a, b, result, cnt);
}

int main(int argc, char * argv[])
{
    if (argc > 1 && strcmp(argv[1], "-help") == 0)
    {
        print_keys();
        return 0;
    }
    else if (argc > 1 && strcmp(argv[1], "-test_i") == 0)
    {
        int index = argv[2][0] - '0';
        long double a = strtold(argv[3], NULL), b = strtold(argv[4], NULL);
        test_integral(index, a, b);
        return 0;
    }
    else if (argc > 1 && strcmp(argv[1], "-test_r") == 0)
    {
        int index1 = argv[2][0] - '0', index2 = argv[3][0] - '0';
        long double a = strtold(argv[4], NULL), b = strtold(argv[5], NULL);
        test_root(index1, index2, a, b);
        return 0;
    }
    long double ans = 0, r1, r2, r3;
    int cnt1, cnt2, cnt3;
    r1 = root(f1, f2, f1_derivative, f2_derivative, 2.5, 4, eps_2, &cnt1); ///3.848
    r2 = root(f3, f2, f3_derivative, f2_derivative, 2.5, 4, eps_2, &cnt2); ///3.244
    r3 = root(f1, f3, f1_derivative, f3_derivative, 0.01, 4, eps_2, &cnt3); ///0.854
    if (argc > 1 && strcmp(argv[1], "-cross") == 0 || argc > 2 && strcmp(argv[2], "-cross") == 0)
    {
        printf("f1 crosses f2 in %Lf\n", r1);
        printf("f2 crosses f3 in %Lf\n", r2);
        printf("f1 crosses f3 in %Lf\n", r3);
    }
    if (argc > 1 && strcmp(argv[1], "-steps") == 0 || argc > 2 && strcmp(argv[2], "-steps") == 0)
    {
        printf("for root1 [f1 cross f2] - %d\n", cnt1);
        printf("for root2 [f2 cross f3] - %d\n", cnt2);
        printf("for root3 [f1 cross f2] - %d\n", cnt3);
    }

    long double result = integral(f1, r3, r1, eps_1);
    result += integral(f2, r1, r2, eps_1);
    result += integral(f3, r2, r3, eps_1);
    printf("ANSWER: %Lf\n", fabsl(result));
    return 0;
}

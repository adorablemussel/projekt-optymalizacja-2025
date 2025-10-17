#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejœciowe:
	// ff - wskaŸnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zak³¹dana dok³adnoœæ rozwi¹zania
	// Nmax - maksymalna liczba wywo³añ funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartoœæ funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi¹zanie z zadan¹ dok³adnoœci¹
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo³añ funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double *p = new double[2] { 0, 0 };
		double x1 = x0 + d;
		int i = 0;
		int fcalls = 0;
		//cout << ff(5.0, ud1, ud2) << endl;
		//cout << ff(4.0, ud1, ud2) << endl;
		//cout << ff(-0.00327869, ud1, ud2) << endl;
		if (ff(x0, ud1, ud2) == ff(x1, ud1, ud2))
		{
			p[0] = x0;
			p[1] = x1;
			return p;
		}
		if (ff(x1, ud1, ud2) > ff(x0, ud1, ud2))
		{
			d = -d;
			x1 = x0 + d;
			if (ff(x1, ud1, ud2) >= ff(x0, ud1, ud2))
			{
				p[0] = x1;
				p[1] = x0 - d;
				return p;
			}
		}
		double xi = x1;
		double xi_next = 0.0;
		while (fcalls <= Nmax)
		{
			i++;
			xi_next = x0 + pow(alpha, i) * d;
			fcalls++;

			if (ff(xi, ud1, ud2) <= ff(xi_next, ud1, ud2))
			{
				if (d > 0)
				{
					p[0] = x0 + pow(alpha, i - 1) * d;
					p[1] = xi_next;
				}
				else
				{
					p[0] = xi_next;
					p[1] = x0 + pow(alpha, i - 1) * d;
				}
				return p;
			}
			xi = xi_next;

		}


		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
       int k = 0;
		double length_ratio = (b - a) / epsilon;
       while (fib_num(k) <= length_ratio)
       {
           k++;
       }
       double Fk1 = (double)fib_num(k - 1);
       double Fk = (double)fib_num(k);

       double a_i = a;
       double b_i = b;
	   double c_i = b_i - (Fk1 / Fk) * (b_i - a_i);
	   double d_i = a_i + b_i - c_i;

       solution Xc(c_i),Xd(d_i);
       Xc.fit_fun(ff, ud1, ud2);
       Xd.fit_fun(ff, ud1, ud2);
       for (int i = 0; i <= k - 3; ++i)
       {
           if (m2d(Xc.y) < m2d(Xd.y))
           {
               b_i = m2d(Xd.x);
               Xd = Xc;

               Fk1 = (double)fib_num(k - i - 2);
               Fk = (double)fib_num(k - i - 1);
               c_i = b_i - (Fk1 / Fk) * (b_i - a_i);
               Xc = solution(c_i);
               Xc.fit_fun(ff, ud1, ud2);
           }
           else
           {
               a_i = m2d(Xc.x);
			   Xc = Xd;
               Fk1= (double)fib_num(k - i - 2);
        		Fk = (double)fib_num(k - i - 1);
               d_i = a_i + b_i - c_i;
               Xd = solution(d_i);
               Xd.fit_fun(ff, ud1, ud2);
           }
       }
       Xopt = Xc;
       Xopt.fit_fun(ff, ud1, ud2);
       Xopt.flag = 1;
       return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;
        auto f = [&](double x)->double {
            matrix X(1, 1, x);
            return m2d(ff(X, ud1, ud2));
        };

        double aa = a, bb = b;
        double c = 0.5 * (aa + bb);

        double fa = f(aa), fb = f(bb), fc = f(c);
        int fcalls = 3;
        double prev_d = NAN;
        double d = NAN;

        while (true)
        {
            if (fcalls > Nmax)
                throw string("lag: przekroczono Nmax");

            double l = fa * (bb*bb - c*c) + fb * (c*c - aa*aa) + fc * (aa*aa - bb*bb);
            double m = fa * (bb - c) + fb * (c - aa) + fc * (aa - bb);

            if (fabs(m) < 1e-12) // zabezpieczenie przed dzieleniem przez zero
            {
                d = c; // ustaw na środek
            }
            else
            {
                d = 0.5 * l / m;
                // ograniczenie d do przedziału [aa, bb]
                if (d <= aa) d = aa + 1e-12;
                if (d >= bb) d = bb - 1e-12;
            }

            double fd = f(d); ++fcalls;

            // aktualizacja przedziału
            if (aa < d && d < c)
            {
                if (fd < fc)
                {
                    bb = c; fb = fc;
                    c = d; fc = fd;
                }
                else
                {
                    aa = d; fa = fd;
                }
            }
            else if (c < d && d < bb)
            {
                if (fd < fc)
                {
                    aa = c; fa = fc;
                    c = d; fc = fd;
                }
                else
                {
                    bb = d; fb = fd;
                }
            }
            else
            {
                // jeśli d wciąż wypadł poza przedział, przypisz do najbliższego krańca
                if (d <= aa) d = aa;
                else if (d >= bb) d = bb;
            }

            // warunki stopu
            if ((bb - aa) < epsilon)
            {
                Xopt.x = matrix(1, 1, 0.5 * (aa + bb));
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                return Xopt;
            }

            if (!isnan(prev_d) && fabs(d - prev_d) < gamma)
            {
                Xopt.x = matrix(1, 1, d);
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                return Xopt;
            }

            prev_d = d;
        }
    }
    catch (const string& ex_info)
    {
        throw string("solution lag(...):\n" + ex_info);
    }
}


solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}

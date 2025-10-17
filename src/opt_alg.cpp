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

       // Krok 2: a(0) = a, b(0) = b.
       // Krok 1: Deklaracja kluczowych zmiennych PRZED pętlą, aby były dostępne w całym zakresie.
       double a_i = a;
       double b_i = b;
       // Założenie dla punktu wewnętrznego c(0) (wg pseudokodu musi być podany)
       // Jeśli nie jest podany, używamy środka:
       double c_i = (a + b) / 2.0;

       // Deklaracja zmiennych d_i i d_i_prev przed pętlą
       double d_i;
       double d_i_prev = c_i; // Wymagane do kryterium |d(i) - d(i-1)| < gamma

       // Deklaracja obiektów solution przed pętlą
       solution Xa(a_i), Xb(b_i), Xc(c_i), Xd;

       // Wymagane obliczenia początkowe dla a, b, c
       Xa.fit_fun(ff, ud1, ud2);
       Xb.fit_fun(ff, ud1, ud2);
       Xc.fit_fun(ff, ud1, ud2);


       // Krok 3: repeat (i = 0 do nieskonczoności, z kryteriami stopu)
       for (int i = 0; solution::f_calls < Nmax; ++i)
       {
           // Sprawdzenie kryterium stopu na początku pętli
           // Krok 39: until b(i) – a(i) < ? or |d(i) – d(i-1)| < ?
           if (i > 0 && ((b_i - a_i < epsilon) || (fabs(d_i - d_i_prev) < gamma)))
           {
               Xopt.flag = 1; // Optymalizacja zakończona sukcesem
               break;
           }

           // Wymagane są aktualne wartości funkcji w punktach a, b, c
           double fa = m2d(Xa.y), fb = m2d(Xb.y), fc = m2d(Xc.y);

           // Krok 4 & 5: Obliczenie l i m
           double l = fa * (pow(b_i, 2) - pow(c_i, 2)) + fb * (pow(c_i, 2) - pow(a_i, 2)) + fc * (pow(a_i, 2) - pow(b_i, 2));
           double m = fa * (b_i - c_i) + fb * (c_i - a_i) + fc * (a_i - b_i);

           // Krok 6 & 7: if m <= 0 then return error
           if (m <= 0.0)
           {
               Xopt.flag = -1;
               Xopt.x = c_i;
               Xopt.fit_fun(ff, ud1, ud2);
               return Xopt;
           }

           // Krok 9: d(i) = 0,5 * l / m
           d_i_prev = d_i; // Zapis d(i-1) przed obliczeniem d(i)
           d_i = 0.5 * l / m;
           Xd.x = d_i;

           // Sprawdzenie, czy d(i) jest w przedziale [a,b]
           if (a_i < d_i && d_i < b_i)
           {
               Xd.fit_fun(ff, ud1, ud2);

               // Krok 10: if a(i) < d(i) < c(i) then
               if (d_i < c_i)
               {
                   // Krok 11-14: if f(d(i)) < f(c(i))
                   if (m2d(Xd.y) < m2d(Xc.y))
                   {
                       // a(i+1)=a(i), c(i+1)=d(i), b(i+1)=c(i)
                       b_i = c_i;
                       c_i = d_i;
                       Xb = Xc;
                       Xc = Xd;
                   }
                   else // Krok 15-18
                   {
                       // a(i+1)=d(i), c(i+1)=c(i), b(i+1)=b(i)
                       a_i = d_i;
                       Xa = Xd;
                   }
               }
               // Krok 21: else if c(i) < d(i) < b(i) then
               else // (d_i > c_i)
               {
                   // Krok 22-25: if f(d(i)) < f(c(i))
                   if (m2d(Xd.y) < m2d(Xc.y))
                   {
                       // a(i+1)=c(i), c(i+1)=d(i), b(i+1)=b(i)
                       a_i = c_i;
                       c_i = d_i;
                       Xa = Xc;
                       Xc = Xd;
                   }
                   else // Krok 26-29
                   {
                       // a(i+1)=a(i), c(i+1)=c(i), b(i+1)=d(i)
                       b_i = d_i;
                       Xb = Xd;
                   }
               }
           }
           // Krok 31-33: else return error (d(i) jest poza [a, b])
           else
           {
               Xopt.flag = -1;
               Xopt.x = c_i;
               Xopt.fit_fun(ff, ud1, ud2);
               return Xopt;
           }

           // Krok 36: Sprawdzenie Nmax
           if (solution::f_calls >= Nmax)
           {
               Xopt.flag = 0; // Oznaczenie osiągnięcia Nmax
               break;
           }
       		//std::cout << i + 1 << " " << (b_i - a_i) << std::endl;
       } // end repeat

       // Krok 40: return x* = d(i)
       Xopt.x = d_i;
       Xopt.fit_fun(ff, ud1, ud2);
       return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
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

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
		//Tu wpisz kod funkcji
		//std::cout << std::setprecision(8);
		std::vector<int> phi = { 1,1 };
		int k = 2;
		while(true) {
			int j = phi[k - 1] + phi[k - 2];
			phi.push_back(j);
			if (j > (b - a) / epsilon) break;
			k++;
		}
		std::vector<double> ai = { a };
		std::vector<double> bi = { b };
		std::vector<double> ci;
		std::vector<double> di;

		ci.push_back(b - double(phi[k - 1]) / phi[k] * (b - a));
		di.push_back(a + b - ci[0]);
		std::cout << ai[0] << " " << bi[0] << " " << ci[0] << " " << di[0] << "\n";
		for (int i = 0; i<=k - 3; i++) {
			if (ff(ci[i],ud1, ud2) < ff(di[i], ud1, ud2)) {
				ai.push_back(ai[i]);
				bi.push_back(di[i]);
			}
			else {
				bi.push_back(bi[i]);
				ai.push_back(ci[i]);
			}
			ci.push_back(bi[i+1] - double(phi[k - i - 2]) / phi[k - i - 1] * (bi[i+1] - ai[i+1]));
			di.push_back(ai[i+1] + bi[i+1] - ci[i+1]);
		};
		Xopt = ci[k - 3 + 1];
		Xopt.fit_fun(ff,ud1,ud2);
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

		// brak argumentu c w sygnaturze - użyjemy środka jako c(0)
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

			// oblicz l i m zgodnie z pseudokodem
			double l = fa * (bb * bb - c * c) + fb * (c * c - aa * aa) + fc * (aa * aa - bb * bb);
			double m = fa * (bb - c) + fb * (c - aa) + fc * (aa - bb);

			if (m <= 0)
				throw string("lag: m <= 0 (blad interpolacji)");

			d = 0.5 * l / m;

			if (aa < d && d < c)
			{
				double fd = f(d); ++fcalls;
				if (fd < fc)
				{
					// a(i+1) = a(i); c(i+1) = d(i); b(i+1) = c(i)
					bb = c; fb = fc;
					c = d; fc = fd;
				}
				else
				{
					// a(i+1) = d(i); c(i+1) = c(i); b(i+1) = b(i)
					aa = d; fa = fd;
				}
			}
			else if (c < d && d < bb)
			{
				double fd = f(d); ++fcalls;
				if (fd < fc)
				{
					// a(i+1) = c(i); c(i+1) = d(i); b(i+1) = b(i)
					aa = c; fa = fc;
					c = d; fc = fd;
				}
				else
				{
					// a(i+1) = a(i); c(i+1) = c(i); b(i+1) = d(i)
					bb = d; fb = fd;
				}
			}
			else
			{
				throw string("lag: d poza przedzialem");
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

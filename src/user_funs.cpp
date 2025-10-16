#include"user_funs.h"
matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera wartoœæ funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wspó³rzêdne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera wartoœæ funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	int n = get_len(Y[0]);									// d³ugoœæ rozwi¹zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad³a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// wartoœæ funkcji celu (ud1 to za³o¿one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
	return y;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * x(0)) * exp(-(0.1 * x(0) - 2 * 3.14) * (0.1 * x(0) - 2 * 3.14)) + 0.002 * (0.1 * x(0) * (0.1 * x(0)));

	return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0(3,1);
	Y0(0)=5.0;	//objętość zbiornika a
	Y0(1)=1.0;	//objętość zbiornika b
	Y0(2)=20.0;	//temperatura zbiornika b
	ud1(0) = m2d(x);
	matrix* Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, ud1, ud2);
	int n = get_len(Y[0]);
	double maxT_B = Y[1](0,2);
	for(int i=1;i<n;i++)
	{
		if(Y[1](i,2)>maxT_B)
		{
			maxT_B=Y[1](i,2);			
		}
	}
	y=abs(maxT_B - m2d(ud1));
	Y[0].~matrix();											
	Y[1].~matrix();
	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(3,1);
	double a = 0.98;
	double b = 0.63; 
	double g = 9.81; // zmiana objętości wody w zbiorniku
	double PA = 0.5;
	double TA = 90.0; // zbiornik a
	double PB = 1.0; // zbiornik b
	double Tin = 20.0;
	double Fin = 0.01; // wlewanie do zbiornika b F in = 10 litrów/s
	//double db = 36.5665; // przekrój otwóru wylewającej się wody
	//double da = 50.0;
	//db=db*1e-4;
	//double da=ud2(0) * 1e-4;
	double DB = 0.00365665;
	double DA = ud1(0);
	double VA=Y(0);
	double VB=Y(1);
	double TB=Y(2);

	double FAout = VA > 0 ? a * b * DA * sqrt(2 * g * VA / PA) : 0;
    double FBout = VB > 0 ? a * b * DB * sqrt(2 * g * VB / PB) : 0;
	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout; 
	dY(2) = (FAout / VB) * (TA - TB) + (Fin / VB) * (Tin - TB);

	return dY;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po³o¿enia to prêdkoœæ
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z prêdkoœci to przyspieszenie
	return dY;
}
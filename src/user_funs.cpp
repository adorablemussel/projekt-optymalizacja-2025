#include"../include/user_funs.h"

#define PI 3.141592653589793
#define E 2.718281828459045


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

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po³o¿enia to prêdkoœæ
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z prêdkoœci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
		matrix y;
		y = -cos(0.1 * x(0)) * pow(E,-(pow((0.1 * x(0) - 2 * PI), 2))) + 0.002 * pow(0.1 * x(0), 2);
		return y;

}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0(3,1);
	Y0(0)=5.0;				//objętość zbiornika a
	Y0(1)=1.0;				//objętość zbiornika b
	Y0(2)=20.0;				//temperatura zbiornika b
	matrix data(1,1);
	data(0)=m2d(x);
	matrix* Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, data, ud2);
	int n = get_len(Y[0]);
	double maxT_B = Y[1](0,2);
	for(int i=1;i<n;i++)
	{
		if(Y[1](i,2)>maxT_B)
		{
			maxT_B=Y[1](i,2);			
		}
	}
	y=fabs(maxT_B - 50);
	Y[0].~matrix();											
	Y[1].~matrix();
	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(3,1);
	double a = 0.98; 		// współczynnik wypływu
	double b = 0.63; 		// współczynnik wypływu
	double g = 9.81; 		// przyspieszenie ziemskie
	double PA = 2.0; 		// ciśnienie w zbiorniku a
	double TA = 95.0; 		// temperatura w zbiorniku a
	double PB = 1.0; 		// ciśnienie w zbiorniku b
	double Tin = 20.0;		// temperatura dopływającej cieczy
	double Fin = 0.01;		// przepływ dopływającej cieczy
	double DB = 0.00365665;	// gęstość cieczy w zbiorniku b
	double DA = m2d(ud1);	// gęstość cieczy w zbiorniku a
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


matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	double x1 = x(0);
	double x2 = x(1);

	double result = pow(x1, 2) + pow(x2, 2) - cos(2.5 * 3.14 * x1) - cos(2.5 * 3.14 * x2) + 2;
	//double result = pow(x1 + 3, 2) + pow(x2 - 2, 2);

	return matrix(1, 1, result);
}

matrix ff2R(matrix x, matrix ud1, matrix ud2) {
	matrix y = 0;
	matrix Y0(2, 1), Yref(2, new double[2] {3.14, 0});
	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Yref, x);
	int n = get_len(Y[0]);
	for (int i = 0; i < n; i++) {
		y = y + 10 * pow(Yref(0) - Y[1](i, 0), 2) + pow(Yref(1) - Y[1](i, 1), 2) + pow(x(0) * (Yref(0) - Y[1](i, 0)) + x(1) * (Yref(1) - Y[1](i, 1)), 2);
	}
	y = y * 0.1;
	return y;
}


matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {
	double mr = 1.0;				//masa ramienia
	double mc = 5.0;				//masa ciezarka
	double l = 1;					//dl. ramienia
	double alfa_ref = ud1(0);		//pi rad
	double omega_ref = ud1(1);		//0 rad/s
	double b = 0.5;					//wsp. tarcia

	double I = (mr * l * l) / 3 + mc * l * l;

	double k1 = ud2(0);			//wsp. wzmocnienia
	double k2 = ud2(1);			//przesylane w ud

	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (k1 * (alfa_ref - Y(0)) + k2 * (omega_ref - Y(1)) - b * Y(1)) / I;
	return dY;
}

matrix ff3T_out(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = (sin(3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)))) / (3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)));

	double a = m2d(ud1);
	double c = m2d(ud2);

	//g1
	if (-x(0) + 1 > 0)
		y = y + c * pow(-x(0) + 1, 2);
	//g2
	if (-x(1) + 1 > 0)
		y = y + c * pow(-x(1) + 1, 2);
	//g3
	if (norm(x) - a > 0)
		y = y + c * pow(norm(x) - a, 2);

	return y;
}

matrix ff3T_in(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = (sin(3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)))) / (3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)));

	double a = m2d(ud1);
	double c = m2d(ud2);

	//g1
	if (-x(0) + 1 > 0)
		y = 1E10;
	else
		y = y - c / (-x(0) + 1);
	//g2
	if (-x(1) + 1 > 0)
		y = 1E10;
	else
		y = y - c / (-x(1) + 1);
	//g3
	if (norm(x) - a > 0)
		y = 1E10;
	else
		y = y - c / (norm(x) - a);

	return y;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    matrix Y0(4, new double[4] {0.0, x(0), 100.0, 0.0}); 
    matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, x(1));
    
    int n = get_len(Y[0]);
    int i50 = 0;
    int i_0 = 0;

    for (int i = 0; i < n; ++i)
    {
        if (abs(Y[1](i, 2) - 50.0) < abs(Y[1](i50, 2) - 50.0))
            i50 = i;
        if (abs(Y[1](i, 2)) < abs(Y[1](i_0, 2)))
            i_0 = i;
    }
    y = -Y[1](i_0, 0);

    // --- FUNKCJE KARY ---
    if (abs(x(0)) - 10 > 0)
        y = y + ud2 * pow(abs(x(0)) - 10, 2);
    if (abs(x(1)) - 10 > 0) 
        y = y + ud2 * pow(abs(x(1)) - 10, 2);
    if (abs(Y[1](i50, 0) - 5.0) - 2.0 > 0)
        y = y + ud2 * pow(abs(Y[1](i50, 0) - 5.0) - 2.0, 2);

    return y;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
	//Wektor zmian po czasie
	matrix dY(4, 1);

	//Zmienne dane
	double x = Y(0);
	double v_x = Y(1);
	double y = Y(2);
	double v_y = Y(3);

	//Dane zadania
	double C = ud1(0);
	double rho = ud1(1);
	double r = ud1(2);
	double m = ud1(3);
	double g = ud1(4);

	double s = 3.14 * pow(r, 2);

	double D_x = (1.0 / 2.0) * C * rho * s * v_x * abs(v_x);
	double D_y = (1.0 / 2.0) * C * rho * s * v_y * abs(v_y);
	double FM_x = rho * v_y * m2d(ud2) * 3.14 * pow(r, 3);
	double FM_y = rho * v_x * m2d(ud2) * 3.14 * pow(r, 3);

	dY(0) = v_x;
	dY(1) = (-D_x - FM_x) / m;
	dY(2) = v_y;
	dY(3) = ((-m * g) - D_y - FM_y) / m;

	return dY;
}

matrix ff4T(matrix x, matrix ud1, matrix ud2)
{
    // zakładamy, że x jest wektorem 2x1: [x1; x2]
    double x1 = x(0);
    double x2 = x(1);

    return (1.0 / 6.0) * pow(x1, 6) - 1.05 * pow(x1, 4) + 2.0 * pow(x1, 2) + pow(x2, 2) + x1 * x2;
}

matrix gf4T(matrix x, matrix ud1, matrix ud2)
{
    matrix g(2, 1);
    double x1 = x(0);
    double x2 = x(1);

    g(0) = pow(x1, 5) - 4.2 * pow(x1, 3) + 4.0 * x1 + x2;
    g(1) = 2.0 * x2 + x1;

    return g;
}

matrix Hf4T(matrix x, matrix ud1, matrix ud2)
{
    matrix H(2, 2);
    double x1 = x(0);

    H(0, 0) = 5.0 * pow(x1, 4) - 12.6 * pow(x1, 2) + 4.0;
    H(0, 1) = 1.0;
    H(1, 0) = 1.0;
    H(1, 1) = 2.0;

    return H;
}


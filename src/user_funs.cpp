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
	matrix ud2_local(1,1);
	ud2_local(0)=x(0);
	matrix* Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, ud1, ud2_local);
	int n = get_len(Y[0]);
	//cout<<*get_size(*Y)<<endl;
	double maxT_B = Y[1](0,0);
	for(int i=1;i<n;i++)
	{
		if(m2d(Y[1](i,0))>maxT_B)
		{
			maxT_B=m2d(Y[1](i,0));
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
	double a = 0.98, b = 0.63, g = 9.81; // zmiana objętości wody w zbiorniku
	double pa = 2, ta = 95; // zbiornik a
	double pb = 1; // zbiornik b
	double tin = 20, vin = 0.01; // wlewanie do zbiornika b F in = 10 litrów/s
	double db = 36.5665; // przekrój otwóru wylewającej się wody
	db=db*1e-4;
	double da=ud2(0) * 1e-4;
	double V_A=Y(0);
	double V_B=Y(1);
	double T_B=Y(2);
	dY(0) = -a * b * da * sqrt(2*g*(V_A/pa));
	dY(1) = -a * b * db * sqrt(2*g*(V_B/pb)) + dY(0) + vin;
	dY(2) = (dY(0)/V_B) * (ta - T_B) + (vin/V_B) * (tin - T_B);

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
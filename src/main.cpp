/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab0();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dok�adne rozwi�zanie optymalne
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Wahadlo
	Nmax = 1000;											// dok�adno��
	epsilon = 1e-2;											// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;											// dolne oraz g�rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("/data/symulacja_lab0.csv");					// definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie�
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
}

void lab1() // Daje zły wynik
{
       // Parametry testu
        double a = 0;            // początek przedziału
        double b = 2;            // koniec przedziału
        double epsilon = 1e-4;   // dokładność
        matrix ud1(1, 1), ud2(1, 1); // dane dodatkowe dla funkcji celu, jeśli potrzebne

        // Wywołanie algorytmu Fibonacciego
        solution opt = fib(ff1T, a, b, epsilon, ud1, ud2);

        // Wypisanie wyniku
        cout << "Rozwiazanie znalezione przez Fibonacci: " << endl;
        cout << opt << endl;
        solution::clear_calls(); // zerowanie licznika wywołań 

}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}

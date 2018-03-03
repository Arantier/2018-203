#include "scherbakovdv.h"

/**
 * Введение в дисциплину
 */
void scherbakovdv::lab1()
{
	printf("Hello world");
}


/**
 * Метод Гаусса с выбором главного элемента - done and passed
 */
void scherbakovdv::lab2()
{
	//A[i][i] - матрица, b[i] - столбец свободных членов
	for (int i = 0; i < N; i++)
		{
			if (!A[i][i])
			{
				double* MaxRow = A[i];
				for (int j = i + 1; j < N; j++)
					if (A[j][i] > MaxRow[i])
					{
						double* tmp = MaxRow;
						MaxRow = A[j];
						A[j] = tmp;
					}
				if (!MaxRow[i])
					return;
			}
			for (int j = N - 1; j > i; j--)
			{
				A[i][j] /= A[i][i];
				for (int k = i + 1; k < N; k++)
					A[k][j] -= A[k][i] * A[i][j];
			}
			b[i] /= A[i][i];
			A[i][i] = 1;
			for (int j = N - 1; j > i; j--)
			{
				b[j] -= b[i] * A[j][i];
				A[j][i] = 0;
			}
		}
		for (int i = N - 1; i >= 0; i--)
		{
			for (int j = i - 1; j >= 0; j--)
			{
				b[j] -= b[i] * A[j][i];
				A[j][i] = 0;
			}
		}
		memcpy(x,b,sizeof(double)*N);
}



/**
 * Метод прогонки - done and passed
 */
void scherbakovdv::lab3()
{
	double *al = new double[N], *bt = new double[N];
	//Прямой ход
	al[0] = -A[0][1]/A[0][0];
	bt[0] = b[0]/A[0][0];
	for (int i=1;i<N-1;i++){
		al[i] = -A[i][i+1]/(A[i][i] + A[i][i-1]*al[i-1]);
		bt[i] = (b[i] - A[i][i-1]*bt[i-1])/(A[i][i] + A[i][i-1]*al[i-1]);
	}
	x[N-1] = (b[N-1]-A[N-1][N-2]*bt[N-1])/(A[N-1][N-1]+A[N-1][N-2]*al[N-2]);
	//Обратный ход
	for (int i=N-2;i!=-1;i--)
		x[i] = al[i]*x[i+1]+bt[i];
	delete[] al;
	delete[] bt;
}



/**
 * Метод простых итераций - working
 */
void scherbakovdv::lab4()
{
	//TODO: узнать у жалнина про разницу между методами якоби и простой итерации
	double diff;
	//Допустимая погрешность
	const double Eps=0.1E-10;
	//Аварийный счётчик
	double counter=0;
	//Переходной массив иксов
	double* xOld = new double[N];
	for (int i=0;i<N;i++)
		xOld[i]=0;
	do {
		for (int i=0;i<N;i++) {
			x[i]=b[i];
			for (int j=0;j<N;j++)
				if (i!=j)
					x[i]-=A[i][j]*xOld[j];
			x[i]/=A[i][i];
		}
		diff=fabs(x[0]-xOld[0]);
		counter++;
		memcpy(xOld,x,sizeof(double)*N);
	} while ((diff>Eps)&&(counter<100));
	delete[] xOld;
	if (diff>Eps) {
		perror("diff>Eps");
	} else if (counter==100) {
		perror("Out of cycle");
	}
}



/**
 * Метод Якоби или Зейделя - done
 */
void scherbakovdv::lab5()
{
	//Достаточное признак сходимости метода Якоби - диагональное преобладание
	//Средняя разность между новыми и старыми значениями
	double diff;
	//Допустимая погрешность
	const double Eps=0.1E-10;
	//Аварийный счётчик
	double counter=0;
	//Переходной массив иксов
	double* xOld = new double[N];
	memcpy(xOld,b,sizeof(double)*N);
	do {
		for (int i=0;i<N;i++) {
			x[i]=b[i];
			for (int j=0;j<N;j++)
				if (i!=j)
					x[i]-=A[i][j]*xOld[j];
			x[i]/=A[i][i];
		}
		diff=fabs(x[0]-xOld[0]);
		counter++;
		memcpy(xOld,x,sizeof(double)*N);
	} while ((diff>Eps)&&(counter<100));
	delete[] xOld;
	if (diff>Eps) {
		perror("diff>Eps");
	} else if (counter==100) {
		perror("Out of cycle");
	}
}



/**
 * Метод минимальных невязок - done
 */
void scherbakovdv::lab6()
{
	//Допустимая погрешность
	const double Eps=0.1E-10;
	printf("LAB 6: Eps is %.2e",Eps);
	//Аварийный счётчик
	double counter=0;
	//vec - вектор невязок (интересно, почему так называется)
	double *vec = new double[N], diff=0;
	//начальное приближение будет b
	memcpy(x,b,sizeof(double)*N);
	do {
		memcpy(vec,b,sizeof(double)*N);
		for (int i=0;i<N;i++)
			for (int j=0;j<N;j++)
				vec[i]-=A[i][j]*x[j];
		//Поиск параметра тау в формуле X^(k+1)=X^(k)-tau*vec
		//tau находится по формуле tau=-num/den
		//num=(A*vec,vec) <- скалярное произведение
		//den=(A8vec,A*vec)
		//Из-за обилия A*vec я счёл разумным добавить вектор Avec
		double tau=0, num=0, den=0, *Avec = new double[N];
		for (int i=0;i<N;i++)
			for (int j=0;j<N;j++)
				Avec[i]=A[i][j]*vec[j];
		for (int i=0;i<N;i++)
		{
			num+=Avec[i]*vec[i];
			//Забавно: я забыл поставить "+=" и поставил просто "=", но результат всё равно сошёлся
			den+=Avec[i]*Avec[i];
		}
		// Прощай, овёс
		delete[] Avec;
		tau=-num/den;
		for (int i=0;i<N;i++)
		{
			x[i]-=tau*vec[i];
			//diff=x[i]-xOld[i], но так экономнее
			diff+=fabs(tau*vec[i]);
		}
		//Для чистоты эксперимента будем проверять среднее арифметическое разности
		diff/=N;
		counter++;
	} while ((diff>Eps)&&(counter<1000));
	delete[] vec;
}



/**
 * Метод сопряженных градиентов - working
 */
void scherbakovdv::lab7()
{

}



/**
 * Метод вращения для нахождения собственных значений матрицы
 */
void scherbakovdv::lab8()
{

}

void scherbakovdv::lab9() {}


std::string scherbakovdv::get_name()
{
  return "Scherbakov D.V.";
}

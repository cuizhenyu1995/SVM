#pragma once
#include "svm_run.h"
#include "svm_in.h"
#include "SystemVue/ModelBuilder.h"
#include "SystemVue/CircularBuffer.h"
#include "SystemVue/Matrix.h"
#include "SystemVue/MatrixCircularBuffer.h"
#include "SystemVue/DFParam.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

using namespace std;
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>

using std::sort;
using std::fabs;


double **x;
double *y;
double *alpha;
double *w;
double b;
double c;
double eps = 1e-6;
struct _E{
	double val;
	int index;
}*E;

bool cmp(const _E & a, const _E & b)
{
	return a.val < b.val;
}

int num_dimension;
int num_samples;

double max(double a, double b)
{
	return a>b ? a : b;
}

double min(double a, double b)
{
	return a>b ? b : a;
}

double kernal(double x1[], double x2[], double dimension)
{
	double ans = 0;
	for (int i = 0; i < dimension; i++)
	{
		ans += x1[i] * x2[i];
	}
	return ans;
}

double target_function()
{
	double ans = 0;
	for (int i = 0; i < num_samples; i++)
	{
		for (int j = 0; j < num_samples; j++)
		{
			ans += alpha[i] * alpha[j] * y[i] * y[j] * kernal(x[i], x[j], num_dimension);
		}
	}

	for (int i = 0; i < num_samples; i++)
	{
		ans -= alpha[i];
	}

	return ans;
}


double g(double _x[], int dimension)
{
	double ans = b;

	for (int i = 0; i < num_samples; i++)
	{
		ans += alpha[i] * y[i] * kernal(x[i], _x, dimension);
	}

	return ans;
}

bool satisfy_constrains(int i, int dimension)
{
	if (alpha[i] == 0)
	{
		if (y[i] * g(x[i], dimension) >= 1)
			return true;
		else
			return false;
	}
	else if (alpha[i] > 0 && alpha[i] < c)
	{
		if (y[i] * g(x[i], dimension) == 1)
			return true;
		else
			return false;
	}
	else
	{
		if (y[i] * g(x[i], dimension) <= 1)
			return true;
		else
			return false;
	}
}


double calE(int i, int dimension)
{
	return g(x[i], dimension) - y[i];
}

void calW()
{
	for (int i = 0; i < num_dimension; i++)
	{
		w[i] = 0;
		for (int j = 0; j < num_samples; j++)
		{
			w[i] += alpha[j] * y[j] * x[j][i];
		}
	}
	return;
}

void calB()
{
	double ans = y[0];
	for (int i = 0; i < num_samples; i++)
	{
		ans -= y[i] * alpha[i] * kernal(x[i], x[0], num_dimension);
	}
	b = ans;
	return;
}


void recalB(int alpha1index, int alpha2index, int dimension, double alpha1old, double alpha2old)
{
	double alpha1new = alpha[alpha1index];
	double alpha2new = alpha[alpha2index];

	alpha[alpha1index] = alpha1old;
	alpha[alpha2index] = alpha2old;

	double e1 = calE(alpha1index, num_dimension);
	double e2 = calE(alpha2index, num_dimension);

	alpha[alpha1index] = alpha1new;
	alpha[alpha2index] = alpha2new;

	double b1new = -e1 - y[alpha1index] * kernal(x[alpha1index], x[alpha1index], dimension)*(alpha1new - alpha1old);
	b1new -= y[alpha2index] * kernal(x[alpha2index], x[alpha1index], dimension)*(alpha2new - alpha2old) + b;

	double b2new = -e2 - y[alpha1index] * kernal(x[alpha1index], x[alpha2index], dimension)*(alpha1new - alpha1old);
	b1new -= y[alpha2index] * kernal(x[alpha2index], x[alpha2index], dimension)*(alpha2new - alpha2old) + b;

	b = (b1new + b2new) / 2;
}

bool optimizehelp(int alpha1index, int alpha2index)
{
	double alpha1new = alpha[alpha1index];
	double alpha2new = alpha[alpha2index];

	double alpha1old = alpha[alpha1index];
	double alpha2old = alpha[alpha2index];

	double H, L;

	if (fabs(y[alpha1index] - y[alpha2index]) > eps)
	{
		L = max(0, alpha2old - alpha1old);
		H = min(c, c + alpha2old - alpha1old);
	}
	else
	{
		L = max(0, alpha2old + alpha1old - c);
		H = min(c, alpha2old + alpha1old);
	}

	//cal new
	double lena = kernal(x[alpha1index], x[alpha1index], num_dimension) + kernal(x[alpha2index], x[alpha2index], num_dimension) - 2 * kernal(x[alpha1index], x[alpha2index], num_dimension);
	alpha2new = alpha2old + y[alpha2index] * (calE(alpha1index, num_dimension) - calE(alpha2index, num_dimension)) / lena;

	if (alpha2new > H)
	{
		alpha2new = H;
	}
	else if (alpha2new < L)
	{
		alpha2new = L;
	}

	alpha1new = alpha1old + y[alpha1index] * y[alpha2index] * (alpha2old - alpha2new);

	double energyold = target_function();

	alpha[alpha1index] = alpha1new;
	alpha[alpha2index] = alpha2new;

	double gap = 0.001;

	recalB(alpha1index, alpha2index, num_dimension, alpha1old, alpha2old);
	return true;
}

bool optimize()
{
	int alpha1index = -1;
	int alpha2index = -1;
	double alpha2new = 0;
	double alpha1new = 0;

	//cal E[]
	for (int i = 0; i < num_samples; i++)
	{
		E[i].val = calE(i, num_dimension);
		E[i].index = i;
	}

	//traverse the alpha1index with 0 < && < c
	for (int i = 0; i < num_samples; i++)
	{
		alpha1new = alpha[i];

		if (alpha1new > 0 && alpha1new < c)
		{

			if (satisfy_constrains(i, num_dimension))
				continue;

			sort(E, E + num_samples, cmp);

			//simply find the maximum or minimun;
			if (alpha1new > 0)
			{
				if (E[0].index == i)
				{
					;
				}
				else
				{
					alpha1index = i;
					alpha2index = E[0].index;
					if (optimizehelp(alpha1index, alpha2index))
					{
						return true;
					}
				}
			}
			else
			{
				if (E[num_samples - 1].index == i)
				{
					;
				}
				else
				{
					alpha1index = i;
					alpha2index = E[num_samples - 1].index;
					if (optimizehelp(alpha1index, alpha2index))
					{
						return true;
					}
				}
			}


			//find the alpha2 > 0 && < c
			for (int j = 0; j < num_samples; j++)
			{
				alpha2new = alpha[j];

				if (alpha2new > 0 && alpha2new < c)
				{
					alpha1index = i;
					alpha2index = j;
					if (optimizehelp(alpha1index, alpha2index))
					{
						return true;
					}
				}
			}

			//find other alpha2
			for (int j = 0; j < num_samples; j++)
			{
				alpha2new = alpha[j];

				if (!(alpha2new > 0 && alpha2new < c))
				{
					alpha1index = i;
					alpha2index = j;
					if (optimizehelp(alpha1index, alpha2index))
					{
						return true;
					}
				}
			}
		}
	}

	//find all alpha1
	for (int i = 0; i < num_samples; i++)
	{
		alpha1new = alpha[i];

		if (!(alpha1new > 0 && alpha1new < c))
		{
			if (satisfy_constrains(i, num_dimension))
				continue;

			sort(E, E + num_samples, cmp);

			//simply find the maximum or minimun;
			if (alpha1new > 0)
			{
				if (E[0].index == i)
				{
					;
				}
				else
				{
					alpha1index = i;
					alpha2index = E[0].index;
					if (optimizehelp(alpha1index, alpha2index))
					{
						return true;
					}
				}
			}
			else
			{
				if (E[num_samples - 1].index == i)
				{
					;
				}
				else
				{
					alpha1index = i;
					alpha2index = E[num_samples - 1].index;
					if (optimizehelp(alpha1index, alpha2index))
					{
						return true;
					}
				}
			}


			//find the alpha2 > 0 && < c
			for (int j = 0; j < num_samples; j++)
			{
				alpha2new = alpha[j];

				if (alpha2new > 0 && alpha2new < c)
				{
					alpha1index = i;
					alpha2index = j;
					if (optimizehelp(alpha1index, alpha2index))
					{
						return true;
					}
				}
			}

			//find other alpha2
			for (int j = 0; j < num_samples; j++)
			{
				alpha2new = alpha[j];

				if (!(alpha2new > 0 && alpha2new < c))
				{
					alpha1index = i;
					alpha2index = j;
					if (optimizehelp(alpha1index, alpha2index))
					{
						return true;
					}
				}
			}
		}
	}

	//for(int i = 0 ; i < num_samples; i++)
	//{
	//    alpha1new = alpha[i];

	//    for(int j = 0 ; j < num_samples; j++)
	//    {
	//        if(1)
	//        {
	//            alpha1index = i;
	//            alpha2index = j;
	//            if(optimizehelp(alpha1index , alpha2index))
	//            {
	//                return true;
	//            }
	//        }
	//    }
	//}
	return false;
}

bool check()
{
	double sum = 0;
	for (int i = 0; i < num_samples; i++)
	{
		sum += alpha[i] * y[i];
		if (!(0 <= alpha[i] && alpha[i] <= c))
		{
			printf("alpha[%d]: %lf wrong\n", i, alpha[i]);
			return false;
		}
		if (!satisfy_constrains(i, num_dimension))
		{
			printf("alpha[%d] not satisfy constrains\n", i);
			return false;
		}
	}

	if (fabs(sum) > eps)
	{
		printf("Sum = %lf\n", sum);
		return false;
	}
	return true;
}
/*
min 1/2*||w||^2
s.t.  (w[i]*x[i] + b[i] - y[i]) >= 0;
*/
/*
step 1: cal alpha[]
step 2: cal w,b
*/

/*
min(para alpha) 1/2*sum(i)sum(j)(alpha[i]*alpha[j]*y[i]*y[j]*x[i]*x[j]) - sum(alpha[i])
s.t. sum(alpha[i] * y[i]) = 0
C>= alpha[i] >= 0
*/

namespace SystemVueModelBuilder {

#ifndef SV_CODE_MRA_RUN
	DEFINE_MODEL_INTERFACE(svm_run)
	{
		SET_MODEL_NAME("svm_run");
		SET_MODEL_DESCRIPTION("");
		SET_MODEL_CATEGORY("svm");
		ADD_MODEL_HEADER_FILE("svm_run.h");
		model.SetModelCodeGenName("svm_run");
		model.SetModelNamespace("SystemVueModelBuilder");

		// Add parameters
		SystemVueModelBuilder::DFParam param = NULL;


		// Add input/output ports
		DFPort IN_N = ADD_MODEL_INPUT(IN_ELEMENT);
		IN_N.SetDescription("Input the number of elements");
		DFPort IN_M = ADD_MODEL_INPUT(IN_ATTRIBUTE);
		IN_M.SetDescription("Input the number of attributes");
		DFPort IN_FILE = ADD_MODEL_INPUT(IN_INPUTFILE);
		IN_FILE.SetDescription("Input the data");


		DFPort OutputPort = ADD_MODEL_OUTPUT(Output);
		OutputPort.SetDescription("Output");

		return true;
	}
#endif

	svm_run::svm_run()
	{
	}

	svm_run::~svm_run()
	{
	}

	bool svm_run::Setup()
	{
		IN_ELEMENT.SetRate(1);
		IN_ATTRIBUTE.SetRate(1);
		IN_INPUTFILE.SetRate(1);
		Output.SetRate(1);

		return true;
	};

	bool svm_run::Initialize()
	{
		return true;
	}

	//bool svm_run::cmp(const _E & a, const _E & b)
	//{
	//	return a.val < b.val;
	//}
	//double svm_run::max(double a, double b)
	//{
	//	return a>b ? a : b;
	//}
	//double svm_run::min(double a, double b)
	//{
	//	return a>b ? b : a;
	//}
	//double svm_run::kernal(double x1[], double x2[], double dimension)
	//{
	//	double ans = 0;
	//	for (int i = 0; i < dimension; i++)
	//	{
	//		ans += x1[i] * x2[i];
	//	}
	//	return ans;
	//}
	//double svm_run::target_function()
	//{
	//	double ans = 0;
	//	for (int i = 0; i < num_samples; i++)
	//	{
	//		for (int j = 0; j < num_samples; j++)
	//		{
	//			ans += alpha[i] * alpha[j] * y[i] * y[j] * kernal(x[i], x[j], num_dimension);
	//		}
	//	}

	//	for (int i = 0; i < num_samples; i++)
	//	{
	//		ans -= alpha[i];
	//	}

	//	return ans;
	//}
	//double svm_run::g(double _x[], int dimension)
	//{
	//	double ans = b;

	//	for (int i = 0; i < num_samples; i++)
	//	{
	//		ans += alpha[i] * y[i] * kernal(x[i], _x, dimension);
	//	}

	//	return ans;
	//}
	//bool svm_run::satisfy_constrains(int i, int dimension)
	//{
	//	if (alpha[i] == 0)
	//	{
	//		if (y[i] * g(x[i], dimension) >= 1)
	//			return true;
	//		else
	//			return false;
	//	}
	//	else if (alpha[i] > 0 && alpha[i] < c)
	//	{
	//		if (y[i] * g(x[i], dimension) == 1)
	//			return true;
	//		else
	//			return false;
	//	}
	//	else
	//	{
	//		if (y[i] * g(x[i], dimension) <= 1)
	//			return true;
	//		else
	//			return false;
	//	}
	//}
	//double svm_run::calE(int i, int dimension)
	//{
	//	return g(x[i], dimension) - y[i];
	//}
	//void svm_run::calW()
	//{
	//	for (int i = 0; i < num_dimension; i++)
	//	{
	//		w[i] = 0;
	//		for (int j = 0; j < num_samples; j++)
	//		{
	//			w[i] += alpha[j] * y[j] * x[j][i];
	//		}
	//	}
	//	return;
	//}
	//void svm_run::calB()
	//{
	//	double ans = y[0];
	//	for (int i = 0; i < num_samples; i++)
	//	{
	//		ans -= y[i] * alpha[i] * kernal(x[i], x[0], num_dimension);
	//	}
	//	b = ans;
	//	return;
	//}
	//void svm_run::recalB(int alpha1index, int alpha2index, int dimension, double alpha1old, double alpha2old)
	//{
	//	double alpha1new = alpha[alpha1index];
	//	double alpha2new = alpha[alpha2index];

	//	alpha[alpha1index] = alpha1old;
	//	alpha[alpha2index] = alpha2old;

	//	double e1 = calE(alpha1index, num_dimension);
	//	double e2 = calE(alpha2index, num_dimension);

	//	alpha[alpha1index] = alpha1new;
	//	alpha[alpha2index] = alpha2new;

	//	double b1new = -e1 - y[alpha1index] * kernal(x[alpha1index], x[alpha1index], dimension)*(alpha1new - alpha1old);
	//	b1new -= y[alpha2index] * kernal(x[alpha2index], x[alpha1index], dimension)*(alpha2new - alpha2old) + b;

	//	double b2new = -e2 - y[alpha1index] * kernal(x[alpha1index], x[alpha2index], dimension)*(alpha1new - alpha1old);
	//	b1new -= y[alpha2index] * kernal(x[alpha2index], x[alpha2index], dimension)*(alpha2new - alpha2old) + b;

	//	b = (b1new + b2new) / 2;
	//}
	//bool svm_run::optimizehelp(int alpha1index, int alpha2index)
	//{
	//	double alpha1new = alpha[alpha1index];
	//	double alpha2new = alpha[alpha2index];

	//	double alpha1old = alpha[alpha1index];
	//	double alpha2old = alpha[alpha2index];

	//	double H, L;

	//	if (fabs(y[alpha1index] - y[alpha2index]) > eps)
	//	{
	//		L = max(0, alpha2old - alpha1old);
	//		H = min(c, c + alpha2old - alpha1old);
	//	}
	//	else
	//	{
	//		L = max(0, alpha2old + alpha1old - c);
	//		H = min(c, alpha2old + alpha1old);
	//	}

	//	//cal new
	//	double lena = kernal(x[alpha1index], x[alpha1index], num_dimension) + kernal(x[alpha2index], x[alpha2index], num_dimension) - 2 * kernal(x[alpha1index], x[alpha2index], num_dimension);
	//	alpha2new = alpha2old + y[alpha2index] * (calE(alpha1index, num_dimension) - calE(alpha2index, num_dimension)) / lena;

	//	if (alpha2new > H)
	//	{
	//		alpha2new = H;
	//	}
	//	else if (alpha2new < L)
	//	{
	//		alpha2new = L;
	//	}

	//	alpha1new = alpha1old + y[alpha1index] * y[alpha2index] * (alpha2old - alpha2new);

	//	double energyold = target_function();

	//	alpha[alpha1index] = alpha1new;
	//	alpha[alpha2index] = alpha2new;

	//	double gap = 0.001;

	//	recalB(alpha1index, alpha2index, num_dimension, alpha1old, alpha2old);
	//	return true;
	//}
	//bool svm_run::optimize()
	//{
	//	int alpha1index = -1;
	//	int alpha2index = -1;
	//	double alpha2new = 0;
	//	double alpha1new = 0;

	//	//cal E[]
	//	for (int i = 0; i < num_samples; i++)
	//	{
	//		E[i].val = calE(i, num_dimension);
	//		E[i].index = i;
	//	}

	//	//traverse the alpha1index with 0 < && < c
	//	for (int i = 0; i < num_samples; i++)
	//	{
	//		alpha1new = alpha[i];

	//		if (alpha1new > 0 && alpha1new < c)
	//		{

	//			if (satisfy_constrains(i, num_dimension))
	//				continue;

	//			sort(E, E + num_samples, cmp);

	//			//simply find the maximum or minimun;
	//			if (alpha1new > 0)
	//			{
	//				if (E[0].index == i)
	//				{
	//					;
	//				}
	//				else
	//				{
	//					alpha1index = i;
	//					alpha2index = E[0].index;
	//					if (optimizehelp(alpha1index, alpha2index))
	//					{
	//						return true;
	//					}
	//				}
	//			}
	//			else
	//			{
	//				if (E[num_samples - 1].index == i)
	//				{
	//					;
	//				}
	//				else
	//				{
	//					alpha1index = i;
	//					alpha2index = E[num_samples - 1].index;
	//					if (optimizehelp(alpha1index, alpha2index))
	//					{
	//						return true;
	//					}
	//				}
	//			}


	//			//find the alpha2 > 0 && < c
	//			for (int j = 0; j < num_samples; j++)
	//			{
	//				alpha2new = alpha[j];

	//				if (alpha2new > 0 && alpha2new < c)
	//				{
	//					alpha1index = i;
	//					alpha2index = j;
	//					if (optimizehelp(alpha1index, alpha2index))
	//					{
	//						return true;
	//					}
	//				}
	//			}

	//			//find other alpha2
	//			for (int j = 0; j < num_samples; j++)
	//			{
	//				alpha2new = alpha[j];

	//				if (!(alpha2new > 0 && alpha2new < c))
	//				{
	//					alpha1index = i;
	//					alpha2index = j;
	//					if (optimizehelp(alpha1index, alpha2index))
	//					{
	//						return true;
	//					}
	//				}
	//			}
	//		}
	//	}

	//	//find all alpha1
	//	for (int i = 0; i < num_samples; i++)
	//	{
	//		alpha1new = alpha[i];

	//		if (!(alpha1new > 0 && alpha1new < c))
	//		{
	//			if (satisfy_constrains(i, num_dimension))
	//				continue;

	//			sort(E, E + num_samples, cmp);

	//			//simply find the maximum or minimun;
	//			if (alpha1new > 0)
	//			{
	//				if (E[0].index == i)
	//				{
	//					;
	//				}
	//				else
	//				{
	//					alpha1index = i;
	//					alpha2index = E[0].index;
	//					if (optimizehelp(alpha1index, alpha2index))
	//					{
	//						return true;
	//					}
	//				}
	//			}
	//			else
	//			{
	//				if (E[num_samples - 1].index == i)
	//				{
	//					;
	//				}
	//				else
	//				{
	//					alpha1index = i;
	//					alpha2index = E[num_samples - 1].index;
	//					if (optimizehelp(alpha1index, alpha2index))
	//					{
	//						return true;
	//					}
	//				}
	//			}


	//			//find the alpha2 > 0 && < c
	//			for (int j = 0; j < num_samples; j++)
	//			{
	//				alpha2new = alpha[j];

	//				if (alpha2new > 0 && alpha2new < c)
	//				{
	//					alpha1index = i;
	//					alpha2index = j;
	//					if (optimizehelp(alpha1index, alpha2index))
	//					{
	//						return true;
	//					}
	//				}
	//			}

	//			//find other alpha2
	//			for (int j = 0; j < num_samples; j++)
	//			{
	//				alpha2new = alpha[j];

	//				if (!(alpha2new > 0 && alpha2new < c))
	//				{
	//					alpha1index = i;
	//					alpha2index = j;
	//					if (optimizehelp(alpha1index, alpha2index))
	//					{
	//						return true;
	//					}
	//				}
	//			}
	//		}
	//	}

	//	//for(int i = 0 ; i < num_samples; i++)
	//	//{
	//	//    alpha1new = alpha[i];

	//	//    for(int j = 0 ; j < num_samples; j++)
	//	//    {
	//	//        if(1)
	//	//        {
	//	//            alpha1index = i;
	//	//            alpha2index = j;
	//	//            if(optimizehelp(alpha1index , alpha2index))
	//	//            {
	//	//                return true;
	//	//            }
	//	//        }
	//	//    }
	//	//}
	//	return false;
	//}
	//bool svm_run::check()
	//{
	//	double sum = 0;
	//	for (int i = 0; i < num_samples; i++)
	//	{
	//		sum += alpha[i] * y[i];
	//		if (!(0 <= alpha[i] && alpha[i] <= c))
	//		{
	//			printf("alpha[%d]: %lf wrong\n", i, alpha[i]);
	//			return false;
	//		}
	//		if (!satisfy_constrains(i, num_dimension))
	//		{
	//			printf("alpha[%d] not satisfy constrains\n", i);
	//			return false;
	//		}
	//	}

	//	if (fabs(sum) > eps)
	//	{
	//		printf("Sum = %lf\n", sum);
	//		return false;
	//	}
	//	return true;
	//}
	bool svm_run::Run()
	{
		num_samples = IN_ELEMENT[0];
		num_dimension = IN_ATTRIBUTE[0];
		x = new double *[num_samples + 5];
		for (int i = 0; i<num_samples + 5; ++i)
			x[i] = new double[num_dimension + 5];
		y = new double[num_samples + 5];
		alpha = new double[num_samples + 5];
		w = new double[num_dimension + 5];

		for (int i = 0; i < num_samples; i++)
		{
			for (int j = 0; j < num_dimension; j++)
			{
				x[i][j]=IN_INPUTFILE[0](i,j);
			}
			y[i] = IN_INPUTFILE[0](i, num_dimension);
		}
		c = 1;
		for (int i = 0; i < num_samples; i++)
		{
			alpha[i] = 0;
		}
		int count = 0;
		while (optimize()){
			calB();
			count++;
		}
		/*printf("%d ", count);*/

		calW();
		calB();

		/*printf("y = ");*/
		out_tmp.Resize(1,NUMBER_OF_ATTRIBUTE + 1);
		out_tmp.Zero();

		for (int i = 0; i < num_dimension; i++)
		{
			out_tmp(0,i)=w[i];
		}
		out_tmp(0, num_dimension) = b;
		Output[0] = out_tmp;
		/*if (!check())
			printf("Not satisfy KKT.\n");
		else
			printf("Satisfy KKT\n");*/
		for (int i = 0; i<num_samples + 5; ++i)
			delete[] x[i];
		delete[] x;
		delete[] y;
		delete[] alpha;
		delete[] w;
		return true;
	}

	bool svm_run::Finalize()
	{
		return true;
	}


}
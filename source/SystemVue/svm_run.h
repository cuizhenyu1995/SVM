#pragma once
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

namespace SystemVueModelBuilder
{
	class svm_run : public SystemVueModelBuilder::DFModel
	{
	public:
		DECLARE_MODEL_INTERFACE(svm_run);
		svm_run();
		virtual ~svm_run();

		virtual bool Setup();
		virtual bool Initialize();
		virtual bool Run();
		virtual bool Finalize();

		/*struct _E{
			double val;
			int index;
		};
		bool cmp(const _E & a, const _E & b);
		double max(double a, double b);
		double min(double a, double b);
		double kernal(double x1[], double x2[], double dimension);
		double target_function();
		double g(double _x[], int dimension);
		bool satisfy_constrains(int i, int dimension);
		double calE(int i, int dimension);
		void calW();
		void calB();
		void recalB(int alpha1index, int alpha2index, int dimension, double alpha1old, double alpha2old);
		bool optimizehelp(int alpha1index, int alpha2index);
		bool optimize();
		bool check();*/


		/// define parameters
		int NUMBER_OF_ELEMENT;
		int NUMBER_OF_ATTRIBUTE;

		/// define ports
		IntCircularBuffer IN_ELEMENT;
		IntCircularBuffer IN_ATTRIBUTE;
		DoubleMatrixCircularBuffer IN_INPUTFILE;

		DoubleMatrixCircularBuffer Output;

	protected:
		DoubleMatrix in_tmp;
		DoubleMatrix out_tmp;

	};
}
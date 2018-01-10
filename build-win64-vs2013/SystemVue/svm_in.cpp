#include "svm_in.h"
#include "svm_run.h"
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

namespace SystemVueModelBuilder {

#ifndef SV_CODE_MRA_IN
	DEFINE_MODEL_INTERFACE(svm_in)
	{
		SET_MODEL_NAME("svm_in");
		SET_MODEL_DESCRIPTION("");
		SET_MODEL_CATEGORY("svm_in");
		ADD_MODEL_HEADER_FILE("svm_in.h");
		model.SetModelCodeGenName("svm_in");
		model.SetModelNamespace("SystemVueModelBuilder");

		// Add parameters
		SystemVueModelBuilder::DFParam param = NULL;
		param = ADD_MODEL_PARAM(NUMBER_OF_ELEMENT);
		param.SetDefaultValue(0);
		param = ADD_MODEL_PARAM(NUMBER_OF_ATTRIBUTE);
		param.SetDefaultValue(0);
		param = ADD_MODEL_PARAM(INPUT_FILE);
		param.SetParamAsFile();

		// Add input/output ports
		DFPort OUT_N = ADD_MODEL_OUTPUT(OUT_ELEMENT);
		OUT_N.SetDescription("Output the number of elements");
		DFPort OUT_M = ADD_MODEL_OUTPUT(OUT_ATTRIBUTE);
		OUT_M.SetDescription("Output the number of attributes");
		DFPort OUT_FILE = ADD_MODEL_OUTPUT(OUT_INPUTFILE);
		OUT_FILE.SetDescription("Output the data");

		return true;
	}
#endif

	svm_in::svm_in()
	{
	}

	svm_in::~svm_in()
	{
	}

	bool svm_in::Setup()
	{
		//setRate的设置针对的是每一个端口，输入端口和输出端口都要进行设置，不设置的话默认为1
		//端口可以认为是一个Buffer, 下游的输入端口会拷贝上游的输出端口的内容，该拷贝为值拷贝或者说是深拷贝
		OUT_ELEMENT.SetRate(1);
		OUT_ATTRIBUTE.SetRate(1);
		OUT_INPUTFILE.SetRate(1);

		return true;
	};

	bool svm_in::Initialize()
	{
		return true;
	}

	//Run函数的终极目的就是产生数据，并赋值到相应的端口
	bool svm_in::Run()
	{
		OUT_ELEMENT[0] = NUMBER_OF_ELEMENT;
		OUT_ATTRIBUTE[0] = NUMBER_OF_ATTRIBUTE;
		OUT_INPUT.Resize(NUMBER_OF_ELEMENT, NUMBER_OF_ATTRIBUTE + 1);
		OUT_INPUT.Zero();
		FILE* in = fopen(INPUT_FILE, "r");

		for (int i = 0; i < NUMBER_OF_ELEMENT; i++)
		{
			for (int j = 0; j < NUMBER_OF_ATTRIBUTE + 1; j++)
			{
				fscanf(in, "%lf", &OUT_INPUT(i, j));
			}
		}
		OUT_INPUTFILE[0] = OUT_INPUT;
		fclose(in);
		return true;
	}


	bool svm_in::Finalize() {

		return true;
	}

}
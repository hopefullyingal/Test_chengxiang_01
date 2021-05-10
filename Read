#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <graphics.h>
#include <conio.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <istream>
#include<string>

using namespace std;
class Read
{


	int i, j, k, s, t;
	double ElecThickness;//电极厚度
	int ElecNumber;//电极总数
	double* DisBetweenAdjacentEle;//相邻电极之间的距离
	double* StepBetweenAdjacentEle;//相邻电极划分的步长
	double* ElecPotential;//电极电位
	double RadiusAperture;//电极内孔径半径
	int GridAperture;//电极内孔径半径等步长划分的网格数
	double RadiusBoundary;//从电极内孔沿到封闭边界处的径向距离
	int GridBoundary;//从电极内孔沿到封闭边界处划分的网格数
	double IterationAccuracy;//迭代精度
	int NST;//输出打印空间电位时网格点间隔数
	double* EOE;//要求扫描搜索等电位线的电位间隔值
	double dov;//要求扫描等位线间隔值
	FILE* ParameterFile;
	FILE* OutputFile;

public:
	struct data{
	public :

		double mElecThickness;//电极厚度
		int mElecNumber;//电极总数
		double* mDisBetweenAdjacentEle;//相邻电极之间的距离
		double* mStepBetweenAdjacentEle;//相邻电极划分的步长
		double* mElecPotential;//电极电位
		double mRadiusAperture;//电极内孔径半径
		int mGridAperture;//电极内孔径半径等步长划分的网格数
		double mRadiusBoundary;//从电极内孔沿到封闭边界处的径向距离
		int mGridBoundary;//从电极内孔沿到封闭边界处划分的网格数
		double mIterationAccuracy;//迭代精度
		int mNST;//输出打印空间电位时网格点间隔数
		double* mEOE;//要求扫描搜索等电位线的电位间隔值
		double mdov;//要求扫描等位线间隔值


	};


	data getValue() {
		/*读入数据*/
		fopen_s(&ParameterFile, "E:\\Learning(School)\\光电成像原理与技术\\光电成像实验\\程序\\1120160852.dat", "r+");//读取相关参数文件
		fopen_s(&OutputFile, "E:\\Learning(School)\\光电成像原理与技术\\光电成像实验\\程序\\11.res", "w+");
		/// <summary>备份
		/// errno_t err = fopen_s(&ParameterFile, "E:\\Learning(School)\\光电成像原理与技术\\光电成像实验\\程序\\1120160852.dat", "r+");//读取相关参数文件
		/// </summary>
		/// <returns></returns>


		if (ParameterFile == NULL)//detect the file 
		{
			printf("Open filefailure!");
			exit(1);
		}


		/*test*/
		////////////////////////////////////////////////////////////////////////
		// using others way to get numbers
		//ifstream in("E:\\Learning(School)\\111.txt");
		//string s;
		//string str[15];

		//int hh = 0;

		//while (getline(in, s))//着行读取数据并存于s中，直至数据全部读取
		//{
		//	 str[hh] = UTF8ToGB(s.c_str()).c_str();
		//	 hh++;
		//	
		//}
		//
			/*TODO: 数组不能开，ElecNumber 不是常量*/
		////////////////////////////////////////////////////////////////
		fscanf_s(ParameterFile, "电极总数：%d；\n", &ElecNumber);
		fscanf_s(ParameterFile, "电极厚度：%lf；\n", &ElecThickness);

		DisBetweenAdjacentEle = (double*)calloc(ElecNumber, sizeof(double));//分配电极距离数组大小
		fscanf_s(ParameterFile, "相邻电极之间的距离：%lf，", &DisBetweenAdjacentEle[0]);//扫描距离数组第一个元素		
		for (i = 1; i < ElecNumber-1; i++)
		{
			fscanf_s(ParameterFile, "%lf，", &DisBetweenAdjacentEle[i]);

		}//把第二个至倒数第二个距离元素都扫描到数组里
		fscanf_s(ParameterFile, "%lf；\n", &DisBetweenAdjacentEle[ElecNumber-1]);//把最后一个距离元素扫描到数组里
		

		StepBetweenAdjacentEle = (double*)calloc(ElecNumber, sizeof(double));//给步长数组分配内存空间
		fscanf_s(ParameterFile, "相邻电极划分的步长：%lf，", &StepBetweenAdjacentEle[0]);//扫描步长数组第一个元素
		for (i = 1; i < ElecNumber-1; i++)
		{
			fscanf_s(ParameterFile, "%lf，", &StepBetweenAdjacentEle[i]);
			
		}//把从第二个至倒数第二个步长元素扫描到数组里
		fscanf_s(ParameterFile, "%lf；\n", &StepBetweenAdjacentEle[ElecNumber-1]);//把最后一个数组元素扫描到数组里



		ElecPotential = (double*)calloc(ElecNumber, sizeof(double));//Allocates memory space to the electronic potential array.
		fscanf_s(ParameterFile, "电极电位：%lf，", &ElecPotential[0]);
		for (i = 1; i < ElecNumber-1; i++)
		{
			fscanf_s(ParameterFile, "%lf，", &ElecPotential[i]);
		}
		fscanf_s(ParameterFile, "%lf；\n", &ElecPotential[ElecNumber-1]);//Put the last element into the array
	
		fscanf_s(ParameterFile, "电极内孔径半径：%lf；\n", &RadiusAperture);//Scan the  radius of the inner aperture of the electrode
		fscanf_s(ParameterFile, "电极内孔径半径等步长划分的网格数：%d；\n", &GridAperture);//扫描电极内孔径半径等步长划分的网格数
		fscanf_s(ParameterFile, "从电极内孔沿到封闭边界处的径向距离：%lf；\n", &RadiusBoundary);//扫描从电极内孔沿到封闭边界处的径向距离
		fscanf_s(ParameterFile, "从电极内孔沿到封闭边界处划分的网格数：%d；\n", &GridBoundary);//扫描从电极内孔沿到封闭边界处划分的网格数
		fscanf_s(ParameterFile, "迭代精度：%lf；\n", &IterationAccuracy);//扫描迭代精度
		fscanf_s(ParameterFile, "输出打印空间电位时网格点间隔数：%d；\n", &NST);//扫描输出打印空间电位时网格点间隔数
		fscanf_s(ParameterFile, "要求扫描搜索等电位线的电位间隔值或电位值:%lf；", &dov);//扫描轴上电位作等距插值时的步长数：



		//TODO EOE 与dov 重复		
		EOE = (double*)calloc(100, sizeof(double));
		//dov = 0;
		i = -1;
		do
		{
			i++;
			fscanf_s(ParameterFile, "%lf；", &EOE[i]);

		} while (EOE[i] != 0);
		if (i == 1)
			dov = EOE[0];

		fclose(ParameterFile);

		
		
		
		data Mdata;
		Mdata.mdov = dov;
		Mdata.mIterationAccuracy = IterationAccuracy;
		Mdata.mEOE = EOE;
		Mdata.mElecThickness = ElecThickness;
		Mdata.mGridAperture = GridAperture;
		Mdata.mGridBoundary = GridBoundary;
		Mdata.mStepBetweenAdjacentEle = StepBetweenAdjacentEle;
		Mdata.mElecNumber = ElecNumber;
		Mdata.mNST = NST;
		Mdata.mRadiusAperture = RadiusAperture;
		Mdata.mRadiusBoundary = RadiusBoundary;
		Mdata.mElecPotential = ElecPotential;
		Mdata.mDisBetweenAdjacentEle = DisBetweenAdjacentEle;


		//setfree();
		return Mdata;




	
	}

	void show(data Mdata ) {

		std::cout << "电极厚度：" << Mdata.mElecThickness << std::endl;
		std::cout << "电极总数：" << Mdata.mElecNumber << std::endl;
		std::cout << "相邻电极之间的距离：" << std::endl;
		for (int i = 0; i < ElecNumber; i++) {
			std::cout<< Mdata.mDisBetweenAdjacentEle[i] << std::endl;
		}


		std::cout << "相邻电极划分的步长：" << std::endl;
		for (int i = 0; i < ElecNumber; i++) {
			std::cout << Mdata.mStepBetweenAdjacentEle[i] << std::endl;
		}
		
		std::cout << "电极电位：" << std::endl;
		for (int i = 0; i < ElecNumber; i++) {
			std::cout << Mdata.mElecPotential[i] << std::endl;
		}

	
		std::cout << "电极内孔径半径：" << Mdata.mRadiusAperture << std::endl;
		std::cout << "电极内孔径半径等步长划分的网格数：" << Mdata.mGridAperture << std::endl;
		std::cout << "从电极内孔沿到封闭边界处的径向距离：" << Mdata.mRadiusBoundary << std::endl;
		std::cout << "从电极内孔沿到封闭边界处划分的网格数:" << Mdata.mGridBoundary << std::endl;
		std::cout << "迭代精度:" << Mdata.mIterationAccuracy << std::endl;
		std::cout << "输出打印空间电位时网格点间隔数:" << Mdata.mNST << std::endl;
		std::cout << "要求扫描等位线间隔值：" << Mdata.mdov << std::endl;

		std::cout << "要求扫描搜索等电位线的电位间隔值:" << std::endl;
		for (int i = 0; i < ElecNumber; i++) {
			std::cout << Mdata.mEOE[i] << std::endl;
		}

		





	}


	/*这个函数是用来测试时输出不乱码*/
	std::string UTF8ToGB(const char* str)
	{
		std::string result;
		WCHAR* strSrc;
		LPSTR szRes;

		//获得临时变量的大小
		int i = MultiByteToWideChar(CP_UTF8, 0, str, -1, NULL, 0);
		strSrc = new WCHAR[i + 1];
		MultiByteToWideChar(CP_UTF8, 0, str, -1, strSrc, i);

		//获得临时变量的大小
		i = WideCharToMultiByte(CP_ACP, 0, strSrc, -1, NULL, 0, NULL, NULL);
		szRes = new CHAR[i + 1];
		WideCharToMultiByte(CP_ACP, 0, strSrc, -1, szRes, i, NULL, NULL);

		result = szRes;
		delete[]strSrc;
		delete[]szRes;

		return result;
	}


	void setfree() {
		free(DisBetweenAdjacentEle);
		free(StepBetweenAdjacentEle);
		free(ElecPotential);
		free(EOE);
	}

};


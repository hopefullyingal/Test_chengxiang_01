#pragma once
#include "Read.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <graphics.h>
#include <conio.h>
#include <iostream>
#include <cstdio>


class inital
{



	double* E;//电场分布
	int row;
	int lin;//行，列
	int* PointGridElectrode;//电极位置对应格点坐标
	double* Z_StepGrid;//z轴方向网络步长划分
	double* Z_PointGrid;//z轴格点位置坐标
	double* Z_PointGridElec;//z轴等电位格点位置坐标
	double* R_StepGrid;//r轴方向网络步长划分
	double* R_PointGrid;//r轴格点位置坐标
	double* TransR_PointGrid;//转换后r轴格点位置坐标
	int* NumEquiPotential ;//等位点个数
	double* R_EquiPotentialPointGrid;//r轴等电位格点位置坐标
	double* TransR_EquiPotentialPointGrid;//转换后r轴等电位格点位置坐标


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

	int i, j, k,s,t;
	


public: void setValue(Read::data Mdata) {
		ElecThickness = Mdata.mElecThickness;
		ElecNumber = Mdata.mElecNumber;
		DisBetweenAdjacentEle = Mdata.mDisBetweenAdjacentEle;
		StepBetweenAdjacentEle = Mdata.mStepBetweenAdjacentEle;
		ElecPotential = Mdata.mElecPotential;
		RadiusAperture = Mdata.mRadiusAperture;
		GridAperture = Mdata.mGridAperture;
		RadiusBoundary = Mdata.mRadiusBoundary;
		GridBoundary = Mdata.mGridBoundary;
		IterationAccuracy = Mdata.mIterationAccuracy;
		NST = Mdata.mNST;
		EOE = Mdata.mEOE;
		dov = Mdata.mdov;

	}


public: struct inital_run
	{
		double* E;//电场分布
		int row, lin;//行，列
		int* PointGridElectrode;//电极位置对应格点坐标
		double* Z_StepGrid;//DisBetweenAdjacentEle轴方向网络步长划分
		double* Z_PointGrid;//DisBetweenAdjacentEle轴格点位置坐标
		double* Z_PointGridElec;//DisBetweenAdjacentEle轴等电位格点位置坐标
		double* R_StepGrid;//r轴方向网络步长划分
		double* R_PointGrid;//r轴格点位置坐标
		double* TransR_PointGrid;//转换后r轴格点位置坐标
		int* NumEquiPotential ;//等位点个数
		double* R_EquiPotentialPointGrid;//r轴等电位格点位置坐标
		double* TransR_EquiPotentialPointGrid;//转换后r轴等电位格点位置坐标
		double dov;

	};


public:	void run() {


	row = GridAperture + GridBoundary + 1;//Calculate the value of the row.
	lin = 0;//Assign a initial value to line so that it's final value can be calculate by a circulation
	for ( i = 0; i < ElecNumber; i++)
	{
		lin = lin +int( *(StepBetweenAdjacentEle + i));
	}
	lin = lin + ElecNumber;
	E = (double*)calloc(lin * row, sizeof(double));
	 s = 0;//列数上界	 
	 t = int(*StepBetweenAdjacentEle) +1 ;//列数下界

	for (int i = 0; i < row; i++)//行循环
	{
		s = 0;//列数上界
		t = int(*StepBetweenAdjacentEle) + 1;//列数下界
		if (i == row - 1)//第一行输入
		{
			for (int j = 0; j < ElecNumber; j++)//电极定位
			{
				if (j == 0)
				{
					int f = 1;//定义最终循环的次数


					for (k = s; k < t; k++)
					{

						//double d = *(ElecPotential + j) / *(StepBetweenAdjacentEle + j);//计算公差移到外面即可
						*(E + i * lin + k) = *(ElecPotential + j) / *(StepBetweenAdjacentEle + j) * (f - 1);
						f = f + 1;
					}
					*(E + i * lin + k - 1) = *(ElecPotential + j);
					*(E + i * lin + k) = *(ElecPotential + j);
					
					s = t + 1;
					t = t + int(*(StepBetweenAdjacentEle + j + 1)) + 1;
				}
				if (j > 0)
				{
					int f = 1;//定义最终循环的次数

					for (k = s; k < t; k++)
					{

						double d = (*(ElecPotential + j) - *(ElecPotential + j - 1)) / *(StepBetweenAdjacentEle + j);//计算公差移到外面即可
						*(E + i * lin + k) = *(ElecPotential + j - 1) + d* f;
						f = f + 1;
					}
					if (k < lin)
					{
						*(E + i * lin + k - 1) = *(ElecPotential + j);
						*(E + i * lin + k) = *(ElecPotential + j);
						
						s = t + 1;
						t = t + int(*(StepBetweenAdjacentEle + j + 1)) + 1;
					}
					if (k == lin)
					{
						*(E + i * lin + lin - 1) = *(ElecPotential + j);
					}

				}


			}
		}
		if (i < row - 1 && i >= GridAperture)
		{
			for (j = 0; j < ElecNumber; j++)//电极定位
			{


				for (k = s; k < t; k++)
				{


					*(E + i * lin + k) = 0;

				}
				//加判断语句规避光屏
				if (k < lin)
				{
					*(E + i * lin + k - 1) = *(ElecPotential + j);
					*(E + i * lin + k) = *(ElecPotential + j);
					
					s = t + 1;
					t = t + int (*(StepBetweenAdjacentEle + j + 1) )+1;
				}
				if (k == lin)
				{
					*(E + i * lin + lin - 1) = *(ElecPotential + j);
				}
			}
		}
		if (i < GridAperture)
		{
			for (j = 0; j < lin; j++)
			{
				*(E + i * lin + j) = 0;
			}
			*(E + i * lin + lin - 1) = *(ElecPotential + ElecNumber - 1);
		}
	}
	/*给电极位置对应格点坐标赋值*/
	PointGridElectrode = (int*)calloc(ElecNumber + 1, sizeof(int));//给电极位置对应格点坐标分配内存空间
	PointGridElectrode[0] = 0;
	for (i = 1; i < ElecNumber + 1; i++)
	{
		if (i == ElecNumber)
		{
			*(PointGridElectrode + i) = lin - 1;//荧光屏格点位置
		}
		else
		{
			*(PointGridElectrode + i) = *(PointGridElectrode + i - 1) + *(StepBetweenAdjacentEle + i - 1) + 1;//电极右格点位置
		}
	}
	/*给z轴方向网络步长划分赋值*/
	Z_StepGrid = (double*)calloc((lin - 1) * (row - 1), sizeof(double));//给z轴方向网络划分步长分配内存空间
	for (i = 0; i < row - 1; i++)
	{
		for (j = 1; j < ElecNumber + 1; j++)
		{
			if (j != ElecNumber)
			{
				for (k = *(PointGridElectrode + j - 1); k < *(PointGridElectrode + j); k++)
					*(Z_StepGrid + i * (lin - 1) + k) = *(DisBetweenAdjacentEle + j - 1) / *(StepBetweenAdjacentEle + j - 1);
				int a = *(PointGridElectrode + j);
				*(Z_StepGrid + i * (lin - 1) + a - 1) = ElecThickness;
			}
			if (j == ElecNumber)
			{
				for (k = *(PointGridElectrode + j - 1); k < *(PointGridElectrode + j); k++)
					*(Z_StepGrid + i * (lin - 1) + k) = *(DisBetweenAdjacentEle + j - 1) / *(StepBetweenAdjacentEle + j - 1);
			}
		}
	}
	/*给z轴格点位置坐标赋值*/
	Z_PointGrid = (double*)calloc(lin, sizeof(double));//给z轴格点位置坐标分配内存空间
	*Z_PointGrid = 0;
	for (i = 1; i < lin; i++)
		*(Z_PointGrid + i) = *(Z_PointGrid + i - 1) + *(Z_StepGrid + i - 1);
	/*给r轴方向网络步长划分赋值*/
	R_StepGrid = (double*)calloc((lin - 1) * (row - 1), sizeof(double));//给R_StepGrid分配内存空间
	for (i = 0; i < row - 1; i++)
	{
		if (i < GridAperture)
		{
			for (j = 0; j < lin - 1; j++)
				*(R_StepGrid + i * (lin - 1) + j) = RadiusAperture / GridAperture;
		}
		else
		{
			for (j = 0; j < lin - 1; j++)
				*(R_StepGrid + i * (lin - 1) + j) = RadiusBoundary / GridBoundary;
		}

	}
	/*给r轴格点位置赋值*/
	R_PointGrid = (double*)calloc(row, sizeof(double));//
	*R_PointGrid = 0;
	for (i = 1; i < row; i++)
		*(R_PointGrid + i) = *(R_PointGrid + i - 1) + *(R_StepGrid + (i - 1) * (lin - 1));

}


public:  inital_run getValue() {
	inital_run Mdata;
	Mdata.E = E;
	Mdata.row = row;
	Mdata.lin = lin;
	Mdata.PointGridElectrode = PointGridElectrode;
	Mdata.Z_StepGrid = Z_StepGrid;
	Mdata.Z_PointGridElec = Z_PointGridElec;
	Mdata.R_StepGrid = R_StepGrid;
	Mdata.R_PointGrid = R_PointGrid;
	Mdata.TransR_PointGrid = TransR_PointGrid;
	Mdata.NumEquiPotential  = NumEquiPotential ;
	Mdata.R_EquiPotentialPointGrid = R_EquiPotentialPointGrid;
	Mdata.TransR_EquiPotentialPointGrid = TransR_EquiPotentialPointGrid;
	Mdata.Z_PointGrid = Z_PointGrid;
	Mdata.dov = dov;



	return Mdata;

}

public: void show(inital_run data) {
	std::cout << "行"<<data.row << std::endl;
	std::cout <<"列"<< data.lin << std::endl;
	std::cout << "电极位置对应格点坐标"<<data.PointGridElectrode << std::endl;
	std::cout <<"z轴方向网络步长划分"<< data.Z_StepGrid << std::endl;
	std::cout << "z轴格点位置坐标"<<data.Z_PointGrid << std::endl;
	std::cout << "z轴等电位格点位置坐标"<<data.Z_PointGridElec << std::endl;
	std::cout << "r轴方向网络步长划分"<<data.R_StepGrid << std::endl;
	std::cout << "r轴格点位置坐标"<<data.R_PointGrid << std::endl;
	std::cout << "转换后r轴格点位置坐标"<<data.TransR_PointGrid << std::endl;
	std::cout << "等位点个数"<<data.NumEquiPotential  << std::endl;
	std::cout <<"r轴等电位格点位置坐标"<<data.R_EquiPotentialPointGrid << std::endl;

	std::cout <<"转换后r轴等电位格点位置坐标"<< data.TransR_EquiPotentialPointGrid << std::endl;
	std::cout << "电场分布"<< data.E << std::endl;



}


public:	  void show_data() {

		  std::cout << "电极厚度：" << ElecThickness << std::endl;
		  std::cout << "电极总数：" << ElecNumber << std::endl;
		  std::cout << "相邻电极之间的距离：" << std::endl;
		  for (int i = 0; i < ElecNumber; i++) {
			  std::cout << DisBetweenAdjacentEle[i] << std::endl;
		  }


		  std::cout << "相邻电极划分的步长：" << std::endl;
		  for (int i = 0; i < ElecNumber - 1; i++) {
			  std::cout <<StepBetweenAdjacentEle[i] << std::endl;
		  }

		  std::cout << *StepBetweenAdjacentEle << std::endl;

		  std::cout << "电极电位：" << std::endl;
		  for (int i = 0; i < ElecNumber; i++) {
			  std::cout << ElecPotential[i] << std::endl;
		  }


		  std::cout << "电极内孔径半径：" << RadiusAperture << std::endl;
		  std::cout << "电极内孔径半径等步长划分的网格数：" << GridAperture << std::endl;
		  std::cout << "从电极内孔沿到封闭边界处的径向距离：" << RadiusBoundary << std::endl;
		  std::cout << "从电极内孔沿到封闭边界处划分的网格数:" << GridBoundary << std::endl;
		  std::cout << "迭代精度:" << IterationAccuracy << std::endl;
		  std::cout << "输出打印空间电位时网格点间隔数:" << NST << std::endl;
		  std::cout << "要求扫描等位线间隔值：" << dov << std::endl;

		/*  std::cout << "要求扫描搜索等电位线的电位间隔值:" << std::endl;
		  for (int i = 0; i < ElecNumber; i++) {
			  std::cout << EOE[i] << std::endl;
		  }*/







	  }
};


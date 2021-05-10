#pragma once
#include "Read.h"
#include "inital.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <graphics.h>
#include <conio.h>
#include <iostream>
#include <cstdio>

class loop
{

public:	struct newData {
	public:
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

		double RadiusAperture;//电极内孔径半径
		double RadiusBoundary;//从电极内孔沿到封闭边界处的径向距离
		int m;
	};


	double c1, c2, c3, c4, c0;    //系数
	double E2;                //迭代后的格点电位
	int turn;                    //迭代轮次
	int times;                 //迭代次数
	int atime;              //迭代总次数

	int flag;

	double w;                 //迭代因子
	double* dV;              //格点电位残差
	double sdV;               //格点电位残差和
	double mdV;               //格点电位最大残差
	double adV;               //格点电位平均残差
	double adV11;              //第二轮12次迭代格点电位平均残差
	double adV12;              //第二轮第13次迭代格点电位平均残差
	double Lambda;             //残差均值比值
	double uLambda;           //中间变量
	double wLambda;
	double wm;                //修正后的w
	double wm1;               //修正后的w
	int    judge;             //判断子


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

	int i, j, k;
	double t, s;
	double a = 0;
	int m;



	FILE* OutputFile;

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

public:newData  loop_run() {
		/*开始计算像管电位分布*/
			/*第一轮迭代*/
	fopen_s(&OutputFile, "E:\\Learning(School)\\光电成像原理与技术\\光电成像实验\\程序\\11.res", "w+");
		turn = 1;
		times = 1;
		atime = 0;
		w = 1;
		SOR();
		atime++;

		/*第二轮迭代*/
		w = 1.375;
		turn = turn + 1;
		for (times = 1; times < 13; times++)
		{
			SOR();
			atime = atime + 1;
			if (times == 11)
				adV11 = adV;
			if (times == 12)
				adV12 = adV;
		}
		Lambda = adV12 / adV11;
		uLambda = (Lambda + w - 1) / (sqrt(Lambda) * w);
		wLambda = 2 / (1 + sqrt(1 - uLambda * uLambda));
		wm = 1.25 * wLambda - 0.5;
		/*迭代若干轮以固定迭代因子*/
		do
		{
			w = wm;
			turn++;
			for (times = 1; times < 13; times++)
			{
				SOR();
				atime = atime + 1;
				if (times == 11)
					adV11 = adV;
				if (times == 12)
					adV12 = adV;
			}
			wm1 = wm;
			Lambda = adV12 / adV11;
			uLambda = (Lambda + w - 1) / (sqrt(Lambda) * w);
			wLambda = 2 / (1 + sqrt(1 - uLambda * uLambda));
			wm = 1.25 * wLambda - 0.5;
		} while (fabs((wm - wm1) / (2 - wm1)) >= 0.05);
		w = wm;
		/*继续多轮迭代以达到精度要求*/
		do
		{
			turn++;
			SOR();   //最佳迭代因子确定后继续迭代
			atime++;
			judge = 0;
			for (i = 0; i < row; i++)
				for (j = 0; j < lin; j++)
				{
					if (*(dV + i * lin + j) >= IterationAccuracy)
						judge++;
				}
		} while (judge > 0);
	
		//for (i = row - 1; i >= 0; i--)
		//{
		//	for (j = 0; j < lin - 1; j++)
		//	{
		//		fprintf(OutputFile, "%3.0lf", *(E + i * lin + j));
		//	}
		//	fprintf(OutputFile, "%3.0lf\ElecNumber", *(E + i * lin + lin - 1));
		//}

		////////////////////////////////////////////
		   /*扫描等位点坐标*/

		if (dov != 0)
		{
			for (m = 1; a < (ElecPotential[ElecNumber - 1] - dov); m++)   //计算电位数量
				a = a + dov;
			/*EOE = (double*)calloc(m - 1, sizeof(double));*/

			*EOE = dov;
			for (i = 1; i < m - 1; i++)
				EOE[i] = EOE[i - 1] + dov;
		}
		else if (dov == 0)
		{
			m = 0;
			for (i = 0; i < 100; i++)
			{
				if (EOE[i] != 0)
					m++;

			}
			m = m + 1;
			
		}
		flag = m;
		//for (m = 1; a < (ElecPotential[ElecNumber - 1] - dov); m++)   //计算电位数量
		//	a = a + dov;
		//EOE = (double*)calloc(m - 1, sizeof(double));
		R_EquiPotentialPointGrid = (double*)calloc((m - 1) * 2 * row, sizeof(double));
		TransR_EquiPotentialPointGrid = (double*)calloc((m - 1) * 2 * row, sizeof(double));
		Z_PointGridElec = (double*)calloc(((m - 1)) * 2 * lin, sizeof(double));
		/**EOE = dov;
		for (i = 1; i < m - 1; i++)
			EOE[i] = EOE[i - 1] + dov;*/
		NumEquiPotential  = (int*)calloc(m - 1, sizeof(int));
		for (i = 0; i < m - 1; i++)
			NumEquiPotential [i] = 0;
		
		for (j = 1; j < lin - 1; j++)
		{
			for (k = 0; k < m - 1; k++)
			{
				for (i = row - 1; i > 0; i--)
				{
					if (((*(E + i * lin + j) < EOE[k]) && (*(E + i * lin + j - 1) > EOE[k])) || ((*(E + i * lin + j) > EOE[k]) && (*(E + i * lin + j - 1) < EOE[k])))
					{
						R_EquiPotentialPointGrid[k * 2 * row + NumEquiPotential [k]] = R_PointGrid[i];
						Z_PointGridElec[k * 2 * lin + NumEquiPotential [k]] = Z_PointGrid[j] + ((Z_PointGrid[j - 1] - Z_PointGrid[j]) / (*(E + i * lin + j - 1) - *(E + i * lin + j))) * (EOE[k] - *(E + i * lin + j));
						NumEquiPotential [k]++;
					}
					if (*(E + i * lin + j) == EOE[k])
					{
						R_EquiPotentialPointGrid[k * 2 * row + NumEquiPotential [k]] = R_PointGrid[i];
						Z_PointGridElec[k * 2 * lin + NumEquiPotential [k]] = Z_PointGrid[j];
						NumEquiPotential [k]++;
					}
					if (((*(E + i * lin + j) < EOE[k]) && (*(E + (i - 1) * lin + j) > EOE[k])) || ((*(E + i * lin + j) > EOE[k]) && (*(E + (i - 1) * lin + j) < EOE[k])))
					{
						Z_PointGridElec[k * 2 * lin + NumEquiPotential [k]] = Z_PointGrid[j];
						R_EquiPotentialPointGrid[k * 2 * row + NumEquiPotential [k]] = R_PointGrid[i] + ((R_PointGrid[i - 1] - R_PointGrid[i]) / (*(E + (i - 1) * lin + j) - *(E + i * lin + j))) * (EOE[k] - *(E + i * lin + j));
						NumEquiPotential [k]++;
					}
				}
				i = 0;
				if (((*(E + i * lin + j) < EOE[k]) && (*(E + i * lin + j - 1) > EOE[k])) || ((*(E + i * lin + j) > EOE[k]) && (*(E + i * lin + j - 1) < EOE[k])))
				{
					R_EquiPotentialPointGrid[k * 2 * row + NumEquiPotential [k]] = R_PointGrid[i];
					Z_PointGridElec[k * 2 * lin + NumEquiPotential [k]] = Z_PointGrid[j] + ((Z_PointGrid[j - 1] - Z_PointGrid[j]) / (*(E + i * lin + j - 1) - *(E + i * lin + j))) * (EOE[k] - *(E + i * lin + j));
					NumEquiPotential [k]++;
				}
				if (*(E + i * lin + j) == EOE[k])
				{
					R_EquiPotentialPointGrid[k * 2 * row + NumEquiPotential [k]] = R_PointGrid[i];
					Z_PointGridElec[k * 2 * lin + NumEquiPotential [k]] = Z_PointGrid[j];
					NumEquiPotential[k]++;
				}
			}
		}

		//
		//for (k = 0; k < m - 1; k++)
		//{
		//	fprintf(OutputFile, "第%d电位,电位值：%lf,\ElecNumber", k + 1, EOE[k]);
		//	for (j = 0; j < NumEquiPotential [k] - 1; j++)
		//	{
		//		fprintf(OutputFile, "(%3.2lf,%3.2lf)\t", Z_PointGridElec[k * 2 * lin + j], R_EquiPotentialPointGrid[k * 2 * row + j]);
		//	}
		//	fprintf(OutputFile, "(%3.2lf,%3.2lf)\t\ElecNumber", Z_PointGridElec[k * 2 * lin + NumEquiPotential [k] - 1], R_EquiPotentialPointGrid[k * 2 * row + NumEquiPotential [k] - 1]);
		//}//输出等位线坐标
		//fclose;

		newData DATA;
		DATA.lin = lin;
		DATA.row = row;
		DATA.NumEquiPotential  = NumEquiPotential ;
		DATA.R_StepGrid = R_StepGrid;
		DATA.Z_StepGrid = Z_StepGrid;
		DATA.PointGridElectrode = PointGridElectrode; 
		DATA.R_PointGrid = R_PointGrid;
		DATA.R_EquiPotentialPointGrid = R_EquiPotentialPointGrid;
		DATA.Z_PointGrid = Z_PointGrid;
		DATA.Z_PointGridElec = Z_PointGridElec;
		DATA.TransR_PointGrid = TransR_PointGrid;
		DATA.TransR_EquiPotentialPointGrid = TransR_EquiPotentialPointGrid;
		
		DATA.RadiusAperture = RadiusAperture;

		DATA.RadiusBoundary = RadiusBoundary;
		DATA.m = flag;

		
		return DATA;
		setfree();
	}

	void SOR( ){
		

		dV = (double*)calloc(row * lin, sizeof(double));
		for (int i = 0; i < row - 1; i++)
		{
			if (i == 0)
			{

				for (int j = 1; j < lin - 1; j++)
				{
					c1 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j - 1)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));   //系数计算
					c2 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));
					c3 = 0;
					c4 = 4 / ((*(R_StepGrid + i * (lin - 1) + j)) * (*(R_StepGrid + i * (lin - 1) + j)));
					c0 = c1 + c2 + c3 + c4;
					E2 = (1 - w) * (*(E + i * lin + j)) + w * (c1 * (*(E + i * lin + j - 1)) + c2 * (*(E + i * lin + j + 1)) + c3 * (*(E + (i - 1) * lin + j)) + c4 * (*(E + (i + 1) * lin + j))) / c0;
					*(dV + i * lin + j) = fabs(E2 - *(E + i * lin + j));
					*(E + i * lin + j) = E2;

				}

			}
			if (i < GridAperture && i>0)
			{
				for (int j = 1; j < lin - 1; j++)
				{
					c1 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j - 1)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));   //系数计算
					c2 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));
					c3 = (2 * (*(R_PointGrid + i)) - *(R_StepGrid + i * (lin - 1) + j)) / ((*(R_PointGrid + i)) * (*(R_StepGrid + (i - 1) * (lin - 1) + j)) * (*(R_StepGrid + (i - 1) * (lin - 1) + j) + (*(R_StepGrid + i * (lin - 1) + j))));
					c4 = (2 * (*(R_PointGrid + i)) + (*(R_StepGrid + (i - 1) * (lin - 1) + j))) / ((*(R_PointGrid + i)) * (*(R_StepGrid + i * (lin - 1) + j)) * (*(R_StepGrid + (i - 1) * (lin - 1) + j) + *(R_StepGrid + i * (lin - 1) + j)));
					c0 = c1 + c2 + c3 + c4;
					E2 = (1 - w) * (*(E + i * lin + j)) + w * (c1 * (*(E + i * lin + j - 1)) + c2 * (*(E + i * lin + j + 1)) + c3 * (*(E + (i - 1) * lin + j)) + c4 * (*(E + (i + 1) * lin + j))) / c0;
					*(dV + i * lin + j) = fabs(E2 - *(E + i * lin + j));
					*(E + i * lin + j) = E2;
				}
			}
			if (i >= GridAperture)
			{
				for (int k = 1; k < ElecNumber + 1; k++)//电极循环
				{
					if (k != ElecNumber)
					{
						for (int j = *(PointGridElectrode + k - 1) + 1; j < *(PointGridElectrode + k) - 1; j++)//列数循环，规避边界点
						{
							c1 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j - 1)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));   //系数计算
							c2 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));
							c3 = (2 * (*(R_PointGrid + i)) - *(R_StepGrid + i * (lin - 1) + j)) / ((*(R_PointGrid + i)) * (*(R_StepGrid + (i - 1) * (lin - 1) + j)) * (*(R_StepGrid + (i - 1) * (lin - 1) + j) + (*(R_StepGrid + i * (lin - 1) + j))));
							c4 = (2 * (*(R_PointGrid + i)) + (*(R_StepGrid + (i - 1) * (lin - 1) + j))) / ((*(R_PointGrid + i)) * (*(R_StepGrid + i * (lin - 1) + j)) * (*(R_StepGrid + (i - 1) * (lin - 1) + j) + *(R_StepGrid + i * (lin - 1) + j)));
							c0 = c1 + c2 + c3 + c4;
							E2 = (1 - w) * (*(E + i * lin + j)) + w * (c1 * (*(E + i * lin + j - 1)) + c2 * (*(E + i * lin + j + 1)) + c3 * (*(E + (i - 1) * lin + j)) + c4 * (*(E + (i + 1) * lin + j))) / c0;
							*(dV + i * lin + j) = fabs(E2 - *(E + i * lin + j));
							*(E + i * lin + j) = E2;
						}
					}
					else
					{
						for (int j = *(PointGridElectrode + k - 1) + 1; j < *(PointGridElectrode + k); j++)
						{
							c1 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j - 1)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));   //系数计算
							c2 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));
							c3 = (2 * (*(R_PointGrid + i)) - *(R_StepGrid + i * (lin - 1) + j)) / ((*(R_PointGrid + i)) * (*(R_StepGrid + (i - 1) * (lin - 1) + j)) * (*(R_StepGrid + (i - 1) * (lin - 1) + j) + (*(R_StepGrid + i * (lin - 1) + j))));
							c4 = (2 * (*(R_PointGrid + i)) + (*(R_StepGrid + (i - 1) * (lin - 1) + j))) / ((*(R_PointGrid + i)) * (*(R_StepGrid + i * (lin - 1) + j)) * (*(R_StepGrid + (i - 1) * (lin - 1) + j) + *(R_StepGrid + i * (lin - 1) + j)));
							c0 = c1 + c2 + c3 + c4;
							E2 = (1 - w) * (*(E + i * lin + j)) + w * (c1 * (*(E + i * lin + j - 1)) + c2 * (*(E + i * lin + j + 1)) + c3 * (*(E + (i - 1) * lin + j)) + c4 * (*(E + (i + 1) * lin + j))) / c0;
							*(dV + i * lin + j) = fabs(E2 - *(E + i * lin + j));
							*(E + i * lin + j) = E2;
						}
					}
				}

			}
		}
		sdV = 0;
		mdV = 0;
		for (int i = 0; i < row * lin; i++)
		{
			sdV = sdV + *(dV + i);
			if (*(dV + i) > mdV)
				mdV = *(dV + i);
		}
		adV = sdV / (row * lin);
		//free(dV);

}

	void setfree() {
		free(DisBetweenAdjacentEle);
		free(StepBetweenAdjacentEle);
		free(ElecPotential);
		free(EOE);
		free(E);
		free(PointGridElectrode);
		free(Z_StepGrid);
		free(Z_PointGrid);
		free(Z_PointGridElec);
		free(R_StepGrid);
		free(R_PointGrid);
		free(TransR_PointGrid);
		free(NumEquiPotential );
		free(R_EquiPotentialPointGrid);
		free(TransR_EquiPotentialPointGrid);

	}

public:	void setValue(inital::inital_run Mdata ,Read::data Mdata1) {
		E = Mdata.E;//电场分布
		row = Mdata.row;
		lin = Mdata.lin;//行，列
		PointGridElectrode = Mdata.PointGridElectrode;//电极位置对应格点坐标
		Z_StepGrid = Mdata.Z_StepGrid;//DisBetweenAdjacentEle轴方向网络步长划分
		Z_PointGrid = Mdata.Z_PointGrid;//DisBetweenAdjacentEle轴格点位置坐标
		Z_PointGridElec = Mdata.Z_PointGridElec;//DisBetweenAdjacentEle轴等电位格点位置坐标
		R_StepGrid = Mdata.R_StepGrid;//r轴方向网络步长划分
		R_PointGrid = Mdata.R_PointGrid;//r轴格点位置坐标
		TransR_PointGrid = Mdata.TransR_PointGrid;//转换后r轴格点位置坐标
		NumEquiPotential  = Mdata.NumEquiPotential ;//等位点个数
		R_EquiPotentialPointGrid = Mdata.R_EquiPotentialPointGrid;//r轴等电位格点位置坐标
		TransR_EquiPotentialPointGrid = Mdata.TransR_EquiPotentialPointGrid;//转换后r轴等电位格点位置坐标


		ElecThickness = Mdata1.mElecThickness;//电极厚度
		ElecNumber = Mdata1.mElecNumber;//电极总数
		DisBetweenAdjacentEle = Mdata1.mDisBetweenAdjacentEle;//相邻电极之间的距离
		StepBetweenAdjacentEle = Mdata1.mStepBetweenAdjacentEle;//相邻电极划分的步长
		ElecPotential = Mdata1.mElecPotential;//电极电位
		RadiusAperture = Mdata1.mRadiusAperture;//电极内孔径半径
		GridAperture = Mdata1.mGridAperture;//电极内孔径半径等步长划分的网格数
		RadiusBoundary = Mdata1.mRadiusBoundary;//从电极内孔沿到封闭边界处的径向距离
		GridBoundary = Mdata1.mGridBoundary;//从电极内孔沿到封闭边界处划分的网格数
		IterationAccuracy = Mdata1.mIterationAccuracy;//迭代精度
		NST = Mdata1.mNST;//输出打印空间电位时网格点间隔数
		EOE = Mdata1.mEOE;//要求扫描搜索等电位线的电位间隔值
		dov = Mdata1.mdov;//要求扫描等位线间隔值


	}

public:	void draw(newData hh) {

		////////为方便制图，开始坐标变换////////////
   /*由于该程序将光轴与阴极面交点设为原点，而画图软件则是将左上角定位原点，因此为了观察习惯，将r轴方向坐标进行变换*/
   /*TransR_PointGrid = (double*)calloc(row, sizeof(double));
   for (i = 0; i < row; i++)
	   TransR_PointGrid[i] = RadiusAperture + RadiusBoundary - R_PointGrid[i];*/


		for (k = 0; k < hh.m - 1; k++)
		{
			for (j = 0; j < hh.NumEquiPotential [k]; j++)
			{
				TransR_EquiPotentialPointGrid[k * 2 * hh.row + j] = hh.RadiusAperture + hh.RadiusBoundary - hh.R_EquiPotentialPointGrid[k * 2 * row + j];

			}
		}
		/*for (j = 0; j < NumEquiPotential [0]; j++)
			printf("%lf\t", R_EquiPotentialPointGrid[j]);
		printf("\ElecNumber");
		for (j = 0; j < NumEquiPotential [0]; j++)
		printf("%lf\t",	TransR_EquiPotentialPointGrid[ j] );*/
		////////////////////////////////////////////

			 /*开始绘制等位线*/
		initgraph(110 + int(10 * hh.Z_PointGrid[hh.lin - 1]), 110 + int(10 * hh.R_PointGrid[hh.row - 1]));
		setbkcolor(BLACK);
		cleardevice();
		rectangle(50, 50, 50 + int(10 * hh.Z_PointGrid[hh.lin - 1]), 50 + int(10 * hh.R_PointGrid[hh.row - 1]));
		for (i = 1; i < ElecNumber; i++)
			bar(50 + (10 * hh.Z_PointGrid[hh.PointGridElectrode[i] - 1]), 50, 50 + (10 * (hh.Z_PointGrid[hh.PointGridElectrode[i]])), 50 + (10 * hh.RadiusBoundary));

		for (k = 0; k < hh.m - 1; k++)
		{
			for (j = 0; j < NumEquiPotential [k] - 1; j++)
				line(50 + 10 * Z_PointGridElec[k * 2 * hh.lin + j], 50 + 10 * (TransR_EquiPotentialPointGrid[k * 2 * row + j]), 50 + 10 * Z_PointGridElec[k * 2 * lin + j + 1], 50 + 10 * (TransR_EquiPotentialPointGrid[k * 2 * row + j + 1]));
			int R = rand() % 255, G = rand() % 255, B = rand() % 255;
			setlinecolor(RGB(R, G, B));
		}

		cout << "test" << endl;
		//_getch();
		closegraph();
		//_getch();
		system("pause");
		
	}


};


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
		double* E;//�糡�ֲ�
		int row, lin;//�У���
		int* PointGridElectrode;//�缫λ�ö�Ӧ�������
		double* Z_StepGrid;//DisBetweenAdjacentEle�᷽�����粽������
		double* Z_PointGrid;//DisBetweenAdjacentEle����λ������
		double* Z_PointGridElec;//DisBetweenAdjacentEle��ȵ�λ���λ������
		double* R_StepGrid;//r�᷽�����粽������
		double* R_PointGrid;//r����λ������
		double* TransR_PointGrid;//ת����r����λ������
		int* NumEquiPotential ;//��λ�����
		double* R_EquiPotentialPointGrid;//r��ȵ�λ���λ������
		double* TransR_EquiPotentialPointGrid;//ת����r��ȵ�λ���λ������

		double RadiusAperture;//�缫�ڿ׾��뾶
		double RadiusBoundary;//�ӵ缫�ڿ��ص���ձ߽紦�ľ������
		int m;
	};


	double c1, c2, c3, c4, c0;    //ϵ��
	double E2;                //������ĸ���λ
	int turn;                    //�����ִ�
	int times;                 //��������
	int atime;              //�����ܴ���

	int flag;

	double w;                 //��������
	double* dV;              //����λ�в�
	double sdV;               //����λ�в��
	double mdV;               //����λ���в�
	double adV;               //����λƽ���в�
	double adV11;              //�ڶ���12�ε�������λƽ���в�
	double adV12;              //�ڶ��ֵ�13�ε�������λƽ���в�
	double Lambda;             //�в��ֵ��ֵ
	double uLambda;           //�м����
	double wLambda;
	double wm;                //�������w
	double wm1;               //�������w
	int    judge;             //�ж���


	double* E;//�糡�ֲ�
	int row, lin;//�У���
	int* PointGridElectrode;//�缫λ�ö�Ӧ�������
	double* Z_StepGrid;//DisBetweenAdjacentEle�᷽�����粽������
	double* Z_PointGrid;//DisBetweenAdjacentEle����λ������
	double* Z_PointGridElec;//DisBetweenAdjacentEle��ȵ�λ���λ������
	double* R_StepGrid;//r�᷽�����粽������
	double* R_PointGrid;//r����λ������
	double* TransR_PointGrid;//ת����r����λ������
	int* NumEquiPotential ;//��λ�����
	double* R_EquiPotentialPointGrid;//r��ȵ�λ���λ������
	double* TransR_EquiPotentialPointGrid;//ת����r��ȵ�λ���λ������

	int i, j, k;
	double t, s;
	double a = 0;
	int m;



	FILE* OutputFile;

	double ElecThickness;//�缫���
	int ElecNumber;//�缫����
	double* DisBetweenAdjacentEle;//���ڵ缫֮��ľ���
	double* StepBetweenAdjacentEle;//���ڵ缫���ֵĲ���
	double* ElecPotential;//�缫��λ
	double RadiusAperture;//�缫�ڿ׾��뾶
	int GridAperture;//�缫�ڿ׾��뾶�Ȳ������ֵ�������
	double RadiusBoundary;//�ӵ缫�ڿ��ص���ձ߽紦�ľ������
	int GridBoundary;//�ӵ缫�ڿ��ص���ձ߽紦���ֵ�������
	double IterationAccuracy;//��������
	int NST;//�����ӡ�ռ��λʱ���������
	double* EOE;//Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ
	double dov;//Ҫ��ɨ���λ�߼��ֵ

public:newData  loop_run() {
		/*��ʼ������ܵ�λ�ֲ�*/
			/*��һ�ֵ���*/
	fopen_s(&OutputFile, "E:\\Learning(School)\\������ԭ���뼼��\\������ʵ��\\����\\11.res", "w+");
		turn = 1;
		times = 1;
		atime = 0;
		w = 1;
		SOR();
		atime++;

		/*�ڶ��ֵ���*/
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
		/*�����������Թ̶���������*/
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
		/*�������ֵ����Դﵽ����Ҫ��*/
		do
		{
			turn++;
			SOR();   //��ѵ�������ȷ�����������
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
		   /*ɨ���λ������*/

		if (dov != 0)
		{
			for (m = 1; a < (ElecPotential[ElecNumber - 1] - dov); m++)   //�����λ����
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
		//for (m = 1; a < (ElecPotential[ElecNumber - 1] - dov); m++)   //�����λ����
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
		//	fprintf(OutputFile, "��%d��λ,��λֵ��%lf,\ElecNumber", k + 1, EOE[k]);
		//	for (j = 0; j < NumEquiPotential [k] - 1; j++)
		//	{
		//		fprintf(OutputFile, "(%3.2lf,%3.2lf)\t", Z_PointGridElec[k * 2 * lin + j], R_EquiPotentialPointGrid[k * 2 * row + j]);
		//	}
		//	fprintf(OutputFile, "(%3.2lf,%3.2lf)\t\ElecNumber", Z_PointGridElec[k * 2 * lin + NumEquiPotential [k] - 1], R_EquiPotentialPointGrid[k * 2 * row + NumEquiPotential [k] - 1]);
		//}//�����λ������
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
					c1 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j - 1)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));   //ϵ������
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
					c1 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j - 1)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));   //ϵ������
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
				for (int k = 1; k < ElecNumber + 1; k++)//�缫ѭ��
				{
					if (k != ElecNumber)
					{
						for (int j = *(PointGridElectrode + k - 1) + 1; j < *(PointGridElectrode + k) - 1; j++)//����ѭ������ܱ߽��
						{
							c1 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j - 1)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));   //ϵ������
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
							c1 = 2 / ((*(Z_StepGrid + i * (lin - 1) + j - 1)) * (*(Z_StepGrid + i * (lin - 1) + j - 1) + *(Z_StepGrid + i * (lin - 1) + j)));   //ϵ������
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
		E = Mdata.E;//�糡�ֲ�
		row = Mdata.row;
		lin = Mdata.lin;//�У���
		PointGridElectrode = Mdata.PointGridElectrode;//�缫λ�ö�Ӧ�������
		Z_StepGrid = Mdata.Z_StepGrid;//DisBetweenAdjacentEle�᷽�����粽������
		Z_PointGrid = Mdata.Z_PointGrid;//DisBetweenAdjacentEle����λ������
		Z_PointGridElec = Mdata.Z_PointGridElec;//DisBetweenAdjacentEle��ȵ�λ���λ������
		R_StepGrid = Mdata.R_StepGrid;//r�᷽�����粽������
		R_PointGrid = Mdata.R_PointGrid;//r����λ������
		TransR_PointGrid = Mdata.TransR_PointGrid;//ת����r����λ������
		NumEquiPotential  = Mdata.NumEquiPotential ;//��λ�����
		R_EquiPotentialPointGrid = Mdata.R_EquiPotentialPointGrid;//r��ȵ�λ���λ������
		TransR_EquiPotentialPointGrid = Mdata.TransR_EquiPotentialPointGrid;//ת����r��ȵ�λ���λ������


		ElecThickness = Mdata1.mElecThickness;//�缫���
		ElecNumber = Mdata1.mElecNumber;//�缫����
		DisBetweenAdjacentEle = Mdata1.mDisBetweenAdjacentEle;//���ڵ缫֮��ľ���
		StepBetweenAdjacentEle = Mdata1.mStepBetweenAdjacentEle;//���ڵ缫���ֵĲ���
		ElecPotential = Mdata1.mElecPotential;//�缫��λ
		RadiusAperture = Mdata1.mRadiusAperture;//�缫�ڿ׾��뾶
		GridAperture = Mdata1.mGridAperture;//�缫�ڿ׾��뾶�Ȳ������ֵ�������
		RadiusBoundary = Mdata1.mRadiusBoundary;//�ӵ缫�ڿ��ص���ձ߽紦�ľ������
		GridBoundary = Mdata1.mGridBoundary;//�ӵ缫�ڿ��ص���ձ߽紦���ֵ�������
		IterationAccuracy = Mdata1.mIterationAccuracy;//��������
		NST = Mdata1.mNST;//�����ӡ�ռ��λʱ���������
		EOE = Mdata1.mEOE;//Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ
		dov = Mdata1.mdov;//Ҫ��ɨ���λ�߼��ֵ


	}

public:	void draw(newData hh) {

		////////Ϊ������ͼ����ʼ����任////////////
   /*���ڸó��򽫹����������潻����Ϊԭ�㣬����ͼ������ǽ����ϽǶ�λԭ�㣬���Ϊ�˹۲�ϰ�ߣ���r�᷽��������б任*/
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

			 /*��ʼ���Ƶ�λ��*/
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


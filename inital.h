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



	double* E;//�糡�ֲ�
	int row;
	int lin;//�У���
	int* PointGridElectrode;//�缫λ�ö�Ӧ�������
	double* Z_StepGrid;//z�᷽�����粽������
	double* Z_PointGrid;//z����λ������
	double* Z_PointGridElec;//z��ȵ�λ���λ������
	double* R_StepGrid;//r�᷽�����粽������
	double* R_PointGrid;//r����λ������
	double* TransR_PointGrid;//ת����r����λ������
	int* NumEquiPotential ;//��λ�����
	double* R_EquiPotentialPointGrid;//r��ȵ�λ���λ������
	double* TransR_EquiPotentialPointGrid;//ת����r��ȵ�λ���λ������


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
	 s = 0;//�����Ͻ�	 
	 t = int(*StepBetweenAdjacentEle) +1 ;//�����½�

	for (int i = 0; i < row; i++)//��ѭ��
	{
		s = 0;//�����Ͻ�
		t = int(*StepBetweenAdjacentEle) + 1;//�����½�
		if (i == row - 1)//��һ������
		{
			for (int j = 0; j < ElecNumber; j++)//�缫��λ
			{
				if (j == 0)
				{
					int f = 1;//��������ѭ���Ĵ���


					for (k = s; k < t; k++)
					{

						//double d = *(ElecPotential + j) / *(StepBetweenAdjacentEle + j);//���㹫���Ƶ����漴��
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
					int f = 1;//��������ѭ���Ĵ���

					for (k = s; k < t; k++)
					{

						double d = (*(ElecPotential + j) - *(ElecPotential + j - 1)) / *(StepBetweenAdjacentEle + j);//���㹫���Ƶ����漴��
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
			for (j = 0; j < ElecNumber; j++)//�缫��λ
			{


				for (k = s; k < t; k++)
				{


					*(E + i * lin + k) = 0;

				}
				//���ж�����ܹ���
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
	/*���缫λ�ö�Ӧ������긳ֵ*/
	PointGridElectrode = (int*)calloc(ElecNumber + 1, sizeof(int));//���缫λ�ö�Ӧ�����������ڴ�ռ�
	PointGridElectrode[0] = 0;
	for (i = 1; i < ElecNumber + 1; i++)
	{
		if (i == ElecNumber)
		{
			*(PointGridElectrode + i) = lin - 1;//ӫ�������λ��
		}
		else
		{
			*(PointGridElectrode + i) = *(PointGridElectrode + i - 1) + *(StepBetweenAdjacentEle + i - 1) + 1;//�缫�Ҹ��λ��
		}
	}
	/*��z�᷽�����粽�����ָ�ֵ*/
	Z_StepGrid = (double*)calloc((lin - 1) * (row - 1), sizeof(double));//��z�᷽�����绮�ֲ��������ڴ�ռ�
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
	/*��z����λ�����긳ֵ*/
	Z_PointGrid = (double*)calloc(lin, sizeof(double));//��z����λ����������ڴ�ռ�
	*Z_PointGrid = 0;
	for (i = 1; i < lin; i++)
		*(Z_PointGrid + i) = *(Z_PointGrid + i - 1) + *(Z_StepGrid + i - 1);
	/*��r�᷽�����粽�����ָ�ֵ*/
	R_StepGrid = (double*)calloc((lin - 1) * (row - 1), sizeof(double));//��R_StepGrid�����ڴ�ռ�
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
	/*��r����λ�ø�ֵ*/
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
	std::cout << "��"<<data.row << std::endl;
	std::cout <<"��"<< data.lin << std::endl;
	std::cout << "�缫λ�ö�Ӧ�������"<<data.PointGridElectrode << std::endl;
	std::cout <<"z�᷽�����粽������"<< data.Z_StepGrid << std::endl;
	std::cout << "z����λ������"<<data.Z_PointGrid << std::endl;
	std::cout << "z��ȵ�λ���λ������"<<data.Z_PointGridElec << std::endl;
	std::cout << "r�᷽�����粽������"<<data.R_StepGrid << std::endl;
	std::cout << "r����λ������"<<data.R_PointGrid << std::endl;
	std::cout << "ת����r����λ������"<<data.TransR_PointGrid << std::endl;
	std::cout << "��λ�����"<<data.NumEquiPotential  << std::endl;
	std::cout <<"r��ȵ�λ���λ������"<<data.R_EquiPotentialPointGrid << std::endl;

	std::cout <<"ת����r��ȵ�λ���λ������"<< data.TransR_EquiPotentialPointGrid << std::endl;
	std::cout << "�糡�ֲ�"<< data.E << std::endl;



}


public:	  void show_data() {

		  std::cout << "�缫��ȣ�" << ElecThickness << std::endl;
		  std::cout << "�缫������" << ElecNumber << std::endl;
		  std::cout << "���ڵ缫֮��ľ��룺" << std::endl;
		  for (int i = 0; i < ElecNumber; i++) {
			  std::cout << DisBetweenAdjacentEle[i] << std::endl;
		  }


		  std::cout << "���ڵ缫���ֵĲ�����" << std::endl;
		  for (int i = 0; i < ElecNumber - 1; i++) {
			  std::cout <<StepBetweenAdjacentEle[i] << std::endl;
		  }

		  std::cout << *StepBetweenAdjacentEle << std::endl;

		  std::cout << "�缫��λ��" << std::endl;
		  for (int i = 0; i < ElecNumber; i++) {
			  std::cout << ElecPotential[i] << std::endl;
		  }


		  std::cout << "�缫�ڿ׾��뾶��" << RadiusAperture << std::endl;
		  std::cout << "�缫�ڿ׾��뾶�Ȳ������ֵ���������" << GridAperture << std::endl;
		  std::cout << "�ӵ缫�ڿ��ص���ձ߽紦�ľ�����룺" << RadiusBoundary << std::endl;
		  std::cout << "�ӵ缫�ڿ��ص���ձ߽紦���ֵ�������:" << GridBoundary << std::endl;
		  std::cout << "��������:" << IterationAccuracy << std::endl;
		  std::cout << "�����ӡ�ռ��λʱ���������:" << NST << std::endl;
		  std::cout << "Ҫ��ɨ���λ�߼��ֵ��" << dov << std::endl;

		/*  std::cout << "Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ:" << std::endl;
		  for (int i = 0; i < ElecNumber; i++) {
			  std::cout << EOE[i] << std::endl;
		  }*/







	  }
};


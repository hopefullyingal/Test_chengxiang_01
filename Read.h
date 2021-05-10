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
	FILE* ParameterFile;
	FILE* OutputFile;

public:
	struct data{
	public :

		double mElecThickness;//�缫���
		int mElecNumber;//�缫����
		double* mDisBetweenAdjacentEle;//���ڵ缫֮��ľ���
		double* mStepBetweenAdjacentEle;//���ڵ缫���ֵĲ���
		double* mElecPotential;//�缫��λ
		double mRadiusAperture;//�缫�ڿ׾��뾶
		int mGridAperture;//�缫�ڿ׾��뾶�Ȳ������ֵ�������
		double mRadiusBoundary;//�ӵ缫�ڿ��ص���ձ߽紦�ľ������
		int mGridBoundary;//�ӵ缫�ڿ��ص���ձ߽紦���ֵ�������
		double mIterationAccuracy;//��������
		int mNST;//�����ӡ�ռ��λʱ���������
		double* mEOE;//Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ
		double mdov;//Ҫ��ɨ���λ�߼��ֵ


	};


	data getValue() {
		/*��������*/
		fopen_s(&ParameterFile, "E:\\Learning(School)\\������ԭ���뼼��\\������ʵ��\\����\\1120160852.dat", "r+");//��ȡ��ز����ļ�
		fopen_s(&OutputFile, "E:\\Learning(School)\\������ԭ���뼼��\\������ʵ��\\����\\11.res", "w+");
		/// <summary>����
		/// errno_t err = fopen_s(&ParameterFile, "E:\\Learning(School)\\������ԭ���뼼��\\������ʵ��\\����\\1120160852.dat", "r+");//��ȡ��ز����ļ�
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

		//while (getline(in, s))//���ж�ȡ���ݲ�����s�У�ֱ������ȫ����ȡ
		//{
		//	 str[hh] = UTF8ToGB(s.c_str()).c_str();
		//	 hh++;
		//	
		//}
		//
			/*TODO: ���鲻�ܿ���ElecNumber ���ǳ���*/
		////////////////////////////////////////////////////////////////
		fscanf_s(ParameterFile, "�缫������%d��\n", &ElecNumber);
		fscanf_s(ParameterFile, "�缫��ȣ�%lf��\n", &ElecThickness);

		DisBetweenAdjacentEle = (double*)calloc(ElecNumber, sizeof(double));//����缫���������С
		fscanf_s(ParameterFile, "���ڵ缫֮��ľ��룺%lf��", &DisBetweenAdjacentEle[0]);//ɨ����������һ��Ԫ��		
		for (i = 1; i < ElecNumber-1; i++)
		{
			fscanf_s(ParameterFile, "%lf��", &DisBetweenAdjacentEle[i]);

		}//�ѵڶ����������ڶ�������Ԫ�ض�ɨ�赽������
		fscanf_s(ParameterFile, "%lf��\n", &DisBetweenAdjacentEle[ElecNumber-1]);//�����һ������Ԫ��ɨ�赽������
		

		StepBetweenAdjacentEle = (double*)calloc(ElecNumber, sizeof(double));//��������������ڴ�ռ�
		fscanf_s(ParameterFile, "���ڵ缫���ֵĲ�����%lf��", &StepBetweenAdjacentEle[0]);//ɨ�貽�������һ��Ԫ��
		for (i = 1; i < ElecNumber-1; i++)
		{
			fscanf_s(ParameterFile, "%lf��", &StepBetweenAdjacentEle[i]);
			
		}//�Ѵӵڶ����������ڶ�������Ԫ��ɨ�赽������
		fscanf_s(ParameterFile, "%lf��\n", &StepBetweenAdjacentEle[ElecNumber-1]);//�����һ������Ԫ��ɨ�赽������



		ElecPotential = (double*)calloc(ElecNumber, sizeof(double));//Allocates memory space to the electronic potential array.
		fscanf_s(ParameterFile, "�缫��λ��%lf��", &ElecPotential[0]);
		for (i = 1; i < ElecNumber-1; i++)
		{
			fscanf_s(ParameterFile, "%lf��", &ElecPotential[i]);
		}
		fscanf_s(ParameterFile, "%lf��\n", &ElecPotential[ElecNumber-1]);//Put the last element into the array
	
		fscanf_s(ParameterFile, "�缫�ڿ׾��뾶��%lf��\n", &RadiusAperture);//Scan the  radius of the inner aperture of the electrode
		fscanf_s(ParameterFile, "�缫�ڿ׾��뾶�Ȳ������ֵ���������%d��\n", &GridAperture);//ɨ��缫�ڿ׾��뾶�Ȳ������ֵ�������
		fscanf_s(ParameterFile, "�ӵ缫�ڿ��ص���ձ߽紦�ľ�����룺%lf��\n", &RadiusBoundary);//ɨ��ӵ缫�ڿ��ص���ձ߽紦�ľ������
		fscanf_s(ParameterFile, "�ӵ缫�ڿ��ص���ձ߽紦���ֵ���������%d��\n", &GridBoundary);//ɨ��ӵ缫�ڿ��ص���ձ߽紦���ֵ�������
		fscanf_s(ParameterFile, "�������ȣ�%lf��\n", &IterationAccuracy);//ɨ���������
		fscanf_s(ParameterFile, "�����ӡ�ռ��λʱ�����������%d��\n", &NST);//ɨ�������ӡ�ռ��λʱ���������
		fscanf_s(ParameterFile, "Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ���λֵ:%lf��", &dov);//ɨ�����ϵ�λ���Ⱦ��ֵʱ�Ĳ�������



		//TODO EOE ��dov �ظ�		
		EOE = (double*)calloc(100, sizeof(double));
		//dov = 0;
		i = -1;
		do
		{
			i++;
			fscanf_s(ParameterFile, "%lf��", &EOE[i]);

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

		std::cout << "�缫��ȣ�" << Mdata.mElecThickness << std::endl;
		std::cout << "�缫������" << Mdata.mElecNumber << std::endl;
		std::cout << "���ڵ缫֮��ľ��룺" << std::endl;
		for (int i = 0; i < ElecNumber; i++) {
			std::cout<< Mdata.mDisBetweenAdjacentEle[i] << std::endl;
		}


		std::cout << "���ڵ缫���ֵĲ�����" << std::endl;
		for (int i = 0; i < ElecNumber; i++) {
			std::cout << Mdata.mStepBetweenAdjacentEle[i] << std::endl;
		}
		
		std::cout << "�缫��λ��" << std::endl;
		for (int i = 0; i < ElecNumber; i++) {
			std::cout << Mdata.mElecPotential[i] << std::endl;
		}

	
		std::cout << "�缫�ڿ׾��뾶��" << Mdata.mRadiusAperture << std::endl;
		std::cout << "�缫�ڿ׾��뾶�Ȳ������ֵ���������" << Mdata.mGridAperture << std::endl;
		std::cout << "�ӵ缫�ڿ��ص���ձ߽紦�ľ�����룺" << Mdata.mRadiusBoundary << std::endl;
		std::cout << "�ӵ缫�ڿ��ص���ձ߽紦���ֵ�������:" << Mdata.mGridBoundary << std::endl;
		std::cout << "��������:" << Mdata.mIterationAccuracy << std::endl;
		std::cout << "�����ӡ�ռ��λʱ���������:" << Mdata.mNST << std::endl;
		std::cout << "Ҫ��ɨ���λ�߼��ֵ��" << Mdata.mdov << std::endl;

		std::cout << "Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ:" << std::endl;
		for (int i = 0; i < ElecNumber; i++) {
			std::cout << Mdata.mEOE[i] << std::endl;
		}

		





	}


	/*�����������������ʱ���������*/
	std::string UTF8ToGB(const char* str)
	{
		std::string result;
		WCHAR* strSrc;
		LPSTR szRes;

		//�����ʱ�����Ĵ�С
		int i = MultiByteToWideChar(CP_UTF8, 0, str, -1, NULL, 0);
		strSrc = new WCHAR[i + 1];
		MultiByteToWideChar(CP_UTF8, 0, str, -1, strSrc, i);

		//�����ʱ�����Ĵ�С
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


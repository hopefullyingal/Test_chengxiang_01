//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <windows.h>
//#include <graphics.h>
//#include <conio.h>
//#include <iostream>
//#include <cstdio>
//
//void test();
//void SOR();
//void setfree();
///////////////
///*����ȫ�ֱ���*/
//
///*����������*/
//double ET;//�缫���
//int n;//�缫����
//double* z;//���ڵ缫֮��ľ���
//double* N;//���ڵ缫���ֵĲ���
//double* V;//�缫��λ
//double r1;//�缫�ڿ׾��뾶
//int M1;//�缫�ڿ׾��뾶�Ȳ������ֵ�������
//double r2;//�ӵ缫�ڿ��ص���ձ߽紦�ľ������
//int M2;//�ӵ缫�ڿ��ص���ձ߽紦���ֵ�������
//double e;//��������
//int NST;//�����ӡ�ռ��λʱ���������
//double* EOE;//Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ
//double dov;//Ҫ��ɨ���λ�߼��ֵ
//
///*����������*/
//double* E;//�糡�ֲ�
//int row, lin;//�У���
//int* pon;//�缫λ�ö�Ӧ�������
//double* loz;//z�᷽�����粽������
//double* poz;//z����λ������
//double* pozE;//z��ȵ�λ���λ������
//double* lor;//r�᷽�����粽������
//double* por;//r����λ������
//double* Tpor;//ת����r����λ������
//int* nopoE;//��λ�����
//double* porE;//r��ȵ�λ���λ������
//double* TporE;//ת����r��ȵ�λ���λ������
//FILE* ParameterFile;
//FILE* OutputFile;
//double c1, c2, c3, c4, c0;    //ϵ��
//double E2;                //������ĸ���λ
//int turn;                    //�����ִ�
//int times;                 //��������
//int atime;              //�����ܴ���
//double w;                 //��������
//double* dV;              //����λ�в�
//double sdV;               //����λ�в��
//double mdV;               //����λ���в�
//double adV;               //����λƽ���в�
//double adV11;              //�ڶ���12�ε�������λƽ���в�
//double adV12;              //�ڶ��ֵ�13�ε�������λƽ���в�
//double Lambda;             //�в��ֵ��ֵ
//double uLambda;           //�м����
//double wLambda;
//double wm;                //�������w
//double wm1;               //�������w
//int    judge;             //�ж���
//
//
//int main()
//{
//	int i, j, k, s, t;
//
//	/*��������*/
//	errno_t err = fopen_s(&ParameterFile, "E:\\Learning(School)\\������ԭ���뼼��\\������ʵ��\\����\\1120160852.dat", "r+");//��ȡ��ز����ļ�
//	fopen_s(&OutputFile, "E:\\Learning(School)\\������ԭ���뼼��\\������ʵ��\\����\\11.res", "w+");
//	if (ParameterFile == NULL)
//	{
//		printf("Open filefailure!");
//		exit(1);
//	}
//	fscanf_s(ParameterFile, "�缫������%d��\n", &n);
//	fscanf_s(ParameterFile, "�缫��ȣ�%lf��\n", &ET);
//	z = (double*)calloc(n, sizeof(double));//����缫���������С
//	fscanf_s(ParameterFile, "���ڵ缫֮��ľ��룺%lf��", &z[0]);//ɨ����������һ��Ԫ��
//	n = n - 1;//ʹ���Ʒ�Χ�������ڶ���Ԫ��
//	for (i = 1; i < n; i++)
//	{
//		fscanf_s(ParameterFile, "%lf��", &z[i]);
//
//	}//�ѵڶ����������ڶ�������Ԫ�ض�ɨ�赽������
//	fscanf_s(ParameterFile, "%lf��\n", &z[n]);//�����һ������Ԫ��ɨ�赽������
//	n = n + 1;//��ԭ��������Ԫ�ظ���
//	N = (double*)calloc(n, sizeof(double));//��������������ڴ�ռ�
//	fscanf_s(ParameterFile, "���ڵ缫���ֵĲ�����%lf��", &N[0]);//ɨ�貽�������һ��Ԫ��
//	n = n - 1;//ʹ���Ʒ�Χ�������ڶ���Ԫ��
//	for (i = 1; i < n; i++)
//	{
//		fscanf_s(ParameterFile, "%lf��", &N[i]);
//	}//�Ѵӵڶ����������ڶ�������Ԫ��ɨ�赽������
//	fscanf_s(ParameterFile, "%lf��\n", &N[n]);//�����һ������Ԫ��ɨ�赽������
//	n = n + 1;//��ԭ�缫����
//	V = (double*)calloc(n, sizeof(double));//Allocates memory space to the electronic potential array.
//	fscanf_s(ParameterFile, "�缫��λ��%lf��", &V[0]);
//	n = n - 1;//ʹ���Ʒ�Χ�������ڶ���Ԫ��
//	for (i = 1; i < n; i++)
//	{
//		fscanf_s(ParameterFile, "%lf��", &V[i]);
//	}
//	fscanf_s(ParameterFile, "%lf��\n", &V[n]);//Put the last element into the array
//	n = n + 1;//��ԭ�缫����
//	fscanf_s(ParameterFile, "�缫�ڿ׾��뾶��%lf��\n", &r1);//Scan the  radius of the inner aperture of the electrode
//	fscanf_s(ParameterFile, "�缫�ڿ׾��뾶�Ȳ������ֵ���������%d��\n", &M1);//ɨ��缫�ڿ׾��뾶�Ȳ������ֵ�������
//	fscanf_s(ParameterFile, "�ӵ缫�ڿ��ص���ձ߽紦�ľ�����룺%lf��\n", &r2);//ɨ��ӵ缫�ڿ��ص���ձ߽紦�ľ������
//	fscanf_s(ParameterFile, "�ӵ缫�ڿ��ص���ձ߽紦���ֵ���������%d��\n", &M2);//ɨ��ӵ缫�ڿ��ص���ձ߽紦���ֵ�������
//	fscanf_s(ParameterFile, "�������ȣ�%lf��\n", &e);//ɨ���������
//	fscanf_s(ParameterFile, "�����ӡ�ռ��λʱ�����������%d��\n", &NST);//ɨ�������ӡ�ռ��λʱ���������
//	fscanf_s(ParameterFile, "Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ���λֵ:", &dov);
//	EOE = (double*)calloc(100, sizeof(double));
//	dov = 0;
//	i = -1;
//	do
//	{
//		i++;
//		fscanf_s(ParameterFile, "%lf��", &EOE[i]);
//
//	} while (EOE[i] != 0);
//	if (i == 1)
//		dov = EOE[0];
//
//
//	//fscanf_s(ParameterFile, "Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ���λֵ:%lf��\n", &dov);//ɨ�����ϵ�λ���Ⱦ��ֵʱ�Ĳ�������
//	fclose(ParameterFile);
//
//	for (i = 1; i < n; i++) {
//
//		std::cout << V[i] << std::endl;
//		
//
//	}
//
//	system("pause");
//	//
//	//	//////////////////////////////////////
//	//		/*�糡��λ��ʼ��*/
//	//	row = M1 + M2 + 1;//Calculate the value of the row.
//	//	lin = 0;//Assign a initial value to line so that it's final value can be calculate by a circulation
//	//	for (i = 0; i < n; i++)
//	//	{
//	//		lin = lin + *(N + i);
//	//	}
//	//	lin = lin + n;
//	//	E = (double*)calloc(lin * row, sizeof(double));
//	//	s = 0;//�����Ͻ�
//	//	t = *N + 1;//�����½�
//	//	for (i = 0; i < row; i++)//��ѭ��
//	//	{
//	//		s = 0;//�����Ͻ�
//	//		t = *N + 1;//�����½�
//	//		if (i == row - 1)//��һ������
//	//		{
//	//			for (j = 0; j < n; j++)//�缫��λ
//	//			{
//	//				if (j == 0)
//	//				{
//	//					int f = 1;//��������ѭ���Ĵ���
//	//
//	//					for (k = s; k < t; k++)
//	//					{
//	//
//	//						double d = *(V + j) / *(N + j);//���㹫���Ƶ����漴��
//	//						*(E + i * lin + k) = d * (f - 1);
//	//						f = f + 1;
//	//					}
//	//					*(E + i * lin + k - 1) = *(V + j);
//	//					*(E + i * lin + k) = *(V + j);
//	//					int o = *(N + j + 1) + 1;
//	//					s = t + 1;
//	//					t = t + o;
//	//				}
//	//				if (j > 0)
//	//				{
//	//					int f = 1;//��������ѭ���Ĵ���
//	//
//	//					for (k = s; k < t; k++)
//	//					{
//	//
//	//						double d = (*(V + j) - *(V + j - 1)) / *(N + j);//���㹫���Ƶ����漴��
//	//						*(E + i * lin + k) = *(V + j - 1) + d * f;
//	//						f = f + 1;
//	//					}
//	//					if (k < lin)
//	//					{
//	//						*(E + i * lin + k - 1) = *(V + j);
//	//						*(E + i * lin + k) = *(V + j);
//	//						int o = *(N + j + 1) + 1;
//	//						s = t + 1;
//	//						t = t + o;
//	//					}
//	//					if (k == lin)
//	//					{
//	//						*(E + i * lin + lin - 1) = *(V + j);
//	//					}
//	//
//	//				}
//	//
//	//
//	//			}
//	//		}
//	//		if (i < row - 1 && i >= M1)
//	//		{
//	//			for (j = 0; j < n; j++)//�缫��λ
//	//			{
//	//
//	//
//	//				for (k = s; k < t; k++)
//	//				{
//	//
//	//
//	//					*(E + i * lin + k) = 0;
//	//
//	//				}
//	//				//���ж�����ܹ���
//	//				if (k < lin)
//	//				{
//	//					*(E + i * lin + k - 1) = *(V + j);
//	//					*(E + i * lin + k) = *(V + j);
//	//					int o = *(N + j + 1) + 1;
//	//					s = t + 1;
//	//					t = t + o;
//	//				}
//	//				if (k == lin)
//	//				{
//	//					*(E + i * lin + lin - 1) = *(V + j);
//	//				}
//	//			}
//	//		}
//	//		if (i < M1)
//	//		{
//	//			for (j = 0; j < lin; j++)
//	//			{
//	//				*(E + i * lin + j) = 0;
//	//			}
//	//			*(E + i * lin + lin - 1) = *(V + n - 1);
//	//		}
//	//	}
//	//	/*���缫λ�ö�Ӧ������긳ֵ*/
//	//	pon = (int*)calloc(n + 1, sizeof(int));//���缫λ�ö�Ӧ�����������ڴ�ռ�
//	//	*pon = 0;
//	//	for (i = 1; i < n + 1; i++)
//	//	{
//	//		if (i == n)
//	//		{
//	//			*(pon + i) = lin - 1;//ӫ�������λ��
//	//		}
//	//		else
//	//		{
//	//			*(pon + i) = *(pon + i - 1) + *(N + i - 1) + 1;//�缫�Ҹ��λ��
//	//		}
//	//	}
//	//	/*��z�᷽�����粽�����ָ�ֵ*/
//	//	loz = (double*)calloc((lin - 1) * (row - 1), sizeof(double));//��z�᷽�����绮�ֲ��������ڴ�ռ�
//	//	for (i = 0; i < row - 1; i++)
//	//	{
//	//		for (j = 1; j < n + 1; j++)
//	//		{
//	//			if (j != n)
//	//			{
//	//				for (k = *(pon + j - 1); k < *(pon + j); k++)
//	//					*(loz + i * (lin - 1) + k) = *(z + j - 1) / *(N + j - 1);
//	//				int a = *(pon + j);
//	//				*(loz + i * (lin - 1) + a - 1) = ET;
//	//			}
//	//			if (j == n)
//	//			{
//	//				for (k = *(pon + j - 1); k < *(pon + j); k++)
//	//					*(loz + i * (lin - 1) + k) = *(z + j - 1) / *(N + j - 1);
//	//			}
//	//		}
//	//	}
//	//	/*��z����λ�����긳ֵ*/
//	//	poz = (double*)calloc(lin, sizeof(double));//��z����λ����������ڴ�ռ�
//	//	*poz = 0;
//	//	for (i = 1; i < lin; i++)
//	//		*(poz + i) = *(poz + i - 1) + *(loz + i - 1);
//	//	/*��r�᷽�����粽�����ָ�ֵ*/
//	//	lor = (double*)calloc((lin - 1) * (row - 1), sizeof(double));//��lor�����ڴ�ռ�
//	//	for (i = 0; i < row - 1; i++)
//	//	{
//	//		if (i < M1)
//	//		{
//	//			for (j = 0; j < lin - 1; j++)
//	//				*(lor + i * (lin - 1) + j) = r1 / M1;
//	//		}
//	//		else
//	//		{
//	//			for (j = 0; j < lin - 1; j++)
//	//				*(lor + i * (lin - 1) + j) = r2 / M2;
//	//		}
//	//
//	//	}
//	//	/*��r����λ�ø�ֵ*/
//	//	por = (double*)calloc(row, sizeof(double));//
//	//	*por = 0;
//	//	for (i = 1; i < row; i++)
//	//		*(por + i) = *(por + i - 1) + *(lor + (i - 1) * (lin - 1));
//	//
//	//	/////////////////////////////////////////////	
//	//		/*��ʼ������ܵ�λ�ֲ�*/
//	//		/*��һ�ֵ���*/
//	//	turn = 1;
//	//	times = 1;
//	//	atime = 0;
//	//	w = 1;
//	//	SOR();
//	//	atime = atime + 1;
//	//
//	//	/*�ڶ��ֵ���*/
//	//	w = 1.375;
//	//	turn = turn + 1;
//	//	for (times = 1; times < 13; times++)
//	//	{
//	//		SOR();
//	//		atime = atime + 1;
//	//		if (times == 11)
//	//			adV11 = adV;
//	//		if (times == 12)
//	//			adV12 = adV;
//	//	}
//	//	Lambda = adV12 / adV11;
//	//	uLambda = (Lambda + w - 1) / (sqrt(Lambda) * w);
//	//	wLambda = 2 / (1 + sqrt(1 - uLambda * uLambda));
//	//	wm = 1.25 * wLambda - 0.5;
//	//	/*�����������Թ̶���������*/
//	//	do
//	//	{
//	//		w = wm;
//	//		turn++;
//	//		for (times = 1; times < 13; times++)
//	//		{
//	//			SOR();
//	//			atime = atime + 1;
//	//			if (times == 11)
//	//				adV11 == adV;
//	//			if (times == 12)
//	//				adV12 == adV;
//	//		}
//	//		wm1 = wm;
//	//		Lambda = adV12 / adV11;
//	//		uLambda = (Lambda + w - 1) / (sqrt(Lambda) * w);
//	//		wLambda = 2 / (1 + sqrt(1 - uLambda * uLambda));
//	//		wm = 1.25 * wLambda - 0.5;
//	//	} while (fabs((wm - wm1) / (2 - wm1)) >= 0.05);
//	//	w = wm;
//	//	/*�������ֵ����Դﵽ����Ҫ��*/
//	//	do
//	//	{
//	//		turn++;
//	//		SOR();   //��ѵ�������ȷ�����������
//	//		atime++;
//	//		judge = 0;
//	//		for (i = 0; i < row; i++)
//	//			for (j = 0; j < lin; j++)
//	//			{
//	//				if (*(dV + i * lin + j) >= e)
//	//					judge++;
//	//			}
//	//	} while (judge > 0);
//	//	for (i = row - 1; i >= 0; i--)
//	//	{
//	//		for (j = 0; j < lin - 1; j++)
//	//		{
//	//			fprintf(OutputFile, "%3.0lf", *(E + i * lin + j));
//	//		}
//	//		fprintf(OutputFile, "%3.0lf\n", *(E + i * lin + lin - 1));
//	//	}
//	//
//	//	////////////////////////////////////////////
//	//	   /*ɨ���λ������*/
//	//	double a = 0;
//	//	int m = 1;
//	//	if (dov != 0)
//	//	{
//	//		for (m = 1; a < (V[n - 1] - dov); m++)   //�����λ����
//	//			a = a + dov;
//	//		/*EOE = (double*)calloc(m - 1, sizeof(double));*/
//	//		*EOE = dov;
//	//		for (i = 1; i < m - 1; i++)
//	//			EOE[i] = EOE[i - 1] + dov;
//	//	}
//	//	else if (dov == 0)
//	//	{
//	//		m = 0;
//	//		for (i = 0; i < 100; i++)
//	//		{
//	//			if (EOE[i] != 0)
//	//				m++;
//	//
//	//		}
//	//		m = m + 1;
//	//		/*EOE = (double*)calloc(m - 1, sizeof(double));*/
//	//	}
//	//	//for (m = 1; a < (V[n - 1] - dov); m++)   //�����λ����
//	//	//	a = a + dov;
//	//	//EOE = (double*)calloc(m - 1, sizeof(double));
//	//	porE = (double*)calloc((m - 1) * 2 * row, sizeof(double));
//	//	TporE = (double*)calloc((m - 1) * 2 * row, sizeof(double));
//	//	pozE = (double*)calloc(((m - 1)) * 2 * lin, sizeof(double));
//	//	/**EOE = dov;
//	//	for (i = 1; i < m - 1; i++)
//	//		EOE[i] = EOE[i - 1] + dov;*/
//	//	nopoE = (int*)calloc(m - 1, sizeof(int));
//	//	for (i = 0; i < m - 1; i++)
//	//		nopoE[i] = 0;
//	//	/*int h = 0;*/
//	//	for (j = 1; j < lin - 1; j++)
//	//	{
//	//		for (k = 0; k < m - 1; k++)
//	//		{
//	//			for (i = row - 1; i > 0; i--)
//	//			{
//	//				if (((*(E + i * lin + j) < EOE[k]) & (*(E + i * lin + j - 1) > EOE[k])) || ((*(E + i * lin + j) > EOE[k]) & (*(E + i * lin + j - 1) < EOE[k])))
//	//				{
//	//					porE[k * 2 * row + nopoE[k]] = por[i];
//	//					pozE[k * 2 * lin + nopoE[k]] = poz[j] + ((poz[j - 1] - poz[j]) / (*(E + i * lin + j - 1) - *(E + i * lin + j))) * (EOE[k] - *(E + i * lin + j));
//	//					nopoE[k]++;
//	//				}
//	//				if (*(E + i * lin + j) == EOE[k])
//	//				{
//	//					porE[k * 2 * row + nopoE[k]] = por[i];
//	//					pozE[k * 2 * lin + nopoE[k]] = poz[j];
//	//					nopoE[k]++;
//	//				}
//	//				if (((*(E + i * lin + j) < EOE[k]) & (*(E + (i - 1) * lin + j) > EOE[k])) || ((*(E + i * lin + j) > EOE[k]) & (*(E + (i - 1) * lin + j) < EOE[k])))
//	//				{
//	//					pozE[k * 2 * lin + nopoE[k]] = poz[j];
//	//					porE[k * 2 * row + nopoE[k]] = por[i] + ((por[i - 1] - por[i]) / (*(E + (i - 1) * lin + j) - *(E + i * lin + j))) * (EOE[k] - *(E + i * lin + j));
//	//					nopoE[k]++;
//	//				}
//	//			}
//	//			i = 0;
//	//			if (((*(E + i * lin + j) < EOE[k]) & (*(E + i * lin + j - 1) > EOE[k])) || ((*(E + i * lin + j) > EOE[k]) & (*(E + i * lin + j - 1) < EOE[k])))
//	//			{
//	//				porE[k * 2 * row + nopoE[k]] = por[i];
//	//				pozE[k * 2 * lin + nopoE[k]] = poz[j] + ((poz[j - 1] - poz[j]) / (*(E + i * lin + j - 1) - *(E + i * lin + j))) * (EOE[k] - *(E + i * lin + j));
//	//				nopoE[k]++;
//	//			}
//	//			if (*(E + i * lin + j) == EOE[k])
//	//			{
//	//				porE[k * 2 * row + nopoE[k]] = por[i];
//	//				pozE[k * 2 * lin + nopoE[k]] = poz[j];
//	//				nopoE[k]++;
//	//			}
//	//		}
//	//	}
//	//	for (k = 0; k < m - 1; k++)
//	//	{
//	//		fprintf(OutputFile, "��%d��λ,��λֵ��%lf,\n", k + 1, EOE[k]);
//	//		for (j = 0; j < nopoE[k] - 1; j++)
//	//		{
//	//			fprintf(OutputFile, "(%3.2lf,%3.2lf)\t", pozE[k * 2 * lin + j], porE[k * 2 * row + j]);
//	//		}
//	//		fprintf(OutputFile, "(%3.2lf,%3.2lf)\t\n", pozE[k * 2 * lin + nopoE[k] - 1], porE[k * 2 * row + nopoE[k] - 1]);
//	//	}//�����λ������
//	//
//	//	 ////////Ϊ������ͼ����ʼ����任////////////
//	//	/*���ڸó��򽫹����������潻����Ϊԭ�㣬����ͼ������ǽ����ϽǶ�λԭ�㣬���Ϊ�˹۲�ϰ�ߣ���r�᷽��������б任*/
//	//	/*Tpor = (double*)calloc(row, sizeof(double));
//	//	for (i = 0; i < row; i++)
//	//		Tpor[i] = r1 + r2 - por[i];*/
//	//	for (k = 0; k < m - 1; k++)
//	//	{
//	//		for (j = 0; j < nopoE[k]; j++)
//	//		{
//	//			TporE[k * 2 * row + j] = r1 + r2 - porE[k * 2 * row + j];
//	//
//	//		}
//	//	}
//	//	/*for (j = 0; j < nopoE[0]; j++)
//	//		printf("%lf\t", porE[j]);
//	//	printf("\n");
//	//	for (j = 0; j < nopoE[0]; j++)
//	//	printf("%lf\t",	TporE[ j] );*/
//	//	////////////////////////////////////////////
//	//		 /*��ʼ���Ƶ�λ��*/
//	//	initgraph(110 + int(10 * poz[lin - 1]), 110 + int(10 * por[row - 1]));
//	//	setbkcolor(BLACK);
//	//	cleardevice();
//	//	rectangle(50, 50, 50 + int(10 * poz[lin - 1]), 50 + int(10 * por[row - 1]));
//	//	for (i = 1; i < n; i++)
//	//		bar(50 + (10 * poz[pon[i] - 1]), 50, 50 + (10 * (poz[pon[i]])), 50 + (10 * r2));
//	//
//	//	for (k = 0; k < m - 1; k++)
//	//	{
//	//		for (j = 0; j < nopoE[k] - 1; j++)
//	//			line(50 + 10 * pozE[k * 2 * lin + j], 50 + 10 * (TporE[k * 2 * row + j]), 50 + 10 * pozE[k * 2 * lin + j + 1], 50 + 10 * (TporE[k * 2 * row + j + 1]));
//	//		int R = rand() % 255, G = rand() % 255, B = rand() % 255;
//	//		setlinecolor(RGB(R, G, B));
//	//	}
//	//	_getch();
//	//	closegraph();
//	//	_getch();
//	//	setfree();
//	//}
//	//
//	//
//	//
//	///*����������й����Ӻ���*/
//	//void test()
//	//{
//	//	int i, j, k, s, t;
//	//	/*printf("%lf\n", v);*/
//	//	printf("%d\n", lin);
//	//	for (i = row - 1; i >= 0; i++)
//	//	{
//	//		for (j = 0; j < lin - 1; j++)
//	//		{
//	//			printf("%2.0lf ", *(E + i * lin + j));
//	//		}
//	//		printf("%3.0lf\n", *(E + i * lin + lin - 1));
//	//	}
//	//	for (i = 0; i < n + 1; i++)
//	//		printf("%lf,", *(pon + i));
//	//	for (i = 0; i < row - 1; i++)
//	//	{
//	//		for (j = 0; j < lin - 2; j++)
//	//			printf("%2.4lf ", *(loz + i * (lin - 1) + j));
//	//		printf("%2.4lf\n", *(loz + i * (lin - 1) + lin - 2));
//	//	}
//	//	for (i = 0; i < lin; i++)
//	//		printf("%2.2lf, ", *(poz + i));
//	//	for (i = 0; i < row - 1; i++)
//	//	{
//	//		for (j = 0; j < lin - 2; j++)
//	//		{
//	//			printf("%2.2lf ", *(lor + i * (lin - 1) + j));
//	//		}
//	//		printf("%2.2lf\n", *(lor + i * (lin - 1) + lin - 2));
//	//	}
//	//	for (i = 1; i < row; i++)
//	//		printf("%2.2lf\n", *(por + i));
//	//}
//	//
//	///*���ųڵ������Ӻ���*/
//	//void SOR()
//	//{
//	//	dV = (double*)calloc(row * lin, sizeof(double));
//	//	for (int i = 0; i < row - 1; i++)
//	//	{
//	//		if (i == 0)
//	//		{
//	//
//	//			for (int j = 1; j < lin - 1; j++)
//	//			{
//	//				c1 = 2 / ((*(loz + i * (lin - 1) + j - 1)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));   //ϵ������
//	//				c2 = 2 / ((*(loz + i * (lin - 1) + j)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));
//	//				c3 = 0;
//	//				c4 = 4 / ((*(lor + i * (lin - 1) + j)) * (*(lor + i * (lin - 1) + j)));
//	//				c0 = c1 + c2 + c3 + c4;
//	//				E2 = (1 - w) * (*(E + i * lin + j)) + w * (c1 * (*(E + i * lin + j - 1)) + c2 * (*(E + i * lin + j + 1)) + c3 * (*(E + (i - 1) * lin + j)) + c4 * (*(E + (i + 1) * lin + j))) / c0;
//	//				*(dV + i * lin + j) = fabs(E2 - *(E + i * lin + j));
//	//				*(E + i * lin + j) = E2;
//	//
//	//			}
//	//
//	//		}
//	//		if (i < M1 && i>0)
//	//		{
//	//			for (int j = 1; j < lin - 1; j++)
//	//			{
//	//				c1 = 2 / ((*(loz + i * (lin - 1) + j - 1)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));   //ϵ������
//	//				c2 = 2 / ((*(loz + i * (lin - 1) + j)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));
//	//				c3 = (2 * (*(por + i)) - *(lor + i * (lin - 1) + j)) / ((*(por + i)) * (*(lor + (i - 1) * (lin - 1) + j)) * (*(lor + (i - 1) * (lin - 1) + j) + (*(lor + i * (lin - 1) + j))));
//	//				c4 = (2 * (*(por + i)) + (*(lor + (i - 1) * (lin - 1) + j))) / ((*(por + i)) * (*(lor + i * (lin - 1) + j)) * (*(lor + (i - 1) * (lin - 1) + j) + *(lor + i * (lin - 1) + j)));
//	//				c0 = c1 + c2 + c3 + c4;
//	//				E2 = (1 - w) * (*(E + i * lin + j)) + w * (c1 * (*(E + i * lin + j - 1)) + c2 * (*(E + i * lin + j + 1)) + c3 * (*(E + (i - 1) * lin + j)) + c4 * (*(E + (i + 1) * lin + j))) / c0;
//	//				*(dV + i * lin + j) = fabs(E2 - *(E + i * lin + j));
//	//				*(E + i * lin + j) = E2;
//	//			}
//	//		}
//	//		if (i >= M1)
//	//		{
//	//			for (int k = 1; k < n + 1; k++)//�缫ѭ��
//	//			{
//	//				if (k != n)
//	//				{
//	//					for (int j = *(pon + k - 1) + 1; j < *(pon + k) - 1; j++)//����ѭ������ܱ߽��
//	//					{
//	//						c1 = 2 / ((*(loz + i * (lin - 1) + j - 1)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));   //ϵ������
//	//						c2 = 2 / ((*(loz + i * (lin - 1) + j)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));
//	//						c3 = (2 * (*(por + i)) - *(lor + i * (lin - 1) + j)) / ((*(por + i)) * (*(lor + (i - 1) * (lin - 1) + j)) * (*(lor + (i - 1) * (lin - 1) + j) + (*(lor + i * (lin - 1) + j))));
//	//						c4 = (2 * (*(por + i)) + (*(lor + (i - 1) * (lin - 1) + j))) / ((*(por + i)) * (*(lor + i * (lin - 1) + j)) * (*(lor + (i - 1) * (lin - 1) + j) + *(lor + i * (lin - 1) + j)));
//	//						c0 = c1 + c2 + c3 + c4;
//	//						E2 = (1 - w) * (*(E + i * lin + j)) + w * (c1 * (*(E + i * lin + j - 1)) + c2 * (*(E + i * lin + j + 1)) + c3 * (*(E + (i - 1) * lin + j)) + c4 * (*(E + (i + 1) * lin + j))) / c0;
//	//						*(dV + i * lin + j) = fabs(E2 - *(E + i * lin + j));
//	//						*(E + i * lin + j) = E2;
//	//					}
//	//				}
//	//				else
//	//				{
//	//					for (int j = *(pon + k - 1) + 1; j < *(pon + k); j++)
//	//					{
//	//						c1 = 2 / ((*(loz + i * (lin - 1) + j - 1)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));   //ϵ������
//	//						c2 = 2 / ((*(loz + i * (lin - 1) + j)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));
//	//						c3 = (2 * (*(por + i)) - *(lor + i * (lin - 1) + j)) / ((*(por + i)) * (*(lor + (i - 1) * (lin - 1) + j)) * (*(lor + (i - 1) * (lin - 1) + j) + (*(lor + i * (lin - 1) + j))));
//	//						c4 = (2 * (*(por + i)) + (*(lor + (i - 1) * (lin - 1) + j))) / ((*(por + i)) * (*(lor + i * (lin - 1) + j)) * (*(lor + (i - 1) * (lin - 1) + j) + *(lor + i * (lin - 1) + j)));
//	//						c0 = c1 + c2 + c3 + c4;
//	//						E2 = (1 - w) * (*(E + i * lin + j)) + w * (c1 * (*(E + i * lin + j - 1)) + c2 * (*(E + i * lin + j + 1)) + c3 * (*(E + (i - 1) * lin + j)) + c4 * (*(E + (i + 1) * lin + j))) / c0;
//	//						*(dV + i * lin + j) = fabs(E2 - *(E + i * lin + j));
//	//						*(E + i * lin + j) = E2;
//	//					}
//	//				}
//	//			}
//	//
//	//		}
//	//	}
//	//	sdV = 0;
//	//	mdV = 0;
//	//	for (int i = 0; i < row * lin; i++)
//	//	{
//	//		sdV = sdV + *(dV + i);
//	//		if (*(dV + i) > mdV)
//	//			mdV = *(dV + i);
//	//	}
//	//	adV = sdV / (row * lin);
//	//}
//	//
//	//void setfree()
//	//{
//	//	free(z);
//	//	free(N);
//	//	free(V);
//	//	free(EOE);
//	//	free(E);
//	//	free(pon);
//	//	free(loz);
//	//	free(poz);
//	//	free(pozE);
//	//	free(lor);
//	//	free(por);
//	//	free(Tpor);
//	//	free(nopoE);
//	//	free(porE);
//	//	free(TporE);
//	//
//	//}
//}
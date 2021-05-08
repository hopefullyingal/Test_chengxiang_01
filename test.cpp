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
///*定义全局变量*/
//
///*定义读入参数*/
//double ET;//电极厚度
//int n;//电极总数
//double* z;//相邻电极之间的距离
//double* N;//相邻电极划分的步长
//double* V;//电极电位
//double r1;//电极内孔径半径
//int M1;//电极内孔径半径等步长划分的网格数
//double r2;//从电极内孔沿到封闭边界处的径向距离
//int M2;//从电极内孔沿到封闭边界处划分的网格数
//double e;//迭代精度
//int NST;//输出打印空间电位时网格点间隔数
//double* EOE;//要求扫描搜索等电位线的电位间隔值
//double dov;//要求扫描等位线间隔值
//
///*定义计算参数*/
//double* E;//电场分布
//int row, lin;//行，列
//int* pon;//电极位置对应格点坐标
//double* loz;//z轴方向网络步长划分
//double* poz;//z轴格点位置坐标
//double* pozE;//z轴等电位格点位置坐标
//double* lor;//r轴方向网络步长划分
//double* por;//r轴格点位置坐标
//double* Tpor;//转换后r轴格点位置坐标
//int* nopoE;//等位点个数
//double* porE;//r轴等电位格点位置坐标
//double* TporE;//转换后r轴等电位格点位置坐标
//FILE* ParameterFile;
//FILE* OutputFile;
//double c1, c2, c3, c4, c0;    //系数
//double E2;                //迭代后的格点电位
//int turn;                    //迭代轮次
//int times;                 //迭代次数
//int atime;              //迭代总次数
//double w;                 //迭代因子
//double* dV;              //格点电位残差
//double sdV;               //格点电位残差和
//double mdV;               //格点电位最大残差
//double adV;               //格点电位平均残差
//double adV11;              //第二轮12次迭代格点电位平均残差
//double adV12;              //第二轮第13次迭代格点电位平均残差
//double Lambda;             //残差均值比值
//double uLambda;           //中间变量
//double wLambda;
//double wm;                //修正后的w
//double wm1;               //修正后的w
//int    judge;             //判断子
//
//
//int main()
//{
//	int i, j, k, s, t;
//
//	/*读入数据*/
//	errno_t err = fopen_s(&ParameterFile, "E:\\Learning(School)\\光电成像原理与技术\\光电成像实验\\程序\\1120160852.dat", "r+");//读取相关参数文件
//	fopen_s(&OutputFile, "E:\\Learning(School)\\光电成像原理与技术\\光电成像实验\\程序\\11.res", "w+");
//	if (ParameterFile == NULL)
//	{
//		printf("Open filefailure!");
//		exit(1);
//	}
//	fscanf_s(ParameterFile, "电极总数：%d；\n", &n);
//	fscanf_s(ParameterFile, "电极厚度：%lf；\n", &ET);
//	z = (double*)calloc(n, sizeof(double));//分配电极距离数组大小
//	fscanf_s(ParameterFile, "相邻电极之间的距离：%lf，", &z[0]);//扫描距离数组第一个元素
//	n = n - 1;//使复制范围到倒数第二个元素
//	for (i = 1; i < n; i++)
//	{
//		fscanf_s(ParameterFile, "%lf，", &z[i]);
//
//	}//把第二个至倒数第二个距离元素都扫描到数组里
//	fscanf_s(ParameterFile, "%lf；\n", &z[n]);//把最后一个距离元素扫描到数组里
//	n = n + 1;//还原距离数组元素个数
//	N = (double*)calloc(n, sizeof(double));//给步长数组分配内存空间
//	fscanf_s(ParameterFile, "相邻电极划分的步长：%lf，", &N[0]);//扫描步长数组第一个元素
//	n = n - 1;//使复制范围到倒数第二个元素
//	for (i = 1; i < n; i++)
//	{
//		fscanf_s(ParameterFile, "%lf，", &N[i]);
//	}//把从第二个至倒数第二个步长元素扫描到数组里
//	fscanf_s(ParameterFile, "%lf；\n", &N[n]);//把最后一个数组元素扫描到数组里
//	n = n + 1;//还原电极个数
//	V = (double*)calloc(n, sizeof(double));//Allocates memory space to the electronic potential array.
//	fscanf_s(ParameterFile, "电极电位：%lf，", &V[0]);
//	n = n - 1;//使复制范围到倒数第二个元素
//	for (i = 1; i < n; i++)
//	{
//		fscanf_s(ParameterFile, "%lf，", &V[i]);
//	}
//	fscanf_s(ParameterFile, "%lf；\n", &V[n]);//Put the last element into the array
//	n = n + 1;//还原电极个数
//	fscanf_s(ParameterFile, "电极内孔径半径：%lf；\n", &r1);//Scan the  radius of the inner aperture of the electrode
//	fscanf_s(ParameterFile, "电极内孔径半径等步长划分的网格数：%d；\n", &M1);//扫描电极内孔径半径等步长划分的网格数
//	fscanf_s(ParameterFile, "从电极内孔沿到封闭边界处的径向距离：%lf；\n", &r2);//扫描从电极内孔沿到封闭边界处的径向距离
//	fscanf_s(ParameterFile, "从电极内孔沿到封闭边界处划分的网格数：%d；\n", &M2);//扫描从电极内孔沿到封闭边界处划分的网格数
//	fscanf_s(ParameterFile, "迭代精度：%lf；\n", &e);//扫描迭代精度
//	fscanf_s(ParameterFile, "输出打印空间电位时网格点间隔数：%d；\n", &NST);//扫描输出打印空间电位时网格点间隔数
//	fscanf_s(ParameterFile, "要求扫描搜索等电位线的电位间隔值或电位值:", &dov);
//	EOE = (double*)calloc(100, sizeof(double));
//	dov = 0;
//	i = -1;
//	do
//	{
//		i++;
//		fscanf_s(ParameterFile, "%lf；", &EOE[i]);
//
//	} while (EOE[i] != 0);
//	if (i == 1)
//		dov = EOE[0];
//
//
//	//fscanf_s(ParameterFile, "要求扫描搜索等电位线的电位间隔值或电位值:%lf；\n", &dov);//扫描轴上电位作等距插值时的步长数：
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
//	//		/*电场电位初始化*/
//	//	row = M1 + M2 + 1;//Calculate the value of the row.
//	//	lin = 0;//Assign a initial value to line so that it's final value can be calculate by a circulation
//	//	for (i = 0; i < n; i++)
//	//	{
//	//		lin = lin + *(N + i);
//	//	}
//	//	lin = lin + n;
//	//	E = (double*)calloc(lin * row, sizeof(double));
//	//	s = 0;//列数上界
//	//	t = *N + 1;//列数下界
//	//	for (i = 0; i < row; i++)//行循环
//	//	{
//	//		s = 0;//列数上界
//	//		t = *N + 1;//列数下界
//	//		if (i == row - 1)//第一行输入
//	//		{
//	//			for (j = 0; j < n; j++)//电极定位
//	//			{
//	//				if (j == 0)
//	//				{
//	//					int f = 1;//定义最终循环的次数
//	//
//	//					for (k = s; k < t; k++)
//	//					{
//	//
//	//						double d = *(V + j) / *(N + j);//计算公差移到外面即可
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
//	//					int f = 1;//定义最终循环的次数
//	//
//	//					for (k = s; k < t; k++)
//	//					{
//	//
//	//						double d = (*(V + j) - *(V + j - 1)) / *(N + j);//计算公差移到外面即可
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
//	//			for (j = 0; j < n; j++)//电极定位
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
//	//				//加判断语句规避光屏
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
//	//	/*给电极位置对应格点坐标赋值*/
//	//	pon = (int*)calloc(n + 1, sizeof(int));//给电极位置对应格点坐标分配内存空间
//	//	*pon = 0;
//	//	for (i = 1; i < n + 1; i++)
//	//	{
//	//		if (i == n)
//	//		{
//	//			*(pon + i) = lin - 1;//荧光屏格点位置
//	//		}
//	//		else
//	//		{
//	//			*(pon + i) = *(pon + i - 1) + *(N + i - 1) + 1;//电极右格点位置
//	//		}
//	//	}
//	//	/*给z轴方向网络步长划分赋值*/
//	//	loz = (double*)calloc((lin - 1) * (row - 1), sizeof(double));//给z轴方向网络划分步长分配内存空间
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
//	//	/*给z轴格点位置坐标赋值*/
//	//	poz = (double*)calloc(lin, sizeof(double));//给z轴格点位置坐标分配内存空间
//	//	*poz = 0;
//	//	for (i = 1; i < lin; i++)
//	//		*(poz + i) = *(poz + i - 1) + *(loz + i - 1);
//	//	/*给r轴方向网络步长划分赋值*/
//	//	lor = (double*)calloc((lin - 1) * (row - 1), sizeof(double));//给lor分配内存空间
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
//	//	/*给r轴格点位置赋值*/
//	//	por = (double*)calloc(row, sizeof(double));//
//	//	*por = 0;
//	//	for (i = 1; i < row; i++)
//	//		*(por + i) = *(por + i - 1) + *(lor + (i - 1) * (lin - 1));
//	//
//	//	/////////////////////////////////////////////	
//	//		/*开始计算像管电位分布*/
//	//		/*第一轮迭代*/
//	//	turn = 1;
//	//	times = 1;
//	//	atime = 0;
//	//	w = 1;
//	//	SOR();
//	//	atime = atime + 1;
//	//
//	//	/*第二轮迭代*/
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
//	//	/*迭代若干轮以固定迭代因子*/
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
//	//	/*继续多轮迭代以达到精度要求*/
//	//	do
//	//	{
//	//		turn++;
//	//		SOR();   //最佳迭代因子确定后继续迭代
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
//	//	   /*扫描等位点坐标*/
//	//	double a = 0;
//	//	int m = 1;
//	//	if (dov != 0)
//	//	{
//	//		for (m = 1; a < (V[n - 1] - dov); m++)   //计算电位数量
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
//	//	//for (m = 1; a < (V[n - 1] - dov); m++)   //计算电位数量
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
//	//		fprintf(OutputFile, "第%d电位,电位值：%lf,\n", k + 1, EOE[k]);
//	//		for (j = 0; j < nopoE[k] - 1; j++)
//	//		{
//	//			fprintf(OutputFile, "(%3.2lf,%3.2lf)\t", pozE[k * 2 * lin + j], porE[k * 2 * row + j]);
//	//		}
//	//		fprintf(OutputFile, "(%3.2lf,%3.2lf)\t\n", pozE[k * 2 * lin + nopoE[k] - 1], porE[k * 2 * row + nopoE[k] - 1]);
//	//	}//输出等位线坐标
//	//
//	//	 ////////为方便制图，开始坐标变换////////////
//	//	/*由于该程序将光轴与阴极面交点设为原点，而画图软件则是将左上角定位原点，因此为了观察习惯，将r轴方向坐标进行变换*/
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
//	//		 /*开始绘制等位线*/
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
//	///*检验程序运行过程子函数*/
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
//	///*超张弛迭代法子函数*/
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
//	//				c1 = 2 / ((*(loz + i * (lin - 1) + j - 1)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));   //系数计算
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
//	//				c1 = 2 / ((*(loz + i * (lin - 1) + j - 1)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));   //系数计算
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
//	//			for (int k = 1; k < n + 1; k++)//电极循环
//	//			{
//	//				if (k != n)
//	//				{
//	//					for (int j = *(pon + k - 1) + 1; j < *(pon + k) - 1; j++)//列数循环，规避边界点
//	//					{
//	//						c1 = 2 / ((*(loz + i * (lin - 1) + j - 1)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));   //系数计算
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
//	//						c1 = 2 / ((*(loz + i * (lin - 1) + j - 1)) * (*(loz + i * (lin - 1) + j - 1) + *(loz + i * (lin - 1) + j)));   //系数计算
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
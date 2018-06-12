/****************************
*领域覆盖算法演示程序
*****************************/

#include <iostream>
#include <algorithm> //算法库
#include <fstream>   // 对文件操作的函数集
#include <string>
#include <vector>
#include <cmath>
#include <cstring> 
#include <cstdlib>
#include <cstdio> 
#include <ctime>   //把日期和时间转换成字符串
#include <cfloat>
#include<math.h>
#include<iomanip>
using namespace std;

typedef struct Sample
{	//样本数据结构
	int dim;			//样本维数
	vector<double>	x;	//样本向量
	int				y;	//样本类别
}Sample;

typedef struct SampleTag
{//样本标号
	int id; 			//样本序号
	int covered;		//这个序号的样本是否被覆盖, 1 表示被覆盖，0表示没有没覆盖 
}SampleTag;

typedef struct SampleSubSet
{	//样本子集标号，它表示同一类所有样本的标号			
	int y; 					//子集的类别 
	vector<SampleTag*> v;		//id的集合
}SampleSubSet;

typedef struct Cover
{	//覆盖的数据结构
	int *sid;			// 样本序号
	int seq;			//覆盖序号
	int cls;			//覆盖的类别
	double r;		 	//覆盖半径
	double* center;	//覆盖的圆心
	double ybs;			//覆盖的样本数 
}Cover;

typedef struct ExpResult
{	//实验结果数据结构
	int a0, a1, a2, a3, a4, a5;
	int trNum;			//训练样本数
	int teNum;			//测试样本数
	int covNum;		//覆盖数
	double trTime;		//训练耗时
	double teTime;		//测试耗时
	int refuse;		//测试拒识样本数
	int guess_corr;	//拒识样本中根据圆心最近猜测正确的样本数
	int correct;		//正确识别样本数
	double corr_rate;	//正确识别率，(corr_rate + guess_corr)/teNum
	int C0N, C1N;//C0,C1类覆盖的数目
}ExpResult;



//########################################################################## 
//全局变量 
int G_lev = 0; //保存最优层组合
int Cn = 0;
vector<int> BND_sample_id;      //保存 边界域 样本 id
vector<int> POS_sample_id;      //保存 正域   样本 id
vector<int> NEG_sample_id;      //保存 负域   样本 id
vector<Sample*> samples;		//样本集合
vector<SampleSubSet*> I;		//样本子集编号
vector<Cover*> C;				//覆盖
vector<Cover*> C1;
vector<Cover*> C2, C3, C4;
vector<ExpResult*> Result;		//实验结果	
int SubSet_POS[4][6] = {
	
	{ 1, 2, 3, 4, 5, 6 },//0

	{ 1, 2, 4, 5, 6 },//0.0015

	//{ 1, 2, 4, 6 },//0.0025

	{ 1, 2, 6 },//0.0095

	{ 1, 6 },//0.01505
};
int SubSet_NEG[5][6] = {
	{ 1, 2, 3, 4, 5, 6},//0.0109
	{ 1, 3, 4, 5, 6 },//0.0449
	{ 1, 4, 5, 6 },//0.0759
	{ 1, 4, 5 },//0.1169
	{ 1, 5 }//0.1279
};


/*********************************************************************
*	这个函数是样本读取的辅助函数，将字符串str用字符spliter分开，
*	并存入vector<string>中
**********************************************************************/
void split(string& str, char spliter, vector<string>& vec)
{
	string::iterator iter1, iter2;  //迭代器(iterator)是一种允许程序员检查容器内元素，并实现元素遍历的数据类型
	if (str == "")
		return;

	iter1 = str.begin();
	iter2 = str.begin();

	while (iter2 != str.end())
	{
		if (*iter2 != spliter)
		{
			iter2++;
			continue;
		}
		else
		{
			if (iter1 != iter2) vec.push_back(string(iter1, iter2));
			iter2++;
			iter1 = iter2;
		}

	}//while

	//处理最后一个字符串
	if (iter1 != str.end() && *iter1 != spliter)
	{
		vec.push_back(string(iter1, iter2));
	}
}

/*********************************************************************
*	这个函数是样本读取的辅助函数，将一个符合UCI格式的一行样本字符串s
*	转换成自定义的sample格式
**********************************************************************/
void dealSample_uci_format(string& s, Sample& sample)
{
	char spliter[] = { ',', ' ', '\t', ';', ':' };   //样本分量分隔符 
	int k = -1;
	for (int i = 0; i < 5; ++i)
	{//自动寻找样本分隔符
		string::size_type loc = s.find(spliter[i], 0);
		if (string::npos == loc)
		{
			continue;
		}
		else
		{
			k = i;
			break;
		}
	}
	if (k < 0)
	{
		printf("样本格式有错！\n");
		exit(0);
	}

	//分割样本
	vector<string> vec_str;
	split(s, spliter[k], vec_str);

	vector<string>::iterator it = vec_str.begin();
	vector<string>::iterator iend = vec_str.end();
	--iend;
	sample.dim = vec_str.size() - 1;
	sample.y = atoi((*iend).c_str());
	while (it != iend)
	{
		double value = atof((*it).c_str());
		sample.x.push_back(value);
		++it;
	}
}

/*********************************************************************
*	加载样本函数，传入样本文件名 sample_set_file，将样本值读取到
*	全局变量samples中
**********************************************************************/
void loadSample(const string& sample_set_file)
{
	if (sample_set_file == "")
	{
		cout << "样本文件不能为空!" << endl;
		getchar();
		exit(0);
	}
	//打开样本文件
	fstream ifs(sample_set_file.c_str(), ios::in); //.c_str()指向文件sample_set_file里面的内容 
	if (!ifs)
	{
		cout << "打开样本文件失败!" << endl;
		getchar();
		exit(0);
	}

	string line = "";
	while (getline(ifs, line)) // getline(cin,s)函数作用：将cin字符串赋值给s
	{
		Sample* sample = new Sample;
		dealSample_uci_format(line, *sample); //处理uci格式的样本
		samples.push_back(sample);
	}
	ifs.close();
}

/*********************************************************************
*	样本排序比较函数，排序函数会利用这个函数根据其目标值
*	对样本排序 ，a和b是样本的指针
**********************************************************************/

bool cmp(Sample* a, Sample* b)
{
	if (a->y > b->y)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*********************************************************************
*	样本排序函数，根据样本的目标值将样本排序，并将排序结果存入
*	全局变量I中，I[k]表示第k类的所有样本序号的集合
**********************************************************************/
void sortSample()
{
	//给样本集合排序，将相同的类别放在一起
	sort(samples.begin(), samples.end(), cmp); //这里已经对样本的类别进行了排序
	SampleSubSet* subset = new SampleSubSet;   //建立一个新的subset子集
	int y = samples[0]->y;
	subset->y = y;                             //将第一个样本的类别 作为子集subset的类别
	SampleTag* tag = new SampleTag;    //建立一个新的样本标号tag,并且给其赋予一个新的id 以及covered
	tag->id = 0;
	tag->covered = 0;
	subset->v.push_back(tag);          //将新建的样本标号tag添加到子集subset中，作为开始
	for (size_t i = 1; i < samples.size(); ++i) // 现在从samples第一个样本开始，以是否被覆盖为评判标准，进行分类
	{
		if (samples[i]->y == y)
		{
			tag = new SampleTag;
			tag->id = i;
			tag->covered = 0;
			subset->v.push_back(tag);
		}
		else
		{
			I.push_back(subset);     //I中存储的是子集种类，也就是有两种类别的子集，被覆盖一类 + 没有被覆盖一类
			subset = new SampleSubSet; //但是I中不止2类，应该说现在有n类，只是有很多类是相同，要是按照是否被覆盖来分类的话，总共就有2类
			y = samples[i]->y;
			subset->y = y;
			tag = new SampleTag;
			tag->id = i;
			tag->covered = 0;
			subset->v.push_back(tag);
		}
	}
	I.push_back(subset);
}

/*********************************************************************
*	求一个样本向量的模长    sqrt()函数是平方根函数
**********************************************************************/
double vlen(Sample& s)
{//求一个向量模长
	double a = 0.0;
	int i;
	for (i = 0; i < s.dim; ++i)
	{
		a += (s.x[i] * s.x[i]);
	}
	return sqrt(a);
}

/**************************************************************************
*	属性的归一化
*	输入 vector<Sample*>& samples
*	输出 归一化的samples
***************************************************************************/
void vec_normalization()
{
	vector<double> vec_max(0);
	vector<double> vec_min(0);
	double max_val;
	double min_val;
	size_t i;
	int j;
	for (j = 0; j < samples[0]->dim; ++j)
	{// 找样本中每维的最大值和最小值
		if (samples[0]->x[j] < samples[1]->x[j])
		{
			min_val = samples[0]->x[j];
			max_val = samples[1]->x[j];
		}
		else
		{
			min_val = samples[1]->x[j];
			max_val = samples[0]->x[j];
		}

		for (i = 2; i < 2 * (samples.size() / 2); i += 2)   //此处用  2*（a/2)的作用是，消除（i+2）会出现越界的情况
		{
			if (samples[i]->x[j] < samples[i + 1]->x[j])
			{
				if (min_val > samples[i]->x[j])
				{
					min_val = samples[i]->x[j];
				}
				if (max_val < samples[i + 1]->x[j])
				{
					max_val = samples[i + 1]->x[j];
				}
			}
			else
			{
				if (min_val > samples[i + 1]->x[j])
				{
					min_val = samples[i + 1]->x[j];
				}
				if (max_val < samples[i]->x[j])
				{
					max_val = samples[i]->x[j];
				}
			}
		}//for
		if (samples.size() % 2 != 0)   //检查样本个数是否为奇数，要是为奇数的话，还得和最后一个没有处理的样本进行比较
		{
			if (min_val > samples[samples.size() - 1]->x[j])
			{
				min_val = samples[samples.size() - 1]->x[j];
			}
			if (max_val < samples[samples.size() - 1]->x[j])
			{
				max_val = samples[samples.size() - 1]->x[j];
			}
		}
		vec_max.push_back(max_val);
		vec_min.push_back(min_val);
	}

	//归一化
	for (j = 0; j < samples[0]->dim; ++j)
	{
		for (size_t i = 0; i < samples.size(); ++i)
		{
			if ((vec_max[j] - vec_min[j]) < 1e-6)
				samples[i]->x[j] = 0;
			else
				samples[i]->x[j] = (samples[i]->x[j] - vec_min[j]) / (vec_max[j] - vec_min[j]);
		}
	}
}

/*********************************************************************
*	将样本投射到球面 ，并进行归一化处理
**********************************************************************/
void projectToSphere()
{
	size_t i;
	double d = vlen(*samples[0]);
	for (i = 1; i < samples.size(); ++i)
	{//找最大模长的向量长度d
		double t = vlen(*samples[i]);
		if (d < t)
			d = t;
	}
	for (i = 0; i < samples.size(); ++i)
	{//增加一维，映射到球面
		double x = vlen(*samples[i]);
		double t = d*d - x*x;
		t = sqrt(t);
		samples[i]->x.push_back(t); //增加一维
		samples[i]->dim += 1;
	}

	for (i = 0; i < samples.size(); ++i)
	{//归一化处理
		double x = vlen(*samples[i]);
		for (int k = 0; k < samples[i]->dim; ++k)
		{
			samples[i]->x[k] /= x;
		}
	}
}

/*********************************************************************
*	求两个样本之间的欧氏距离
**********************************************************************/
double inner_product(Sample& s1, Sample& s2)
{//求欧氏距离 
	if (s1.dim != s2.dim)
	{
		cout << "s1,s2维度不相等！" << endl;
		return 0.0;
	}
	double a = 0.0;
	int i;
	for (i = 0; i < s1.dim; ++i)
	{

		a += s1.x[i] * s1.x[i] + s2.x[i] * s2.x[i] - 2 * s1.x[i] * s2.x[i];
	}
	return  sqrt(a);
}
/*********************************************************************
	求两个样本之间的欧氏距离
    在划分好的属性层上进行计算 样本之间的欧式距离
**********************************************************************/
double inner_product_subset(Sample& s1, double* x, int y,int c)
{   //在属性子集的范围内
	//计算两个样本之间的欧式距离
	
	double a = 0.0;
	if (c == 0)
	{
		int t = (sizeof(SubSet_POS[y]) / sizeof(int));  // 计算SutSet数组有多少列
		for (int j = 0; j < s1.dim; j++)
		{
			int flag = 0;
			for (int i = 0; i < t; i++)
			{
				if (j == SubSet_POS[y][i])
					flag = 1;
			}
			if (flag == 1)
				a += s1.x[j] * s1.x[j] + x[j] * x[j] - 2 * s1.x[j] * x[j];
		}
		
	}
	if (c == 1)
	{
		int t = (sizeof(SubSet_NEG[y]) / sizeof(int));  // 计算SutSet数组有多少列
		for (int j = 0; j < s1.dim; j++)
		{
			int flag = 0;
			for (int i = 0; i < t; i++)
			{
				if (j == SubSet_NEG[y][i])
					flag = 1;
			}
			if (flag == 1)
				a += s1.x[j] * s1.x[j] + x[j] * x[j] - 2 * s1.x[j] * x[j];
		}
		
	}
	return sqrt(a);
}
/*********************************************************************
*		求样本s1和数组x的欧氏距离
**********************************************************************/
double inner_product_x(Sample& s1, double* x)
{//求欧氏距离
	
	double a = 0.0;
	for (int i = 0; i < s1.dim; ++i)
	{
		//a += (s1.x[i] - x[i])*(s1.x[i] - x[i]);
		a += s1.x[i] * s1.x[i] + x[i] * x[i] - 2 * s1.x[i] * x[i];
		//cout << x[i] << "," << s1.x[i]<<"    ";

	}
	//cout << a<<"****************888888888888888888888888888888888888888888" << endl;
	return sqrt(a);
}

/*********************************************************************
*		求向量x1和向量x2的欧氏距离
**********************************************************************/
double inner_product(double* x1, double* x2, int dim)
{//求欧氏距离 
	double a = 0.0;
	for (int i = 0; i < dim; ++i)
	{
		a += (x1[i] - x2[i])*(x1[i] - x2[i]);
	}
	return sqrt(a);
}
/*********************************************************************
*	求第t类样本中s点距异类点最近值（最小欧氏距离），s是所求点所在子集id
**********************************************************************/
double find_diffMin_d(int s, int t)                                                                                                    
{
	int t1 = (t + 1) % I.size();  //t1表示（t+1)类样本子集（SampleSubSet)
	int k = I[t1]->v[0]->id;  //k表示的是：（t+1)类->样本->序号   此处样本是从0开始的
	int a = I[t]->v[s]->id;   //a表示的是 t类->s样本->序号

	double diffMin_d = inner_product(*samples[a], *samples[k]);  //将t类中s样本与其他类中的样本（随机选取）计算了二者的欧氏距离，作为最初的diffMin_d值
	for (size_t i = 0; i < I.size(); ++i)
	{
		if ((int)i != t)
		{
			for (size_t j = 0; j < I[i]->v.size(); ++j)
			{
				k = I[i]->v[j]->id; //实际id
				double x = inner_product(*samples[a], *samples[k]);
				if (x < diffMin_d)
				{
					diffMin_d = x;
				}
			}//for
		}//else
	}//for
	return diffMin_d;
}

/*********************************************************************
*	求第t类样本中s点距同类最远距离（最大欧氏距离），s是所求点所在子集id
**********************************************************************/
double find_SameMax_d(int s, int t, double diffMin_d)   //求的这个同类中最大距离，一定要比异类中的最小距离小才符合
{
	int a = I[t]->v[s]->id; //原点实际id号
	double sameMax_d = 0;//和运用内积计算之间的差别在此 
	for (size_t j = 0; j < I[t]->v.size(); ++j)
	{
		int k = I[t]->v[j]->id; //实际id号
		double x = inner_product(*samples[a], *samples[k]);
		if (x < diffMin_d)
		{
			if (x>sameMax_d)
			{
				sameMax_d = x;
			}
		}
		

	}
	return sameMax_d;
}

/*********************************************************************
*	根据中心点s，样本类别号t，覆盖半径d，做覆盖，并返回覆盖的点数，
*	函数执行过程中修改样本子集中样本号的covered属性，以标明被覆盖
**********************************************************************/
int cover_sample(int s, int t, double d, Cover *c)
{
	int cov_num = 0;
	int a = I[t]->v[s]->id;
	c->sid = new int[I[t]->v.size()];
	//int j = 0;
	for (size_t i = 0; i < I[t]->v.size(); ++i)
	{
		int k = I[t]->v[i]->id; //实际id
		if (k == a)
		{//圆心本身
			I[t]->v[i]->covered = 1;
			++cov_num;
			c->sid[cov_num - 1] = k;
			continue;
		}
		if (inner_product(*samples[a], *samples[k]) <= d&&I[t]->v[i]->covered == 0)
		{
			I[t]->v[i]->covered = 1;
			++cov_num;
			c->sid[cov_num - 1] = k;
		}
	}
	return cov_num; //返回覆盖的样本数
}

/*********************************************************************
*	训练完成后将覆盖的结果保存成文件，文件名为 model_file
**********************************************************************/
void save_model_to_file(const string& model_file)
{
	if (model_file == "")
	{
		cout << "模型文件不能为空!" << endl;
		getchar();
		exit(0);
	}
	ofstream ofs(model_file.c_str(), ios::out);
	if (!ofs == NULL)
	{
		cout << "创建模型文件失败!" << endl;
		getchar();
		exit(0);
	}
	int dim = samples[0]->dim;
	//开始写文件
	ofs << C.size() << endl;					//写入覆盖数
	ofs << dim << endl;						//写入覆盖维度，即样本投射后的维度
	ofs.precision(10);						//经测试，精度太小，保存的模型会有误差，同样本再测试时会有拒识现象
	for (size_t i = 0; i < C.size(); ++i)
	{
		ofs << C[i]->seq << ","			//写入覆盖序号
			<< C[i]->cls << ","			//写入覆盖类别
			<< C[i]->r << ",";			//写入覆盖半径
		ofs << "\t\t\t\t";
		for (int k = 0; k < dim; ++k)		//循环写入覆盖圆心
		{
			ofs << fixed << C[i]->center[k] << ", ";
		}
		ofs << endl;
	}
	ofs.close();
}

/*********************************************************************
*	从文件中读取训练过的覆盖数据， 文件名为 model_file
**********************************************************************/
void load_model_from_file(const string& model_file)
{
	if (model_file == "")
	{
		cout << "模型文件不能为空!" << endl;
		getchar();
		exit(0);
	}
	ifstream ifs(model_file.c_str(), ios::in);
	if (!ifs == NULL)
	{
		cout << "打开模型文件失败!" << endl;
		getchar();
		exit(0);
	}
	int covnum = 0;
	int dim = 0;
	string line;
	//读取覆盖数据
	ifs >> covnum;	//读取覆盖数目
	ifs >> dim;		//读取覆盖维度
	getline(ifs, line);
	for (int i = 0; i < covnum; ++i)
	{
		Cover* c = new Cover();
		c->center = new double[dim];
		vector<string> vec_str;

		getline(ifs, line);
		split(line, ',', vec_str);

		c->seq = atoi(vec_str[0].c_str());	//覆盖序号
		c->cls = atoi(vec_str[1].c_str());	//覆盖类别
		c->r = atof(vec_str[2].c_str());	//覆盖半径

		for (int k = 0; k < dim; ++k)
		{
			c->center[k] = atof(vec_str[k + 3].c_str());
		}
		C.push_back(c);
	}
	ifs.close();
}

//根据覆盖中的样本数对覆盖进行排序 
void sortCover(Cover * cover, int N)  // 冒泡排序算法
{
	int x, y;
	Cover temp;
	for (y = 0; y < N - 1; y++)
	{
		for (x = 1; x<N - y; x++)  //这里不是应该对x进行++么？为啥变成y++了
		{
			if (cover[x].ybs>cover[x - 1].ybs)
			{
				temp = cover[x - 1];
				cover[x - 1] = cover[x];
				cover[x] = temp;
			}
		}
	}
}


/*********************************************************************
*	样本训练函数
**********************************************************************/
void sample_train(double a, double b)
{//求覆盖
	clock_t t1;

	t1 = clock();
	sortSample();	//给样本按类别排序
	int seq = 0;
	int coved_num = 0;
	//srand((unsigned)time(NULL));
	for (size_t t = 0; t < I.size(); ++t)
	{
		vector<int> v;
		for (size_t i = 0; i < I[t]->v.size(); ++i)
			v.push_back(i);  //在一个临时的集合v上操作
		int uncovedn = I[t]->v.size();
		while (!v.empty())
		{
			random_shuffle(v.begin(), v.end());
			int s = v[0];
			if (I[t]->v[s]->covered == 1)
			{
				v.erase(v.begin());
				continue;
			}
			double d1 = find_diffMin_d(s, t); //找异类最近点 内积最大,欧式距离最小 d1 ;
			double d2 = find_SameMax_d(s, t, d1); //找同类最远点 内积最小，欧式距离最大 d2
			//double d = (d1+d2)/2 ;
			double d = d2;


			Cover* c = new Cover;

			//根据半径对样本做覆盖标记
			coved_num = cover_sample(s, t, d, c); //s是圆心，t是第t类样本，d是覆盖半径

			//构造一个覆盖

			int k = I[t]->v[s]->id;
			int dim = samples[k]->dim;
			//c->sid = k;
			c->seq = seq;
			c->center = new double[dim];
			for (int i = 0; i < dim; ++i)
			{
				c->center[i] = samples[k]->x[i];
			}
			c->cls = samples[I[t]->v[0]->id]->y;
			c->r = d;
			c->ybs = coved_num;

			if (c->cls == 0)
				C1.push_back(c);//正类	
			else
				C2.push_back(c);//负类 


			//修改未覆盖数
			uncovedn -= coved_num;
			//printf("%d\t第%d类\t本次子覆盖:%d\t覆盖半径：%f\n",seq,t,coved_num,d);
			seq++;
			v.erase(v.begin());
		}


	}//for

	//对C1,C2进行排序 
	//sortCover(C1[0],C1.size());
	//sortCover(C2[0],C2.size());

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//将C1类和C2类进行处理，根据a,b值进行释放。（a,b） ;//重新形成覆盖集C3,C4 
	//int i;
	//for(i=0;i<C1.size()*(1-a);i++)
	//C3.push_back(C1[i]);
	//for(i=0;i<C2.size()*(1-b);i++)
	//C4.push_back(C2[i]);
	//
	////存入C中 
	//for(i=0;i<C3.size();i++)
	//C.push_back(C3[i]);
	//for(i=0;i<C4.size();i++)
	//C.push_back(C4[i]); 
	//cout<<"C.size="<<C.size()<<".";
	//
	//	t2 = clock() ;
	//cout<<" T="<<(t2-t1)/CLOCKS_PER_SEC;
	//	//存储实验数据
	//	vector<ExpResult*>::reverse_iterator it = Result.rbegin() ;
	//	(*it)->covNum = C.size() ; //覆盖数
	//	(*it)->trNum = samples.size() ; //训练样本数
	//	(*it)->trTime = (double)(t2-t1)/CLOCKS_PER_SEC ;
	//	(*it)->C0N= C3.size();
	//	(*it)->C1N= C4.size();
	//
	//cout<<"C3="<<C3.size()<<" ";
	//cout<<"C4="<<C4.size()<<" ";
	//C1.clear();
	//C2.clear();
	//C3.clear();
	//C4.clear();
}



/*********************************************************************
*	样本测试函数
**********************************************************************/
void sample_test()
{//测试
	//int a[6] = { 0 };
	int a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0;
	int correct = 0;		//可识别且正确分类数
	int refuse = 0;		//拒识样本数
	int uc = 0;			//样本拒识时根据离覆盖最近 来划分样本类别的情况下 能正确分类数
	int total_BND = 0;
	int total_correct_BND = 0;
	double correct_index = 0.0;
	//clock_t t1, t2;
	//t1 = clock();

	for (size_t i = 0; i < BND_sample_id.size(); ++i)
	{//对于边界域样本测试，能否正确测试

		double cnt_nearest = DBL_MAX; //中心最近点
		//double cnt_nearest1=DBL_MAX;
		int k = -1; //记录中心最近点的覆盖的下标
		//int k1= -1;
		size_t j = 0;

		while (j<C1.size())
		{//每个覆盖
			
			double r = C1[j]->r;

			if (r <= 0.000001) 
			{
				++j; continue;
			}
			else
			{
			
			//cout <<j<<","<< C1[j]->center[0] << "***************************" << endl;
			 double d = inner_product_x(*samples[BND_sample_id[i]], C1[j]->center);

         
			 if (cnt_nearest>(d - r))
			{
				cnt_nearest = d - r;
				k = j;
			}
			/*	if(cnt_nearest1>d)
			{
			cnt_nearest1 = d;
			k1 = j 
			}*/
			++j;
			}
		}//while



		if (cnt_nearest > 0)
		{//此点距最近的的覆盖都在半径之外，则为拒识样本点
			++refuse;
			if (samples[BND_sample_id[i]]->y == C1[k]->cls) //根据距离覆盖中心最近来给拒识样本划分
				++uc;
			if (samples[BND_sample_id[i]]->y == 0)
				++a4;//++a[4];           //a4记录没有被覆盖的样本数
			if (samples[BND_sample_id[i]]->y == 1)
				++a5;//++a[5];		     //a5记录被覆盖的样本数

		}
		else
		{//可识样本
			if (samples[BND_sample_id[i]]->y == C1[k]->cls)
			{

				++correct;
				if (samples[BND_sample_id[i]]->y == 0)//正类 
					++a0;//++a[0];
				if (samples[BND_sample_id[i]]->y == 1)
					++a1;//++a[1];

			}
			else
			{//  边界域
				if (samples[BND_sample_id[i]]->y == 0)
					++a2;//++a[2];
				if (samples[BND_sample_id[i]]->y == 1)
					++a3;//++a[3];		
			}
		}


	}//for

	for (int i = 0; i < BND_sample_id.size(); i++)
	{
		total_BND++;
		if (samples[BND_sample_id[i]]->y == 0)
			total_correct_BND++;

	}
	correct_index = a0 / total_correct_BND;
	cout << "边界域中样本总数：" << total_BND << endl;
	cout << "边界域中样本被分到正域中，且被正确划分的数量为：" << a0 << endl;
	cout << "正确率为:" << correct_index << endl;


	

	//t2 = clock();
	//cout << " l to l:" << a0 << "," << "l to s:" << a2 << "," << "l to B:" << a4 << "," << "s to l:" << a3 << "," << "s to s:" << a1 << "," << "s to B:" << a5 << endl;

	////存储实验数据
	//vector<ExpResult*>::reverse_iterator it = Result.rbegin();
	//(*it)->a0 = a0;
	//(*it)->a1 = a1;
	//(*it)->a2 = a2;
	//(*it)->a3 = a3;
	//(*it)->a4 = a4;
	//(*it)->a5 = a5;
	//(*it)->correct = correct; //测试正确数
	//(*it)->refuse = refuse; //训练样本数
	//(*it)->guess_corr = uc;
	//(*it)->teNum = samples.size();
	//(*it)->teTime = (double)(t2 - t1) / CLOCKS_PER_SEC;
	//(*it)->corr_rate = ((float)(correct + uc)) / samples.size();
}




/************************************************************************/
/* 输出正类，负类，边界域数据                                              */
/************************************************************************/
void output_three_way_data(int bint)
{
	
	
	ofstream ofBND;
	ofBND.open("Car边界域.txt");

	cout << "正在生成正类数据..." << endl;
	ofstream ofPOS;
	ofPOS.open("Car正域.txt");
	/*int sum = 0;
	for (vector<Cover *>::iterator iter = C1.begin(); iter != (C1.end()); iter++)
	{
		
		sum = sum + (*iter)->ybs;
		
				
	}
	cout << "zhenglei=" << sum;
	for (vector<Cover *>::iterator iter = C2.begin(); iter != (C2.end()); iter++)
	{

		sum = sum + (*iter)->ybs;

	}
	cout << "zhenglei+fulei=" << sum;*/

	for (vector<Cover *>::iterator iter = C1.begin(); iter != (C1.end()); iter++)
	{
		for (int i = 0; i < (*iter)->ybs; i++)
		{
			int sid = (*iter)->sid[i];

				// 边界类
				if ((*iter)->ybs < bint)
				{
					//将边界域中的样本序号保存在 BND_sample_id数组中
					BND_sample_id.push_back(sid);

					// 获取到单个样本中归一化的数字
					vector<double>	x1 = samples[sid]->x;
					//ofBND << sid << ",";        //这里是将样本序号添加在 第一位了，导致输出的样本变成了 37 维
					for (int i = 0; i < (x1.size() - 1); i++)
					{
						ofBND << x1[i] << ",";
					}
					ofBND << (*iter)->cls;
					ofBND << endl;
					//(*iter)->ybs = 0;
					//(*iter)->center = 0;
				}
				// 正类
				else {
					// 获取到单个样本中归一化的数字
					POS_sample_id.push_back(sid); //将分到正域的样本id保存在 POS_sample_id 中
					//vector<double>	x2 = samples[sid]->x;
					vector<double>	x2 = samples[sid]->x;
					for (int i = 0; i < (x2.size() - 1); i++)
					{
						ofPOS << x2[i] << ",";
					}
					ofPOS << (*iter)->cls;
					ofPOS << endl;
				}
			
		}

		if ((*iter)->ybs < bint)
		{
			(*iter)->ybs = 0;
			(*iter)->r = 0;
		}
		
	}
	ofPOS.close();

	cout << "正在生成负类数据..." << endl;
	ofstream ofNEG;
	ofNEG.open("Car负域.txt");
	for (vector<Cover *>::iterator iter = C2.begin(); iter != (C2.end()); iter++)
	{
		for (int i = 0; i < (*iter)->ybs; i++)
		{
			int sid = (*iter)->sid[i];
			
				// 边界类
				if ((*iter)->ybs < bint)
				{
					//将边界域中的样本序号保存在 BND_sample_id数组中
					BND_sample_id.push_back(sid);
					// 获取到单个样本中归一化的数字
					vector<double>	x1 = samples[sid]->x;
					//ofBND << sid << ",";
					for (int i = 0; i < x1.size() - 1; i++){
						ofBND << x1[i] << ",";
					}
					ofBND << (*iter)->cls;
					ofBND << endl;

					//(*iter)->ybs = 0;
					//(*iter)->center = 0;
				}
				// 负类
				else {
					// 获取到单个样本中归一化的数字
					NEG_sample_id.push_back(sid);   //将负域中的样本id 保存在  NEG_sample_id 中
					vector<double>	x2 = samples[sid]->x;
					for (int i = 0; i < x2.size() - 1; i++)
					{
						ofNEG << x2[i] << ",";
					}
					ofNEG << (*iter)->cls;
					ofNEG << endl;
				}
			
			
			
		}
		if ((*iter)->ybs < bint)
		{
			(*iter)->ybs = 0;
			(*iter)->r = 0;
		}
	}
	ofNEG.close();

	cout << "正在生成边界类数据..." << endl;

	//for (size_t i = 0; i < I.size() - 1; ++i){ 

	//	vector<SampleTag*> st = I[i]->v;

	//	for (vector<SampleTag*>::iterator iter = st.begin(); iter != (st.end() - 1); iter++)
	//	{
	//		// 边界域数据
	//		if ((*iter)->covered == 0)
	//		{
	//			vector<double>	xx = samples[(*iter)->id]->x;
	//			for (int i = 0; i < xx.size(); i++){
	//				ofBND << xx[i] << ",";
	//			}
	//			ofBND << I[i]->y;
	//			ofBND << endl;
	//		}
	//	}
	//}


	ofBND.close();
	

}

/*
 边界域样本 划分类别函数
 在 BND_test()中直接调用
*/
int sample_test_all(Sample& s)
{
	size_t j1 = 0;
	size_t j2 = 0;
	double d1 = 0.0;  //记录边界域样本与覆盖样本之间的距离，并且保存最小距离
	double d2 = 0.0;
	int m1 = -1;      //记录中心最近点的覆盖下标
	int m2 = -1;
	int y = -1;        //样本被划分后的类别



	double cnt1 = DBL_MAX;
	//double cnt1_nearest = DBL_MAX;
	double cnt2 = DBL_MAX;
	//cout << "************2**********" << endl;
	while (j1 < C1.size())
	{
		//cout << "************3**********" << endl;
		double r = C1[j1]->r;
		if (r <= 0.000001)
		{
			++j1;
			continue;
		}
		else
		{
			//cout << "************4**********" << endl;
			double d = inner_product_x(s, C1[j1]->center);
			//cout << "************2**********" << endl;
			if (cnt1 > (d-r))
			{
				//cout << "************11**********" << endl;
				cnt1 = d-r;
				m1 = j1;
				d1 = d;
			}
			++j1;
		}
	}
	while (j2 < C2.size())
	{
		double r = C2[j2]->r;
		//cout << "************5**********" << endl;
		if (r <= 0.000001)
		{
			++j2;
			continue;
		}
		else
		{
			double d = inner_product_x(s, C2[j2]->center);
			//cout << "************6**********" << endl;
			if (cnt2 > (d-r))
			{
				cnt2 = d-r;
				m2 = j2;
				d2 = d;
			}
			++j2;
		}
	}
	//cout <<endl<< "cnt1=" << cnt1 << ",cnt2=" << cnt2 << endl;
	//int x1 = m1;
	//int x2 = m2;
	if (cnt1< cnt2)
	{
		y = 0;
		return y;
	}
	else
	{
		y = 1;
		return y;
	}
}

int sample_test_sub(Sample& s,int x,int c)   //x表示数组的第几层   c表示属于正域数组层，还是负域数组层
{//x表示Sub+(Sub-)的第几层，c表示：属于Sub+  or 属于Sub-
	size_t j1 = 0;
	size_t j2 = 0;
	double d1 = 0.0;  //记录边界域样本与覆盖样本之间的距离，并且保存最小距离
	double d2 = 0.0;
	//int m1 = -1;      //记录中心最近点的覆盖下标
	//int m2 = -1;
	int y = -1;        //样本被划分后的类别

	

	double cnt1 = DBL_MAX;
	//double cnt1_nearest = DBL_MAX;
	double cnt2 = DBL_MAX;
	//cout << "************2**********" << endl;
	while (j1 < C1.size())
	{
		//cout << "************3**********" << endl;
		double r = C1[j1]->r;
		if (r <= 0.000001)
		{
			++j1;
			continue;
		}
		else
		{
			//cout << "************4**********" << endl;
			double d = inner_product_subset(s, C1[j1]->center,x,c);
			//cout << "************2**********" << endl;
			if (cnt1 > d)
			{
				//cout << "************11**********" << endl;
				cnt1 = d;
				//m1 = j1;
				//d1 = d;
			}
			++j1;
		}
	}
	while (j2 < C2.size())
	{
		double r = C2[j2]->r;
		//cout << "************5**********" << endl;
		if (r <= 0.000001)
		{
			++j2;
			continue;
		}
		else
		{
			double d = inner_product_subset(s,C2[j2]->center, x,c);
			
			if (cnt2 > d)
			{
				cnt2 = d;
				//m2 = j2;
				d2 = d;
			}
			++j2;
		}
	}
	//cout <<endl<< "cnt1=d1-r1=" << cnt1 << "," << "cnt2=d2-r2" << cnt2;

	if (cnt1< cnt2)
	{
		y = 0;
		return y;
	}
	else
	{
		y = 1;
		return y;
	}
}

void BND_test()
{/*将最优层赋值给 全局变量G_lev*/

	/*
	对每一个边界域样本进行 全属性分类+层次属性分类
	*/
	//output_three_way_data(bint);
	int a[1000] = { 0 };
	for (size_t i = 0; i < BND_sample_id.size(); i++)
	{
		/*  在全属性条件下，判断边界域样本属于哪一类

		y = 1和0的数量差不多，与y1,y2相差比较明显
		*/
		int y = sample_test_all(*samples[BND_sample_id[i]]);
		//cout << "样本序号 i= " << i << endl;
		//cout << "--全： " << y << endl;
		//int a1 = 0;

		for (size_t j = 0; j < (sizeof(SubSet_POS) / sizeof(int) / (sizeof(SubSet_POS[0]) / sizeof(int))); j++)
		{
			/*在正域划分的属性条件下，判断边界域样本属于哪一类*/

			/*for (int i = 0; i < sizeof(SubSet_POS[j]); i++)
			cout << SubSet_POS[j][i]<<",";
			cout << endl << "***********" << endl;*/

			int y1 = sample_test_sub(*samples[BND_sample_id[i]], j, 0);  //j表示数组的第几层   0表示在正域的数组层里

			//cout << "--正： " << y1 << ",";
			//cout << "******" << y1 << "_________" << endl;
			/*
			y1 基本上都为1，也就是正域选择的层，大都将边界域样本划分为 负类
			*/
			for (size_t k = 0; k < (sizeof(SubSet_NEG) / sizeof(int) / (sizeof(SubSet_NEG[0]) / sizeof(int))); k++)
			{
				/*在负域划分的属性条件下，判断边界域样本属于哪一类

				y2基本上都为1，也就是负域划分的层，大都将边界域样本划分为 负类
				*/
				int y2 = sample_test_sub(*samples[BND_sample_id[i]], k, 1);
				//cout << "--负：" << y2 << ",";
				if (y1 = y2 && samples[BND_sample_id[i]]->y == y2 || y == y1&&samples[BND_sample_id[i]]->y == y || y == y2&&samples[BND_sample_id[i]]->y == y)
				{
					int a1 = j*((sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int))) + k;
					a[a1] = a[a1] + 1;
				}
			}
			//cout << endl<< "--------*---------" << endl;

		}

		//cout << endl<<"_______#______"<<endl;
	}

	int len1 = (sizeof(SubSet_POS) / sizeof(int)) / (sizeof(SubSet_POS[0]) / sizeof(int));
	int len2 = (sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int));
	int len = len1*len2;

	int max = DBL_MIN;
	int i_max = 0;
	cout << endl;
	for (int i = 0; i < len; i++)
	{
		cout << "---------------" << i << "--------------" << endl;
		cout << a[i] << endl;
		if (a[i]>max)
		{
			max = a[i];
			i_max = i;
		}
	}

	G_lev = i_max;  //G_lev保存最优层组合
	int count_BND = BND_sample_id.size();
	double rate = 0.0;
	rate = (double)max / (double)count_BND;
	int lev_1 = G_lev / ((sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int)));
	int lev_2 = G_lev % ((sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int)));

	cout << endl << "**************************最优层组合结果如下***********************" << endl << endl;
	cout << "最优层组合序号 = " << i_max << ",Sub+ = " << lev_1 << ",Sub- = " << lev_2 << endl << endl;
	cout << "最优层组合 正确分类个数 = " << max << endl << endl;

	cout << "********************* BND_test() 运行结束****************************" << endl << endl;
}

void samples_test()
{

	int nc1 = 0;  //计算通过 覆盖算法 将测试样本划分为 正域且正确 的样本数量
	int nc2 = 0;  //计算通过 覆盖算法 将测试样本划分为 负域且正确 的样本数量

	int nc1_n = 0;//通过覆盖算法，将样本划分为 正域，但是与样本本身的类别不同，即划分错误的个数
	int nc2_n = 0;//通过覆盖算法，将样本划分为 负域，但是与样本本身的类别不同，即划分错误的个数



	int err = 0;

	int bnd = 0;//保存边界域样本总数；
	int bnd_1 = 0;//保存边界域样本 正域个数
	int bnd_2 = 0;//保存边界域样本 负域个数

	int all = 0;
	int all_1 = 0;
	int all_2 = 0;

	int sub_10 = 0;
	int sub_11 = 0;
	int sub_12 = 0;

	int sub_20 = 0; //保存最优Sub-测试边界域样本，正确分类的个数
	int sub_21 = 0;//保存 正类个数
	int sub_22 = 0;

	int zuhe = 0;//保存最优属性组合 测试边界域样本，正确分类个数
	int zuhe_1 = 0;
	int zuhe_2 = 0;

	int b_concret = 0;//利用求最优层的组合方式，对边界域样本用最优层组合 进行分类
	int b_concret_1 = 0;
	int b_concret_2 = 0;

	int great_lev = G_lev;  //用最优层，对边界域样本进行分类
	for (int i = 0; i < samples.size(); i++)
	{
		int flag1 = 0; // 用来标志 样本是否被 覆盖
		//int flag2 = 0;

		int j1 = 0; //对C1进行循环操作
		int j2 = 0; //对C2进行循环操作
		double d1 = DBL_MAX; //保存 样本与 C1覆盖中心 的 最小距离距离差
		double d2 = DBL_MAX; //保存 样本与 C2覆盖中心 的 最小距离距离差

		//int x1 = 0;
		//int x2 = 0;
		while (j1 < C1.size())
		{ // 计算样本在 C1 覆盖中的最小距离
			double r = C1[j1]->r;
			if (r < 0.000001)
			{
				j1++;
				continue;
			}
			else
			{
				double d = inner_product_x(*samples[i], C1[j1]->center);
				if (d1>(d - r))
				{
					d1 = d - r;
					//x1 = j1;
				}
				j1++;
			}
		}

		while (j2 < C2.size())
		{
			double r = C2[j2]->r;
			if (r < 0.000001)
			{
				j2++;
				continue;
			}
			else
			{
				double d = inner_product_x(*samples[i], C2[j2]->center);
				if (d2>(d - r))
				{
					d2 = d - r;
					//x2 = j2;
				}
				j2++;
			}
		}
		if ((d1 - d2) < 0)
		{
			if (d1 < 0 && samples[i]->y == 0)
			{
				nc1++;
				flag1 = 1;
			}
			if (d1 < 0 && samples[i]->y == 1)
			{
				nc1_n++;
				flag1 = 1;
			}
		}
		if ((d2 - d1) < 0)
		{
			if (d2 < 0 && samples[i]->y == 1)
			{
				nc2++;
				flag1 = 1;
			}
			if (d2 < 0 && samples[i]->y == 0)
			{
				nc2_n++;
				flag1 = 1;
			}
		}


		/*
		上面 C1 , C2 循环，是将样本进行二分类
		经过上面分类过程后， 未对样本进行的样本 进一步分类，
		这里用 训练找到的最优层 对样本 进一步二分类
		*/
		if (flag1 == 0)
		{
			int lev_1 = great_lev / ((sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int)));
			int lev_2 = great_lev % ((sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int)));

			bnd = bnd + 1; //保存边界域样本个数
			if (samples[i]->y == 0)
			{
				bnd_1 = bnd_1 + 1; //边界域样本中 正域样本个数
			}
			if (samples[i]->y == 1)
			{
				bnd_2 = bnd_2 + 1; //边界域样本中 负域样本个数
			}

			int y = sample_test_all(*samples[i]); //根据全属性，对样本进行分类
			int y1 = sample_test_sub(*samples[i], lev_1, 0);//用最优Sub+属性层，对样本进行分类
			int y2 = sample_test_sub(*samples[i], lev_2, 1);//用最优Sub-属性层，对样本进行分类

			//cout << "y=" << y << ",";
			if (y == samples[i]->y) // 全属性测试
			{
				//cout << "y==samples[i]->,";
				all = all + 1;
				if (y == 0)
				{
					all_1 = all_1 + 1;
				}
				if (y == 1)
				{
					all_2 = all_2 + 1;
				}
			}

			//cout << "y1=" << y1 << ",";
			if (y1 == samples[i]->y) //最优Sub+对 边界域样本进行测试
			{
				//cout<<"y1=samples[i]->y,";
				sub_10 = sub_10 + 1;
				if (y1 == 0)
				{
					sub_11 = sub_11 + 1;
				}
				if (y1 == 1)
				{
					sub_12 = sub_12 + 1;
				}
			}

			//cout << "y2=" << y2 << ",";
			if (y2 == samples[i]->y)//最优Sub-对 边界域样本进行测试
			{
				//cout << "y2==samples[i]->y,";
				sub_20 = sub_20 + 1;
				if (y2 == 0)
				{
					sub_21 = sub_21 + 1;
				}
				if (y2 == 1)
				{
					sub_22 = sub_22 + 1;
				}
			}

			if (y1 = y2 && samples[i]->y == y2 || y == y1&&samples[i]->y == y || y == y2&&samples[i]->y == y)
			{
				b_concret = b_concret + 1;
				if (samples[i]->y == 0)
					b_concret_1 = b_concret_1 + 1;
				if (samples[i]->y == 1)
					b_concret_2 = b_concret_2 + 1;
			}

			j1 = 0;
			j2 = 0;
			int flag3 = 0;
			int o1 = -1;   //保存样本的类别
			int o2 = -1;
			int o3 = -1;

			double min11 = DBL_MAX;
			double min21 = DBL_MAX;
			double min31 = DBL_MAX;

			double min12 = DBL_MAX;
			double min22 = DBL_MAX;
			double min32 = DBL_MAX;




			//cout << endl << "great_lev=" << great_lev << "-------------" << (sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int))<<endl;


			//cout <<endl<< "lev_1=" << lev_1 << "----------   lev_2" << lev_2 << endl;
			while (j1 < C1.size())
			{
				double d11 = inner_product_subset(*samples[i], C1[j1]->center, lev_1, 0);
				if (min11>d11)
					min11 = d11;

				double d21 = inner_product_subset(*samples[i], C1[j1]->center, lev_2, 1);
				if (min21 > d21)
					min21 = d21;

				double d31 = inner_product_x(*samples[i], C1[j1]->center);
				double r = C1[j1]->r;
				if (min31 > (d31 - r))
					min31 = d31 - r;
				j1++;
			}

			while (j2 < C2.size())
			{
				double d12 = inner_product_subset(*samples[i], C2[j2]->center, lev_1, 0);
				if (min12>d12)
					min12 = d12;

				double d22 = inner_product_subset(*samples[i], C2[j2]->center, lev_2, 1);
				if (min22 > d22)
					min22 = d22;

				double d32 = inner_product_x(*samples[i], C2[j2]->center);
				double r = C2[j2]->r;
				if (min32 > (d32 - r))
					min32 = d32 - r;
				j2++;
			}
			if (min11 < min12)
				o1 = 0;
			else
				o1 = 1;

			if (min21 < min22)
				o2 = 0;
			else
				o2 = 1;

			if (min31 < min32)
				o3 = 0;
			else
				o3 = 1;

			if (o1 == o2 && samples[i]->y == o1 || o1 == o3  && samples[i]->y == o1 || o2 == o3  && samples[i]->y == o2)
			{
				zuhe = zuhe + 1;
				if (samples[i]->y == 0)
				{
					zuhe_1 = zuhe_1 + 1;
				}
				if (samples[i]->y == 1)
				{
					zuhe_2 = zuhe_2 + 1;
				}
				//b2++;
				flag3 = 1;
			}
			if (flag3 == 0)
			{
				err++;
			}



		}
		//cout << endl;
	}
	int correct = 0;
	int fault = 0;

	correct = nc1 + nc2 + zuhe;
	fault = nc1_n + nc2_n + err;

	double rate;//总体正确率
	rate = double(correct) / double(samples.size());

	double rate_all = 0;//全属性正确率
	rate_all = double(all) / double(bnd);

	double rate_sub1 = 0;// Sub+正确率
	rate_sub1 = double(sub_10) / double(bnd);

	double rate_sub2 = 0;// Sub-正确率
	rate_sub2 = double(sub_20) / double(bnd);

	double rate_zuhe = 0;// 组合正确率
	rate_zuhe = double(zuhe) / double(bnd);

	cout << endl;
	cout << "+++++++++++最优组合 正确分类个数 = " << b_concret << ",正域个数 = " << b_concret_1 << ",负域个数 = " << b_concret_2 << endl << "正确率 = " << double(b_concret) / double(bnd) << endl << endl << endl;

	cout << endl << "----------边界域 = " << bnd << ",其中 正域样本个数 = " << bnd_1 << ",负域样本个数 = " << bnd_2 << endl << endl;
	cout << "根据 全属性对边界域样本进行分类，正确分类个数 = " << all << ",其中 正类样本个数 = " << all_1 << ",负类样本个数 = " << all_2 << endl;
	cout << "-----------正确率 = " << rate_all << endl << endl;
	cout << "根据 最优Sub+属性层对边界域样本进行分类，正确分类个数 = " << sub_10 << ",其中 正类样本个数 = " << sub_11 << ",负类样本个数 = " << sub_12 << endl;
	cout << "-----------正确率 = " << rate_sub1 << endl << endl;
	cout << "根据 最优Sub-属性层对边界域样本进行分类，正确分类个数 = " << sub_20 << ",其中 正类样本个数 = " << sub_21 << ",负类样本个数 = " << sub_22 << endl;
	cout << "-----------正确率 = " << rate_sub2 << endl << endl;
	cout << "根据 最优Sub+和最优Sub-的组合 对边界域样本进行分类，正确分类个数 = " << zuhe << ",其中 正类样本个数 = " << zuhe_1 << ",负类样本个数 = " << zuhe_2 << endl;
	cout << "-----------正确率 = " << rate_zuhe << endl << endl;

	cout << endl;
	cout << "----------通过覆盖对样本进行划分，正确划分的个数 = " << nc1 + nc2 << ",其中，正类样本个数 = " << nc1 << ",负类样本个数 = " << nc2 << endl;
	cout << "----------错误划分的个数 = " << nc1_n + nc2_n << "--------------" << endl << endl;

	cout << "----------通过最优层组合，正确分类的个数 =" << zuhe;
	cout << "----------错误分类的个数 = " << err << "---------------" << endl << endl;

	cout << "----------正确分类总数为 = nc1+nc2+zuhe=" << correct << "------------" << endl << endl;
	cout << "----------错误分类总数为 =nc1_n+nc2_n+err= " << fault << "------------" << endl << endl;

	cout << "----------样本总数为 = samples.size()=" << samples.size() << "------------" << endl << endl;
	cout << "----------正确率 = " << rate << "-------------" << endl;
}


/*********************************************************************
*	对样本文件执行十交叉验证
**********************************************************************/
void nfoldCrossTest(double a, double b)
{
	const int N = 10;
	int sampleCnt = samples.size();	//样本总数
	int unit = sampleCnt / N;		//n fold 每份样本数
	int sampleSelectedCnt = N*unit;	//样本取整，为了计算公平，去掉整分后的零头

	vector<Sample*> samplesBak(samples);
	random_shuffle(samplesBak.begin(), samplesBak.end()); //打乱样本顺序
	vector<Sample*>::iterator itstart = samplesBak.begin();
	vector<Sample*>::iterator iend = samplesBak.end();

	iend -= (sampleCnt - sampleSelectedCnt);  // 去掉整分后的零头 数目
	vector<Sample*>::iterator it1 = itstart;
	vector<Sample*>::iterator it2 = itstart;
	for (int i = 1; i <= N; ++i)
	{
		cout << i << ".";
		samples.clear();
		I.clear();
		C.clear();
		it2 = it1 + unit;
		samples.insert(samples.begin(), itstart, it1);
		samples.insert(samples.begin(), it2, iend);

		ExpResult *rel = new ExpResult;
		Result.push_back(rel);	//实验数据空间准备
		sample_train(a, b);

		samples.clear();
		samples.insert(samples.begin(), it1, it2);
		
		samples_test();
		/*cout << endl;
		size_t x = samples.size;
		cout << x << endl;*/

		it1 = it2;
	}
	cout << endl;
}

/*********************************************************************
*	显示样本训练结果
**********************************************************************/
void print_train_result()
{
	printf(
		"训练结果：\n"
		"训练样本为：%d个，维度：%d，类别：%d\n"
		"得到覆盖数：%d个\n"
		"训练耗时:%.4fs\n"
		, samples.size(), samples[0]->dim - 1, I.size(), Result[0]->covNum, Result[0]->trTime
		); //samples[0]->dim - 1，因为样本在投射到球面时增加了一维
}

/*********************************************************************
*	显示样本测试结果
**********************************************************************/
void print_test_result()
{
	printf("测试结果：\n");
	printf("-----------------------------------------------------------------------");
	printf("-----------------------------------------------------------------------\n");
	printf("测试样本数\t可识样本数\t可识正确数\t可识正确率\t拒识样本数\t拒识正确数\t拒识比\t拒识正确率\t总正确率\t耗时\n");
	printf("-----------------------------------------------------------------------");
	printf("-----------------------------------------------------------------------\n");
	int t1 = 0, t2 = 0, t3 = 0; double t4 = 0.0; int t5 = 0, t6 = 0; double t7 = 0.0, t8 = 0.0, t9 = 0.0, t10 = 0.0;
	t1 = Result[0]->teNum;
	t2 = Result[0]->teNum - Result[0]->refuse;
	t3 = Result[0]->correct;
	(t2 == 0) ? (t4 = 0.0) : (t4 = (double)t3 / t2); t4 *= 100;
	t5 = Result[0]->refuse;
	t6 = Result[0]->guess_corr;
	(t1 == 0) ? (t7 = 0.0) : (t7 = (double)t5 / t1); t7 *= 100;
	(t5 == 0) ? (t8 = 0.0) : (t8 = (double)t6 / t5); t8 *= 100;
	t9 = Result[0]->corr_rate * 100;
	t10 = Result[0]->teTime;
	printf("%d\t\t%d\t\t%d\t\t%.2f%%\t\t%d\t\t%d\t\t%.2f%%\t%.2f%%\t\t%.2f%%\t\t%.4f\n",
		t1, t2, t3, t4, t5, t6, t7, t8, t9, t10);
}

/*********************************************************************
*	显示十交叉验证结果
**********************************************************************/
void print_nfold_result()
{
	printf("10交叉验证:\n");
	printf("-------------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");
	printf("训练样本数\t覆盖数\t测试样本数\t可识样本数\t可识正确数\t可识正确率\t拒识样本数\t拒识正确数\t拒识比\t拒识正确率\t总正确率\t耗时\t\t正to正\t\t正to负\t正to边界\t负to正\t负to负\t负to边界\t\tC0Num\tC1Num\n");
	printf("-------------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");

	int s0 = 0, s1 = 0, s2 = 0; double s5 = 0.0, s4 = 0.0, s6 = 0.0, s3 = 0.0, s7 = 0.0; double s8 = 0.0, s9 = 0.0, s10 = 0.0, s11 = 0.0; int s12 = 0, s13 = 0, s14 = 0, s15 = 0, s16 = 0, s17 = 0; int s21 = 0, s22 = 0;
	int t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0; double t5 = 0.0; int t6 = 0, t7 = 0; double t8 = 0.0, t9 = 0.0, t10 = 0.0, t11 = 0.0; int t12 = 0, t13 = 0, t14 = 0, t15 = 0, t16 = 0, t17 = 0; int t21 = 0, t22 = 0;
	vector<ExpResult*>::iterator it = Result.begin();
	while (it != Result.end())
	{
		t0 = 0; t1 = 0; t2 = 0; t3 = 0; t4 = 0; t5 = 0.0; t6 = 0; t7 = 0; t8 = 0.0; t9 = 0.0; t10 = 0.0; t11 = 0.0; t12 = 0; t13 = 0; t14 = 0; t15 = 0; t16 = 0; t17 = 0, t21 = 0, t22 = 0;
		t0 = (*it)->trNum;
		t1 = (*it)->covNum;
		t2 = (*it)->teNum;
		t3 = (*it)->teNum - (*it)->refuse;
		t4 = (*it)->correct;
		(t3 == 0) ? (t5 = 0.0) : (t5 = (double)t4 / t3); t5 *= 100;
		t6 = (*it)->refuse;
		t7 = (*it)->guess_corr;
		(t2 == 0) ? (t8 = 0.0) : (t8 = (double)t6 / t2); t8 *= 100;
		(t6 == 0) ? (t9 = 0.0) : (t9 = (double)t7 / t6); t9 *= 100;
		t10 = (*it)->corr_rate * 100;
		t11 = (*it)->trTime + (*it)->teTime;
		t12 = (*it)->a0;
		t13 = (*it)->a1;
		t14 = (*it)->a2;
		t15 = (*it)->a3;
		t16 = (*it)->a4;
		t17 = (*it)->a5;
		t21 = (*it)->C0N;
		t22 = (*it)->C1N;
		printf("%d\t\t%d\t%d\t\t%d\t\t%d\t\t%.2f%%\t\t%d\t\t%d\t\t%.2f%%\t%.2f%%\t\t%.2f%%\t\t%.4fs\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t%d\t%d\n",
			t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t14, t16, t15, t13, t17, t21, t22);
		s0 += t0; s1 += t1; s2 += t2; s3 += t3; s4 += t4; s5 += t5; s6 += t6; s7 += t7; s8 += t8; s9 += t9; s10 += t10; s11 += t11; s12 += t12; s13 += t13; s14 += t14; s15 += t15; s16 += t16; s17 += t17; s21 += t21; s22 += t22;
		++it;
	}
	int n = Result.size();
	s0 /= n; s1 /= n; s2 /= n; s3 /= n; s4 /= n; s5 /= n; s6 /= n; s7 /= n; s8 /= n; s9 /= n; s10 /= n; s11 /= n; s12 /= n; s13 /= n; s14 /= n; s15 /= n; s16 /= n; s17 /= n; s21 /= n; s22 /= n;
	printf("-------------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");
	printf("(平均值)\n");
	printf("%d\t\t%d\t%d\t\t%f\t\t%f\t\t%.2f%%\t\t%f\t\t%f\t\t%.2f%%\t%.2f%%\t\t%.2f%%\t\t%.4fs\t\%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
		s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s14, s16, s15, s13, s17, s21, s22);
	printf("-------------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");
}



int main()
{
	//char btnums[4];
	//cout << "请给出划分边界域的范围（格式：1或者1,2）：";

	//cin >> btnums;

	cout << "输入 1 ，表示进行训练样本，对样本进行 正域，边界域，负域 分类，并且以txt文件输出" << endl;
	cout << "输入 2 ，表示选取正域里的最优属性层，负域里的最优属性层" << endl;
	cout << "输入 3 ，表示进行 样本总体测试,以及 是否进行十折交叉测试" << endl;
	//cout << "输入 4 ，表示进行十折交叉 测试" << endl;
	cout << "请输入（1,2,3）:";
	
	int choice;
	cin >> choice;
	if (choice == 1)
	{
		int bint;
		cout << "请输入边界域的最大范围：";
		cin >> bint;
		cout << "正在加载训练样本..." << endl;
		loadSample("数据集car.txt");
		cout << "正在进行正规化处理..." << endl;
		vec_normalization();
		projectToSphere();

		cout << "正在训练样本分类..." << endl;
		sample_train(0, 0);

		cout << "输出：正域 边界域 负域 ..." << endl;
		output_three_way_data(bint);
		cout << endl << "************分类完成**********" << endl;
		
		int np = POS_sample_id.size();
		int nb = BND_sample_id.size();
		int nn = NEG_sample_id.size();
		int n = samples.size();
		
		cout << "正域样本数 = " << np << endl;
		cout << "边界域样本数 = " << nb << endl;
		cout << "负域样本数 = " << nn << endl;
		cout << "正域+负域+边界域 = " << np + nb + nn << endl;
		if ((np + nb + nn) == n)
		{
			cout << "正域+负域+边界域 = 样本总数 = " << n << endl;
		}
		else
			cout << "分类出现错误！" << endl;
		cin >> choice;
	}
	else if (choice == 2)
	{
		int bint;
		cout << "请输入边界域的最大范围：";
		cin >> bint;
		cout << "正在加载训练样本..." << endl;
		loadSample("数据集car.txt");
		cout << "正在进行正规化处理..." << endl;
		vec_normalization();
		projectToSphere();

		cout << "正在训练样本分类..." << endl;
		sample_train(0, 0);

		//cout << "输出：正域 边界域 负域 ..." << endl;
		output_three_way_data(bint);
		cout << endl << "************分类完成**********" << endl;

		int np = POS_sample_id.size();
		int nb = BND_sample_id.size();
		int nn = NEG_sample_id.size();
		int n = samples.size();

		cout << "正域样本数 = " << np << endl;
		cout << "边界域样本数 = " << nb << endl;
		cout << "负域样本数 = " << nn << endl;
		cout << "正域+负域+边界域 = " << np + nb + nn << endl;
		if ((np + nb + nn) == n)
		{
			cout << "正域+负域+边界域 = 样本总数 = " << n << endl;
		}
		else
			cout << "分类出现错误！" << endl;

		cout << "----------输出最优层----------"<<endl;
		BND_test();  //输出最优层       BND_tes(bint)里面调用了output_three_way_data(bint)函数
		cout <<endl<< "************完成************" << endl;
		cin >> choice;
	}
	else
	if (choice == 3)
	{
		int bint;
		cout << "请输入边界域的最大范围：";
		cin >> bint;
		cout << "正在加载训练样本..." << endl;
		loadSample("数据集car.txt");
		cout << "正在进行正规化处理..." << endl;
		vec_normalization();
		projectToSphere();

		cout << "正在训练样本分类..." << endl;
		sample_train(0, 0);

		cout << "输出：正域 边界域 负域 ..." << endl;
		output_three_way_data(bint);
		cout << endl << "************分类完成**********" << endl;

		int np = POS_sample_id.size();
		int nb = BND_sample_id.size();
		int nn = NEG_sample_id.size();
		int n = samples.size();

		cout << "正域样本数 = " << np << endl;
		cout << "边界域样本数 = " << nb << endl;
		cout << "负域样本数 = " << nn << endl;
		cout << "正域+负域+边界域 = " << np + nb + nn << endl;
		if ((np + nb + nn) == n)
		{
			cout << "正域+负域+边界域 = 样本总数 = " << n << endl;
		}
		else
			cout << "分类出现错误！" << endl;

		cout << "----------输出最优层----------" << endl;
		BND_test();

		cout << "对总体样本进行测试..." << endl;
		samples_test();//用覆盖先对确定的样本进行划分 ，再用选取的最优层，对边界域样本进行划分
		
		cout << endl << "************总体样本测试完成************" << endl<<endl;
		
		cout << "继续进行 十折交叉测试：" << endl;
		nfoldCrossTest(0, 0);
		cin >> choice;

		/*cout << "是否继续进行 十折交叉测试？（输入1表示进行，否则结束）：";
		int select;
		cin >> select;
		if (select == 1)
		{
			nfoldCrossTest(0, 0);
			cin >> choice;
		}
		else
			cin >> choice;*/
		
	}
	//if (choice == 4)
	//{
	//	int bint;
	//	cout << "请输入边界域的最大范围：";
	//	cin >> bint;
	//	cout << "正在加载训练样本..." << endl;
	//	loadSample("数据集spambase.data");
	//	cout << "正在进行正规化处理..." << endl;
	//	vec_normalization();
	//	projectToSphere();
	//	sample_train(0, 0);

	//	nfoldCrossTest(0, 0);
	//	//print_nfold_result();
	//	cout <<endl<< "************完成************" << endl;
	//	cin >> choice;
	//}
	//int bint;
	//cout << "请输入边界域的最大范围：";
	//cin >> bint;

	//cout << "正在加载训练样本..." << endl;
	//loadSample("chess.txt"); //加载训练样本
	//cout << "正在进行正规化处理..." << endl;
	//vec_normalization();
	//cout << "正在将样本投射到球面..." << endl;
	//projectToSphere();//样本投射到球面

	////nfoldCrossTest(0, 0);
	////cout << "正在训练样本分类..." << endl;
	//
	////sample_train(0, 0);
	//
	//// print_nfold_result();
	////cout << "正在输出结果..." << endl;
	//// 输出三类数据
	////output_three_way_data(bint);
	//cout << "*************************" << endl << endl;

	////samples_test();
	////BND_test();
	//nfoldCrossTest(0, 0);
	//print_nfold_result();
	//cin >> bint;


	return 0;
}

/*********************************************************************
*	主函数
**********************************************************************/
/*int main()
{
printf("1. 训练\n");
printf("2. 测试\n");
printf("3. 十交叉验证\n");
printf("Choose:\n");
int choice=0 ;
cin>>choice ;
getchar() ;

if(choice==1)
{
	double a,b;
	cout<<"输入a="<<endl;
	cin>>a;
	loadSample("iris.txt"); //加载训练样本
	projectToSphere() ;//样本投射到球面
	ExpResult* rs = new ExpResult ;//实验数据空间准备
	memset(rs,0,sizeof(ExpResult)) ;
	Result.push_back(rs) ;
	sample_train(a,b) ;
	save_model_to_file("model.txt") ;	//保存训练结果数据
	print_train_result() ; //显示训练结果
}

else 
	if(choice==2)
	{
		loadSample("iris.txt"); //加载训练样本
		projectToSphere() ;//样本投射到球面
		ExpResult *rs = new ExpResult ;
		memset(rs,0,sizeof(ExpResult)) ;
		Result.push_back(rs) ;		//实验数据空间准备
		load_model_from_file("model.txt") ;	//加载模型数据
		sample_test() ;
		print_test_result() ; //显示测试结果
	}
	else 
		if(choice==3)
		{
			double a,b;
			cout<<"a=";
			cin>>a;
			cout<<"b=";
			cin>>b;
			loadSample("chess.txt"); //加载训练样本
			vec_normalization();
			projectToSphere() ;//样本投射到球面

			nfoldCrossTest(a,b) ;
			print_nfold_result() ;
		}
		else
		{
			cout<<"错误!"<<endl ;
		}

getchar() ;
return 0 ;
}

*/

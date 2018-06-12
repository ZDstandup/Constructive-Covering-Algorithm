/****************************
*���򸲸��㷨��ʾ����
*****************************/

#include <iostream>
#include <algorithm> //�㷨��
#include <fstream>   // ���ļ������ĺ�����
#include <string>
#include <vector>
#include <cmath>
#include <cstring> 
#include <cstdlib>
#include <cstdio> 
#include <ctime>   //�����ں�ʱ��ת�����ַ���
#include <cfloat>
#include<math.h>
#include<iomanip>
using namespace std;

typedef struct Sample
{	//�������ݽṹ
	int dim;			//����ά��
	vector<double>	x;	//��������
	int				y;	//�������
}Sample;

typedef struct SampleTag
{//�������
	int id; 			//�������
	int covered;		//�����ŵ������Ƿ񱻸���, 1 ��ʾ�����ǣ�0��ʾû��û���� 
}SampleTag;

typedef struct SampleSubSet
{	//�����Ӽ���ţ�����ʾͬһ�����������ı��			
	int y; 					//�Ӽ������ 
	vector<SampleTag*> v;		//id����
}SampleSubSet;

typedef struct Cover
{	//���ǵ����ݽṹ
	int *sid;			// �������
	int seq;			//�������
	int cls;			//���ǵ����
	double r;		 	//���ǰ뾶
	double* center;	//���ǵ�Բ��
	double ybs;			//���ǵ������� 
}Cover;

typedef struct ExpResult
{	//ʵ�������ݽṹ
	int a0, a1, a2, a3, a4, a5;
	int trNum;			//ѵ��������
	int teNum;			//����������
	int covNum;		//������
	double trTime;		//ѵ����ʱ
	double teTime;		//���Ժ�ʱ
	int refuse;		//���Ծ�ʶ������
	int guess_corr;	//��ʶ�����и���Բ������²���ȷ��������
	int correct;		//��ȷʶ��������
	double corr_rate;	//��ȷʶ���ʣ�(corr_rate + guess_corr)/teNum
	int C0N, C1N;//C0,C1�า�ǵ���Ŀ
}ExpResult;



//########################################################################## 
//ȫ�ֱ��� 
int G_lev = 0; //�������Ų����
int Cn = 0;
vector<int> BND_sample_id;      //���� �߽��� ���� id
vector<int> POS_sample_id;      //���� ����   ���� id
vector<int> NEG_sample_id;      //���� ����   ���� id
vector<Sample*> samples;		//��������
vector<SampleSubSet*> I;		//�����Ӽ����
vector<Cover*> C;				//����
vector<Cover*> C1;
vector<Cover*> C2, C3, C4;
vector<ExpResult*> Result;		//ʵ����	
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
*	���������������ȡ�ĸ������������ַ���str���ַ�spliter�ֿ���
*	������vector<string>��
**********************************************************************/
void split(string& str, char spliter, vector<string>& vec)
{
	string::iterator iter1, iter2;  //������(iterator)��һ���������Ա���������Ԫ�أ���ʵ��Ԫ�ر�������������
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

	//�������һ���ַ���
	if (iter1 != str.end() && *iter1 != spliter)
	{
		vec.push_back(string(iter1, iter2));
	}
}

/*********************************************************************
*	���������������ȡ�ĸ�����������һ������UCI��ʽ��һ�������ַ���s
*	ת�����Զ����sample��ʽ
**********************************************************************/
void dealSample_uci_format(string& s, Sample& sample)
{
	char spliter[] = { ',', ' ', '\t', ';', ':' };   //���������ָ��� 
	int k = -1;
	for (int i = 0; i < 5; ++i)
	{//�Զ�Ѱ�������ָ���
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
		printf("������ʽ�д�\n");
		exit(0);
	}

	//�ָ�����
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
*	�����������������������ļ��� sample_set_file��������ֵ��ȡ��
*	ȫ�ֱ���samples��
**********************************************************************/
void loadSample(const string& sample_set_file)
{
	if (sample_set_file == "")
	{
		cout << "�����ļ�����Ϊ��!" << endl;
		getchar();
		exit(0);
	}
	//�������ļ�
	fstream ifs(sample_set_file.c_str(), ios::in); //.c_str()ָ���ļ�sample_set_file��������� 
	if (!ifs)
	{
		cout << "�������ļ�ʧ��!" << endl;
		getchar();
		exit(0);
	}

	string line = "";
	while (getline(ifs, line)) // getline(cin,s)�������ã���cin�ַ�����ֵ��s
	{
		Sample* sample = new Sample;
		dealSample_uci_format(line, *sample); //����uci��ʽ������
		samples.push_back(sample);
	}
	ifs.close();
}

/*********************************************************************
*	��������ȽϺ������������������������������Ŀ��ֵ
*	���������� ��a��b��������ָ��
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
*	����������������������Ŀ��ֵ���������򣬲�������������
*	ȫ�ֱ���I�У�I[k]��ʾ��k�������������ŵļ���
**********************************************************************/
void sortSample()
{
	//�������������򣬽���ͬ��������һ��
	sort(samples.begin(), samples.end(), cmp); //�����Ѿ���������������������
	SampleSubSet* subset = new SampleSubSet;   //����һ���µ�subset�Ӽ�
	int y = samples[0]->y;
	subset->y = y;                             //����һ����������� ��Ϊ�Ӽ�subset�����
	SampleTag* tag = new SampleTag;    //����һ���µ��������tag,���Ҹ��丳��һ���µ�id �Լ�covered
	tag->id = 0;
	tag->covered = 0;
	subset->v.push_back(tag);          //���½����������tag��ӵ��Ӽ�subset�У���Ϊ��ʼ
	for (size_t i = 1; i < samples.size(); ++i) // ���ڴ�samples��һ��������ʼ�����Ƿ񱻸���Ϊ���б�׼�����з���
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
			I.push_back(subset);     //I�д洢�����Ӽ����࣬Ҳ���������������Ӽ���������һ�� + û�б�����һ��
			subset = new SampleSubSet; //����I�в�ֹ2�࣬Ӧ��˵������n�ֻ࣬���кܶ�������ͬ��Ҫ�ǰ����Ƿ񱻸���������Ļ����ܹ�����2��
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
*	��һ������������ģ��    sqrt()������ƽ��������
**********************************************************************/
double vlen(Sample& s)
{//��һ������ģ��
	double a = 0.0;
	int i;
	for (i = 0; i < s.dim; ++i)
	{
		a += (s.x[i] * s.x[i]);
	}
	return sqrt(a);
}

/**************************************************************************
*	���ԵĹ�һ��
*	���� vector<Sample*>& samples
*	��� ��һ����samples
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
	{// ��������ÿά�����ֵ����Сֵ
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

		for (i = 2; i < 2 * (samples.size() / 2); i += 2)   //�˴���  2*��a/2)�������ǣ�������i+2�������Խ������
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
		if (samples.size() % 2 != 0)   //������������Ƿ�Ϊ������Ҫ��Ϊ�����Ļ������ú����һ��û�д�����������бȽ�
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

	//��һ��
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
*	������Ͷ�䵽���� �������й�һ������
**********************************************************************/
void projectToSphere()
{
	size_t i;
	double d = vlen(*samples[0]);
	for (i = 1; i < samples.size(); ++i)
	{//�����ģ������������d
		double t = vlen(*samples[i]);
		if (d < t)
			d = t;
	}
	for (i = 0; i < samples.size(); ++i)
	{//����һά��ӳ�䵽����
		double x = vlen(*samples[i]);
		double t = d*d - x*x;
		t = sqrt(t);
		samples[i]->x.push_back(t); //����һά
		samples[i]->dim += 1;
	}

	for (i = 0; i < samples.size(); ++i)
	{//��һ������
		double x = vlen(*samples[i]);
		for (int k = 0; k < samples[i]->dim; ++k)
		{
			samples[i]->x[k] /= x;
		}
	}
}

/*********************************************************************
*	����������֮���ŷ�Ͼ���
**********************************************************************/
double inner_product(Sample& s1, Sample& s2)
{//��ŷ�Ͼ��� 
	if (s1.dim != s2.dim)
	{
		cout << "s1,s2ά�Ȳ���ȣ�" << endl;
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
	����������֮���ŷ�Ͼ���
    �ڻ��ֺõ����Բ��Ͻ��м��� ����֮���ŷʽ����
**********************************************************************/
double inner_product_subset(Sample& s1, double* x, int y,int c)
{   //�������Ӽ��ķ�Χ��
	//������������֮���ŷʽ����
	
	double a = 0.0;
	if (c == 0)
	{
		int t = (sizeof(SubSet_POS[y]) / sizeof(int));  // ����SutSet�����ж�����
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
		int t = (sizeof(SubSet_NEG[y]) / sizeof(int));  // ����SutSet�����ж�����
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
*		������s1������x��ŷ�Ͼ���
**********************************************************************/
double inner_product_x(Sample& s1, double* x)
{//��ŷ�Ͼ���
	
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
*		������x1������x2��ŷ�Ͼ���
**********************************************************************/
double inner_product(double* x1, double* x2, int dim)
{//��ŷ�Ͼ��� 
	double a = 0.0;
	for (int i = 0; i < dim; ++i)
	{
		a += (x1[i] - x2[i])*(x1[i] - x2[i]);
	}
	return sqrt(a);
}
/*********************************************************************
*	���t��������s�����������ֵ����Сŷ�Ͼ��룩��s������������Ӽ�id
**********************************************************************/
double find_diffMin_d(int s, int t)                                                                                                    
{
	int t1 = (t + 1) % I.size();  //t1��ʾ��t+1)�������Ӽ���SampleSubSet)
	int k = I[t1]->v[0]->id;  //k��ʾ���ǣ���t+1)��->����->���   �˴������Ǵ�0��ʼ��
	int a = I[t]->v[s]->id;   //a��ʾ���� t��->s����->���

	double diffMin_d = inner_product(*samples[a], *samples[k]);  //��t����s�������������е����������ѡȡ�������˶��ߵ�ŷ�Ͼ��룬��Ϊ�����diffMin_dֵ
	for (size_t i = 0; i < I.size(); ++i)
	{
		if ((int)i != t)
		{
			for (size_t j = 0; j < I[i]->v.size(); ++j)
			{
				k = I[i]->v[j]->id; //ʵ��id
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
*	���t��������s���ͬ����Զ���루���ŷ�Ͼ��룩��s������������Ӽ�id
**********************************************************************/
double find_SameMax_d(int s, int t, double diffMin_d)   //������ͬ���������룬һ��Ҫ�������е���С����С�ŷ���
{
	int a = I[t]->v[s]->id; //ԭ��ʵ��id��
	double sameMax_d = 0;//�������ڻ�����֮��Ĳ���ڴ� 
	for (size_t j = 0; j < I[t]->v.size(); ++j)
	{
		int k = I[t]->v[j]->id; //ʵ��id��
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
*	�������ĵ�s����������t�����ǰ뾶d�������ǣ������ظ��ǵĵ�����
*	����ִ�й������޸������Ӽ��������ŵ�covered���ԣ��Ա���������
**********************************************************************/
int cover_sample(int s, int t, double d, Cover *c)
{
	int cov_num = 0;
	int a = I[t]->v[s]->id;
	c->sid = new int[I[t]->v.size()];
	//int j = 0;
	for (size_t i = 0; i < I[t]->v.size(); ++i)
	{
		int k = I[t]->v[i]->id; //ʵ��id
		if (k == a)
		{//Բ�ı���
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
	return cov_num; //���ظ��ǵ�������
}

/*********************************************************************
*	ѵ����ɺ󽫸��ǵĽ��������ļ����ļ���Ϊ model_file
**********************************************************************/
void save_model_to_file(const string& model_file)
{
	if (model_file == "")
	{
		cout << "ģ���ļ�����Ϊ��!" << endl;
		getchar();
		exit(0);
	}
	ofstream ofs(model_file.c_str(), ios::out);
	if (!ofs == NULL)
	{
		cout << "����ģ���ļ�ʧ��!" << endl;
		getchar();
		exit(0);
	}
	int dim = samples[0]->dim;
	//��ʼд�ļ�
	ofs << C.size() << endl;					//д�븲����
	ofs << dim << endl;						//д�븲��ά�ȣ�������Ͷ����ά��
	ofs.precision(10);						//�����ԣ�����̫С�������ģ�ͻ�����ͬ�����ٲ���ʱ���о�ʶ����
	for (size_t i = 0; i < C.size(); ++i)
	{
		ofs << C[i]->seq << ","			//д�븲�����
			<< C[i]->cls << ","			//д�븲�����
			<< C[i]->r << ",";			//д�븲�ǰ뾶
		ofs << "\t\t\t\t";
		for (int k = 0; k < dim; ++k)		//ѭ��д�븲��Բ��
		{
			ofs << fixed << C[i]->center[k] << ", ";
		}
		ofs << endl;
	}
	ofs.close();
}

/*********************************************************************
*	���ļ��ж�ȡѵ�����ĸ������ݣ� �ļ���Ϊ model_file
**********************************************************************/
void load_model_from_file(const string& model_file)
{
	if (model_file == "")
	{
		cout << "ģ���ļ�����Ϊ��!" << endl;
		getchar();
		exit(0);
	}
	ifstream ifs(model_file.c_str(), ios::in);
	if (!ifs == NULL)
	{
		cout << "��ģ���ļ�ʧ��!" << endl;
		getchar();
		exit(0);
	}
	int covnum = 0;
	int dim = 0;
	string line;
	//��ȡ��������
	ifs >> covnum;	//��ȡ������Ŀ
	ifs >> dim;		//��ȡ����ά��
	getline(ifs, line);
	for (int i = 0; i < covnum; ++i)
	{
		Cover* c = new Cover();
		c->center = new double[dim];
		vector<string> vec_str;

		getline(ifs, line);
		split(line, ',', vec_str);

		c->seq = atoi(vec_str[0].c_str());	//�������
		c->cls = atoi(vec_str[1].c_str());	//�������
		c->r = atof(vec_str[2].c_str());	//���ǰ뾶

		for (int k = 0; k < dim; ++k)
		{
			c->center[k] = atof(vec_str[k + 3].c_str());
		}
		C.push_back(c);
	}
	ifs.close();
}

//���ݸ����е��������Ը��ǽ������� 
void sortCover(Cover * cover, int N)  // ð�������㷨
{
	int x, y;
	Cover temp;
	for (y = 0; y < N - 1; y++)
	{
		for (x = 1; x<N - y; x++)  //���ﲻ��Ӧ�ö�x����++ô��Ϊɶ���y++��
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
*	����ѵ������
**********************************************************************/
void sample_train(double a, double b)
{//�󸲸�
	clock_t t1;

	t1 = clock();
	sortSample();	//���������������
	int seq = 0;
	int coved_num = 0;
	//srand((unsigned)time(NULL));
	for (size_t t = 0; t < I.size(); ++t)
	{
		vector<int> v;
		for (size_t i = 0; i < I[t]->v.size(); ++i)
			v.push_back(i);  //��һ����ʱ�ļ���v�ϲ���
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
			double d1 = find_diffMin_d(s, t); //����������� �ڻ����,ŷʽ������С d1 ;
			double d2 = find_SameMax_d(s, t, d1); //��ͬ����Զ�� �ڻ���С��ŷʽ������� d2
			//double d = (d1+d2)/2 ;
			double d = d2;


			Cover* c = new Cover;

			//���ݰ뾶�����������Ǳ��
			coved_num = cover_sample(s, t, d, c); //s��Բ�ģ�t�ǵ�t��������d�Ǹ��ǰ뾶

			//����һ������

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
				C1.push_back(c);//����	
			else
				C2.push_back(c);//���� 


			//�޸�δ������
			uncovedn -= coved_num;
			//printf("%d\t��%d��\t�����Ӹ���:%d\t���ǰ뾶��%f\n",seq,t,coved_num,d);
			seq++;
			v.erase(v.begin());
		}


	}//for

	//��C1,C2�������� 
	//sortCover(C1[0],C1.size());
	//sortCover(C2[0],C2.size());

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//��C1���C2����д�������a,bֵ�����ͷš���a,b�� ;//�����γɸ��Ǽ�C3,C4 
	//int i;
	//for(i=0;i<C1.size()*(1-a);i++)
	//C3.push_back(C1[i]);
	//for(i=0;i<C2.size()*(1-b);i++)
	//C4.push_back(C2[i]);
	//
	////����C�� 
	//for(i=0;i<C3.size();i++)
	//C.push_back(C3[i]);
	//for(i=0;i<C4.size();i++)
	//C.push_back(C4[i]); 
	//cout<<"C.size="<<C.size()<<".";
	//
	//	t2 = clock() ;
	//cout<<" T="<<(t2-t1)/CLOCKS_PER_SEC;
	//	//�洢ʵ������
	//	vector<ExpResult*>::reverse_iterator it = Result.rbegin() ;
	//	(*it)->covNum = C.size() ; //������
	//	(*it)->trNum = samples.size() ; //ѵ��������
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
*	�������Ժ���
**********************************************************************/
void sample_test()
{//����
	//int a[6] = { 0 };
	int a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0;
	int correct = 0;		//��ʶ������ȷ������
	int refuse = 0;		//��ʶ������
	int uc = 0;			//������ʶʱ�����븲����� ������������������� ����ȷ������
	int total_BND = 0;
	int total_correct_BND = 0;
	double correct_index = 0.0;
	//clock_t t1, t2;
	//t1 = clock();

	for (size_t i = 0; i < BND_sample_id.size(); ++i)
	{//���ڱ߽����������ԣ��ܷ���ȷ����

		double cnt_nearest = DBL_MAX; //���������
		//double cnt_nearest1=DBL_MAX;
		int k = -1; //��¼���������ĸ��ǵ��±�
		//int k1= -1;
		size_t j = 0;

		while (j<C1.size())
		{//ÿ������
			
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
		{//�˵������ĵĸ��Ƕ��ڰ뾶֮�⣬��Ϊ��ʶ������
			++refuse;
			if (samples[BND_sample_id[i]]->y == C1[k]->cls) //���ݾ��븲���������������ʶ��������
				++uc;
			if (samples[BND_sample_id[i]]->y == 0)
				++a4;//++a[4];           //a4��¼û�б����ǵ�������
			if (samples[BND_sample_id[i]]->y == 1)
				++a5;//++a[5];		     //a5��¼�����ǵ�������

		}
		else
		{//��ʶ����
			if (samples[BND_sample_id[i]]->y == C1[k]->cls)
			{

				++correct;
				if (samples[BND_sample_id[i]]->y == 0)//���� 
					++a0;//++a[0];
				if (samples[BND_sample_id[i]]->y == 1)
					++a1;//++a[1];

			}
			else
			{//  �߽���
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
	cout << "�߽���������������" << total_BND << endl;
	cout << "�߽������������ֵ������У��ұ���ȷ���ֵ�����Ϊ��" << a0 << endl;
	cout << "��ȷ��Ϊ:" << correct_index << endl;


	

	//t2 = clock();
	//cout << " l to l:" << a0 << "," << "l to s:" << a2 << "," << "l to B:" << a4 << "," << "s to l:" << a3 << "," << "s to s:" << a1 << "," << "s to B:" << a5 << endl;

	////�洢ʵ������
	//vector<ExpResult*>::reverse_iterator it = Result.rbegin();
	//(*it)->a0 = a0;
	//(*it)->a1 = a1;
	//(*it)->a2 = a2;
	//(*it)->a3 = a3;
	//(*it)->a4 = a4;
	//(*it)->a5 = a5;
	//(*it)->correct = correct; //������ȷ��
	//(*it)->refuse = refuse; //ѵ��������
	//(*it)->guess_corr = uc;
	//(*it)->teNum = samples.size();
	//(*it)->teTime = (double)(t2 - t1) / CLOCKS_PER_SEC;
	//(*it)->corr_rate = ((float)(correct + uc)) / samples.size();
}




/************************************************************************/
/* ������࣬���࣬�߽�������                                              */
/************************************************************************/
void output_three_way_data(int bint)
{
	
	
	ofstream ofBND;
	ofBND.open("Car�߽���.txt");

	cout << "����������������..." << endl;
	ofstream ofPOS;
	ofPOS.open("Car����.txt");
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

				// �߽���
				if ((*iter)->ybs < bint)
				{
					//���߽����е�������ű����� BND_sample_id������
					BND_sample_id.push_back(sid);

					// ��ȡ�����������й�һ��������
					vector<double>	x1 = samples[sid]->x;
					//ofBND << sid << ",";        //�����ǽ������������� ��һλ�ˣ������������������� 37 ά
					for (int i = 0; i < (x1.size() - 1); i++)
					{
						ofBND << x1[i] << ",";
					}
					ofBND << (*iter)->cls;
					ofBND << endl;
					//(*iter)->ybs = 0;
					//(*iter)->center = 0;
				}
				// ����
				else {
					// ��ȡ�����������й�һ��������
					POS_sample_id.push_back(sid); //���ֵ����������id������ POS_sample_id ��
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

	cout << "�������ɸ�������..." << endl;
	ofstream ofNEG;
	ofNEG.open("Car����.txt");
	for (vector<Cover *>::iterator iter = C2.begin(); iter != (C2.end()); iter++)
	{
		for (int i = 0; i < (*iter)->ybs; i++)
		{
			int sid = (*iter)->sid[i];
			
				// �߽���
				if ((*iter)->ybs < bint)
				{
					//���߽����е�������ű����� BND_sample_id������
					BND_sample_id.push_back(sid);
					// ��ȡ�����������й�һ��������
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
				// ����
				else {
					// ��ȡ�����������й�һ��������
					NEG_sample_id.push_back(sid);   //�������е�����id ������  NEG_sample_id ��
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

	cout << "�������ɱ߽�������..." << endl;

	//for (size_t i = 0; i < I.size() - 1; ++i){ 

	//	vector<SampleTag*> st = I[i]->v;

	//	for (vector<SampleTag*>::iterator iter = st.begin(); iter != (st.end() - 1); iter++)
	//	{
	//		// �߽�������
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
 �߽������� ���������
 �� BND_test()��ֱ�ӵ���
*/
int sample_test_all(Sample& s)
{
	size_t j1 = 0;
	size_t j2 = 0;
	double d1 = 0.0;  //��¼�߽��������븲������֮��ľ��룬���ұ�����С����
	double d2 = 0.0;
	int m1 = -1;      //��¼���������ĸ����±�
	int m2 = -1;
	int y = -1;        //���������ֺ�����



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

int sample_test_sub(Sample& s,int x,int c)   //x��ʾ����ĵڼ���   c��ʾ������������㣬���Ǹ��������
{//x��ʾSub+(Sub-)�ĵڼ��㣬c��ʾ������Sub+  or ����Sub-
	size_t j1 = 0;
	size_t j2 = 0;
	double d1 = 0.0;  //��¼�߽��������븲������֮��ľ��룬���ұ�����С����
	double d2 = 0.0;
	//int m1 = -1;      //��¼���������ĸ����±�
	//int m2 = -1;
	int y = -1;        //���������ֺ�����

	

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
{/*�����Ų㸳ֵ�� ȫ�ֱ���G_lev*/

	/*
	��ÿһ���߽����������� ȫ���Է���+������Է���
	*/
	//output_three_way_data(bint);
	int a[1000] = { 0 };
	for (size_t i = 0; i < BND_sample_id.size(); i++)
	{
		/*  ��ȫ���������£��жϱ߽�������������һ��

		y = 1��0��������࣬��y1,y2���Ƚ�����
		*/
		int y = sample_test_all(*samples[BND_sample_id[i]]);
		//cout << "������� i= " << i << endl;
		//cout << "--ȫ�� " << y << endl;
		//int a1 = 0;

		for (size_t j = 0; j < (sizeof(SubSet_POS) / sizeof(int) / (sizeof(SubSet_POS[0]) / sizeof(int))); j++)
		{
			/*�����򻮷ֵ����������£��жϱ߽�������������һ��*/

			/*for (int i = 0; i < sizeof(SubSet_POS[j]); i++)
			cout << SubSet_POS[j][i]<<",";
			cout << endl << "***********" << endl;*/

			int y1 = sample_test_sub(*samples[BND_sample_id[i]], j, 0);  //j��ʾ����ĵڼ���   0��ʾ��������������

			//cout << "--���� " << y1 << ",";
			//cout << "******" << y1 << "_________" << endl;
			/*
			y1 �����϶�Ϊ1��Ҳ��������ѡ��Ĳ㣬�󶼽��߽�����������Ϊ ����
			*/
			for (size_t k = 0; k < (sizeof(SubSet_NEG) / sizeof(int) / (sizeof(SubSet_NEG[0]) / sizeof(int))); k++)
			{
				/*�ڸ��򻮷ֵ����������£��жϱ߽�������������һ��

				y2�����϶�Ϊ1��Ҳ���Ǹ��򻮷ֵĲ㣬�󶼽��߽�����������Ϊ ����
				*/
				int y2 = sample_test_sub(*samples[BND_sample_id[i]], k, 1);
				//cout << "--����" << y2 << ",";
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

	G_lev = i_max;  //G_lev�������Ų����
	int count_BND = BND_sample_id.size();
	double rate = 0.0;
	rate = (double)max / (double)count_BND;
	int lev_1 = G_lev / ((sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int)));
	int lev_2 = G_lev % ((sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int)));

	cout << endl << "**************************���Ų���Ͻ������***********************" << endl << endl;
	cout << "���Ų������� = " << i_max << ",Sub+ = " << lev_1 << ",Sub- = " << lev_2 << endl << endl;
	cout << "���Ų���� ��ȷ������� = " << max << endl << endl;

	cout << "********************* BND_test() ���н���****************************" << endl << endl;
}

void samples_test()
{

	int nc1 = 0;  //����ͨ�� �����㷨 ��������������Ϊ ��������ȷ ����������
	int nc2 = 0;  //����ͨ�� �����㷨 ��������������Ϊ ��������ȷ ����������

	int nc1_n = 0;//ͨ�������㷨������������Ϊ ���򣬵�����������������ͬ�������ִ���ĸ���
	int nc2_n = 0;//ͨ�������㷨������������Ϊ ���򣬵�����������������ͬ�������ִ���ĸ���



	int err = 0;

	int bnd = 0;//����߽�������������
	int bnd_1 = 0;//����߽������� �������
	int bnd_2 = 0;//����߽������� �������

	int all = 0;
	int all_1 = 0;
	int all_2 = 0;

	int sub_10 = 0;
	int sub_11 = 0;
	int sub_12 = 0;

	int sub_20 = 0; //��������Sub-���Ա߽�����������ȷ����ĸ���
	int sub_21 = 0;//���� �������
	int sub_22 = 0;

	int zuhe = 0;//��������������� ���Ա߽�����������ȷ�������
	int zuhe_1 = 0;
	int zuhe_2 = 0;

	int b_concret = 0;//���������Ų����Ϸ�ʽ���Ա߽������������Ų���� ���з���
	int b_concret_1 = 0;
	int b_concret_2 = 0;

	int great_lev = G_lev;  //�����Ų㣬�Ա߽����������з���
	for (int i = 0; i < samples.size(); i++)
	{
		int flag1 = 0; // ������־ �����Ƿ� ����
		//int flag2 = 0;

		int j1 = 0; //��C1����ѭ������
		int j2 = 0; //��C2����ѭ������
		double d1 = DBL_MAX; //���� ������ C1�������� �� ��С��������
		double d2 = DBL_MAX; //���� ������ C2�������� �� ��С��������

		//int x1 = 0;
		//int x2 = 0;
		while (j1 < C1.size())
		{ // ���������� C1 �����е���С����
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
		���� C1 , C2 ѭ�����ǽ��������ж�����
		�������������̺� δ���������е����� ��һ�����࣬
		������ ѵ���ҵ������Ų� ������ ��һ��������
		*/
		if (flag1 == 0)
		{
			int lev_1 = great_lev / ((sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int)));
			int lev_2 = great_lev % ((sizeof(SubSet_NEG) / sizeof(int)) / (sizeof(SubSet_NEG[0]) / sizeof(int)));

			bnd = bnd + 1; //����߽�����������
			if (samples[i]->y == 0)
			{
				bnd_1 = bnd_1 + 1; //�߽��������� ������������
			}
			if (samples[i]->y == 1)
			{
				bnd_2 = bnd_2 + 1; //�߽��������� ������������
			}

			int y = sample_test_all(*samples[i]); //����ȫ���ԣ����������з���
			int y1 = sample_test_sub(*samples[i], lev_1, 0);//������Sub+���Բ㣬���������з���
			int y2 = sample_test_sub(*samples[i], lev_2, 1);//������Sub-���Բ㣬���������з���

			//cout << "y=" << y << ",";
			if (y == samples[i]->y) // ȫ���Բ���
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
			if (y1 == samples[i]->y) //����Sub+�� �߽����������в���
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
			if (y2 == samples[i]->y)//����Sub-�� �߽����������в���
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
			int o1 = -1;   //�������������
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

	double rate;//������ȷ��
	rate = double(correct) / double(samples.size());

	double rate_all = 0;//ȫ������ȷ��
	rate_all = double(all) / double(bnd);

	double rate_sub1 = 0;// Sub+��ȷ��
	rate_sub1 = double(sub_10) / double(bnd);

	double rate_sub2 = 0;// Sub-��ȷ��
	rate_sub2 = double(sub_20) / double(bnd);

	double rate_zuhe = 0;// �����ȷ��
	rate_zuhe = double(zuhe) / double(bnd);

	cout << endl;
	cout << "+++++++++++������� ��ȷ������� = " << b_concret << ",������� = " << b_concret_1 << ",������� = " << b_concret_2 << endl << "��ȷ�� = " << double(b_concret) / double(bnd) << endl << endl << endl;

	cout << endl << "----------�߽��� = " << bnd << ",���� ������������ = " << bnd_1 << ",������������ = " << bnd_2 << endl << endl;
	cout << "���� ȫ���ԶԱ߽����������з��࣬��ȷ������� = " << all << ",���� ������������ = " << all_1 << ",������������ = " << all_2 << endl;
	cout << "-----------��ȷ�� = " << rate_all << endl << endl;
	cout << "���� ����Sub+���Բ�Ա߽����������з��࣬��ȷ������� = " << sub_10 << ",���� ������������ = " << sub_11 << ",������������ = " << sub_12 << endl;
	cout << "-----------��ȷ�� = " << rate_sub1 << endl << endl;
	cout << "���� ����Sub-���Բ�Ա߽����������з��࣬��ȷ������� = " << sub_20 << ",���� ������������ = " << sub_21 << ",������������ = " << sub_22 << endl;
	cout << "-----------��ȷ�� = " << rate_sub2 << endl << endl;
	cout << "���� ����Sub+������Sub-����� �Ա߽����������з��࣬��ȷ������� = " << zuhe << ",���� ������������ = " << zuhe_1 << ",������������ = " << zuhe_2 << endl;
	cout << "-----------��ȷ�� = " << rate_zuhe << endl << endl;

	cout << endl;
	cout << "----------ͨ�����Ƕ��������л��֣���ȷ���ֵĸ��� = " << nc1 + nc2 << ",���У������������� = " << nc1 << ",������������ = " << nc2 << endl;
	cout << "----------���󻮷ֵĸ��� = " << nc1_n + nc2_n << "--------------" << endl << endl;

	cout << "----------ͨ�����Ų���ϣ���ȷ����ĸ��� =" << zuhe;
	cout << "----------�������ĸ��� = " << err << "---------------" << endl << endl;

	cout << "----------��ȷ��������Ϊ = nc1+nc2+zuhe=" << correct << "------------" << endl << endl;
	cout << "----------�����������Ϊ =nc1_n+nc2_n+err= " << fault << "------------" << endl << endl;

	cout << "----------��������Ϊ = samples.size()=" << samples.size() << "------------" << endl << endl;
	cout << "----------��ȷ�� = " << rate << "-------------" << endl;
}


/*********************************************************************
*	�������ļ�ִ��ʮ������֤
**********************************************************************/
void nfoldCrossTest(double a, double b)
{
	const int N = 10;
	int sampleCnt = samples.size();	//��������
	int unit = sampleCnt / N;		//n fold ÿ��������
	int sampleSelectedCnt = N*unit;	//����ȡ����Ϊ�˼��㹫ƽ��ȥ�����ֺ����ͷ

	vector<Sample*> samplesBak(samples);
	random_shuffle(samplesBak.begin(), samplesBak.end()); //��������˳��
	vector<Sample*>::iterator itstart = samplesBak.begin();
	vector<Sample*>::iterator iend = samplesBak.end();

	iend -= (sampleCnt - sampleSelectedCnt);  // ȥ�����ֺ����ͷ ��Ŀ
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
		Result.push_back(rel);	//ʵ�����ݿռ�׼��
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
*	��ʾ����ѵ�����
**********************************************************************/
void print_train_result()
{
	printf(
		"ѵ�������\n"
		"ѵ������Ϊ��%d����ά�ȣ�%d�����%d\n"
		"�õ���������%d��\n"
		"ѵ����ʱ:%.4fs\n"
		, samples.size(), samples[0]->dim - 1, I.size(), Result[0]->covNum, Result[0]->trTime
		); //samples[0]->dim - 1����Ϊ������Ͷ�䵽����ʱ������һά
}

/*********************************************************************
*	��ʾ�������Խ��
**********************************************************************/
void print_test_result()
{
	printf("���Խ����\n");
	printf("-----------------------------------------------------------------------");
	printf("-----------------------------------------------------------------------\n");
	printf("����������\t��ʶ������\t��ʶ��ȷ��\t��ʶ��ȷ��\t��ʶ������\t��ʶ��ȷ��\t��ʶ��\t��ʶ��ȷ��\t����ȷ��\t��ʱ\n");
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
*	��ʾʮ������֤���
**********************************************************************/
void print_nfold_result()
{
	printf("10������֤:\n");
	printf("-------------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");
	printf("ѵ��������\t������\t����������\t��ʶ������\t��ʶ��ȷ��\t��ʶ��ȷ��\t��ʶ������\t��ʶ��ȷ��\t��ʶ��\t��ʶ��ȷ��\t����ȷ��\t��ʱ\t\t��to��\t\t��to��\t��to�߽�\t��to��\t��to��\t��to�߽�\t\tC0Num\tC1Num\n");
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
	printf("(ƽ��ֵ)\n");
	printf("%d\t\t%d\t%d\t\t%f\t\t%f\t\t%.2f%%\t\t%f\t\t%f\t\t%.2f%%\t%.2f%%\t\t%.2f%%\t\t%.4fs\t\%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
		s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s14, s16, s15, s13, s17, s21, s22);
	printf("-------------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");
}



int main()
{
	//char btnums[4];
	//cout << "��������ֱ߽���ķ�Χ����ʽ��1����1,2����";

	//cin >> btnums;

	cout << "���� 1 ����ʾ����ѵ������������������ ���򣬱߽��򣬸��� ���࣬������txt�ļ����" << endl;
	cout << "���� 2 ����ʾѡȡ��������������Բ㣬��������������Բ�" << endl;
	cout << "���� 3 ����ʾ���� �����������,�Լ� �Ƿ����ʮ�۽������" << endl;
	//cout << "���� 4 ����ʾ����ʮ�۽��� ����" << endl;
	cout << "�����루1,2,3��:";
	
	int choice;
	cin >> choice;
	if (choice == 1)
	{
		int bint;
		cout << "������߽�������Χ��";
		cin >> bint;
		cout << "���ڼ���ѵ������..." << endl;
		loadSample("���ݼ�car.txt");
		cout << "���ڽ������滯����..." << endl;
		vec_normalization();
		projectToSphere();

		cout << "����ѵ����������..." << endl;
		sample_train(0, 0);

		cout << "��������� �߽��� ���� ..." << endl;
		output_three_way_data(bint);
		cout << endl << "************�������**********" << endl;
		
		int np = POS_sample_id.size();
		int nb = BND_sample_id.size();
		int nn = NEG_sample_id.size();
		int n = samples.size();
		
		cout << "���������� = " << np << endl;
		cout << "�߽��������� = " << nb << endl;
		cout << "���������� = " << nn << endl;
		cout << "����+����+�߽��� = " << np + nb + nn << endl;
		if ((np + nb + nn) == n)
		{
			cout << "����+����+�߽��� = �������� = " << n << endl;
		}
		else
			cout << "������ִ���" << endl;
		cin >> choice;
	}
	else if (choice == 2)
	{
		int bint;
		cout << "������߽�������Χ��";
		cin >> bint;
		cout << "���ڼ���ѵ������..." << endl;
		loadSample("���ݼ�car.txt");
		cout << "���ڽ������滯����..." << endl;
		vec_normalization();
		projectToSphere();

		cout << "����ѵ����������..." << endl;
		sample_train(0, 0);

		//cout << "��������� �߽��� ���� ..." << endl;
		output_three_way_data(bint);
		cout << endl << "************�������**********" << endl;

		int np = POS_sample_id.size();
		int nb = BND_sample_id.size();
		int nn = NEG_sample_id.size();
		int n = samples.size();

		cout << "���������� = " << np << endl;
		cout << "�߽��������� = " << nb << endl;
		cout << "���������� = " << nn << endl;
		cout << "����+����+�߽��� = " << np + nb + nn << endl;
		if ((np + nb + nn) == n)
		{
			cout << "����+����+�߽��� = �������� = " << n << endl;
		}
		else
			cout << "������ִ���" << endl;

		cout << "----------������Ų�----------"<<endl;
		BND_test();  //������Ų�       BND_tes(bint)���������output_three_way_data(bint)����
		cout <<endl<< "************���************" << endl;
		cin >> choice;
	}
	else
	if (choice == 3)
	{
		int bint;
		cout << "������߽�������Χ��";
		cin >> bint;
		cout << "���ڼ���ѵ������..." << endl;
		loadSample("���ݼ�car.txt");
		cout << "���ڽ������滯����..." << endl;
		vec_normalization();
		projectToSphere();

		cout << "����ѵ����������..." << endl;
		sample_train(0, 0);

		cout << "��������� �߽��� ���� ..." << endl;
		output_three_way_data(bint);
		cout << endl << "************�������**********" << endl;

		int np = POS_sample_id.size();
		int nb = BND_sample_id.size();
		int nn = NEG_sample_id.size();
		int n = samples.size();

		cout << "���������� = " << np << endl;
		cout << "�߽��������� = " << nb << endl;
		cout << "���������� = " << nn << endl;
		cout << "����+����+�߽��� = " << np + nb + nn << endl;
		if ((np + nb + nn) == n)
		{
			cout << "����+����+�߽��� = �������� = " << n << endl;
		}
		else
			cout << "������ִ���" << endl;

		cout << "----------������Ų�----------" << endl;
		BND_test();

		cout << "�������������в���..." << endl;
		samples_test();//�ø����ȶ�ȷ�����������л��� ������ѡȡ�����Ų㣬�Ա߽����������л���
		
		cout << endl << "************���������������************" << endl<<endl;
		
		cout << "�������� ʮ�۽�����ԣ�" << endl;
		nfoldCrossTest(0, 0);
		cin >> choice;

		/*cout << "�Ƿ�������� ʮ�۽�����ԣ�������1��ʾ���У������������";
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
	//	cout << "������߽�������Χ��";
	//	cin >> bint;
	//	cout << "���ڼ���ѵ������..." << endl;
	//	loadSample("���ݼ�spambase.data");
	//	cout << "���ڽ������滯����..." << endl;
	//	vec_normalization();
	//	projectToSphere();
	//	sample_train(0, 0);

	//	nfoldCrossTest(0, 0);
	//	//print_nfold_result();
	//	cout <<endl<< "************���************" << endl;
	//	cin >> choice;
	//}
	//int bint;
	//cout << "������߽�������Χ��";
	//cin >> bint;

	//cout << "���ڼ���ѵ������..." << endl;
	//loadSample("chess.txt"); //����ѵ������
	//cout << "���ڽ������滯����..." << endl;
	//vec_normalization();
	//cout << "���ڽ�����Ͷ�䵽����..." << endl;
	//projectToSphere();//����Ͷ�䵽����

	////nfoldCrossTest(0, 0);
	////cout << "����ѵ����������..." << endl;
	//
	////sample_train(0, 0);
	//
	//// print_nfold_result();
	////cout << "����������..." << endl;
	//// �����������
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
*	������
**********************************************************************/
/*int main()
{
printf("1. ѵ��\n");
printf("2. ����\n");
printf("3. ʮ������֤\n");
printf("Choose:\n");
int choice=0 ;
cin>>choice ;
getchar() ;

if(choice==1)
{
	double a,b;
	cout<<"����a="<<endl;
	cin>>a;
	loadSample("iris.txt"); //����ѵ������
	projectToSphere() ;//����Ͷ�䵽����
	ExpResult* rs = new ExpResult ;//ʵ�����ݿռ�׼��
	memset(rs,0,sizeof(ExpResult)) ;
	Result.push_back(rs) ;
	sample_train(a,b) ;
	save_model_to_file("model.txt") ;	//����ѵ���������
	print_train_result() ; //��ʾѵ�����
}

else 
	if(choice==2)
	{
		loadSample("iris.txt"); //����ѵ������
		projectToSphere() ;//����Ͷ�䵽����
		ExpResult *rs = new ExpResult ;
		memset(rs,0,sizeof(ExpResult)) ;
		Result.push_back(rs) ;		//ʵ�����ݿռ�׼��
		load_model_from_file("model.txt") ;	//����ģ������
		sample_test() ;
		print_test_result() ; //��ʾ���Խ��
	}
	else 
		if(choice==3)
		{
			double a,b;
			cout<<"a=";
			cin>>a;
			cout<<"b=";
			cin>>b;
			loadSample("chess.txt"); //����ѵ������
			vec_normalization();
			projectToSphere() ;//����Ͷ�䵽����

			nfoldCrossTest(a,b) ;
			print_nfold_result() ;
		}
		else
		{
			cout<<"����!"<<endl ;
		}

getchar() ;
return 0 ;
}

*/